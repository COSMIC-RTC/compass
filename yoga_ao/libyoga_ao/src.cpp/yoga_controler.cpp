#include <yoga_controler.h>
#include <string>


yoga_controler::yoga_controler(yoga_context *context, long nvalid, long nactu, long delay, int device,const char *typec)
{
  this->current_context=context;
  //long *dims_data2 = new long[3]; dims_data2[0] = 2; 
  //long *dims_data3 = new long[4]; dims_data3[0] = 3;

  this->nvalid      = nvalid;
  this->nactu       = nactu;
  this->delay       = delay;
  this->device      = device;
  this->typec       = string(typec);
  this->gain        = 0.0f;

  this->nstreams = 1;//nvalid/10;
  while(nactu % this->nstreams!=0) nstreams--;
  cerr << "controler uses " << nstreams << " streams" << endl;
  streams = new yoga_streams(nstreams);

  long *dims_data2 = new long[3]; dims_data2[0] = 2;
  dims_data2[1]    = 2*nvalid; dims_data2[2] = nactu;
  this->d_imat     = new yoga_obj<float>(context, dims_data2);
  this->h_imat     = new yoga_host_obj<float>(dims_data2,MA_PAGELOCK);

  dims_data2[1] = nactu; dims_data2[2] = 2*nvalid;
  this->d_cmat  = new yoga_obj<float>(context, dims_data2);
  this->h_cmat  = new yoga_host_obj<float>(dims_data2,MA_PAGELOCK);

  dims_data2[1]   = 2*nvalid; dims_data2[2] = 2*nvalid;
  this->h_mes2mod = new yoga_host_obj<float>(dims_data2,MA_PAGELOCK);
  this->d_mes2mod = new yoga_obj<float>(context, dims_data2);

  this->d_cenbuff = 0L;

  if (delay > 0) {
    dims_data2[2] = delay+1;
    this->d_cenbuff = new yoga_obj<float>(context, dims_data2);
  }

  dims_data2[1]   = nactu; dims_data2[2] = nactu;
  this->h_mod2act = new yoga_host_obj<float>(dims_data2,MA_PAGELOCK);
  this->d_mod2act = new yoga_obj<float>(context, dims_data2);

  long *dims_data1  = new long[2]; dims_data1[0] = 1;
  dims_data1[1]     = 2*nvalid < nactu ? 2*nvalid : nactu;
  this->d_eigenvals = new yoga_obj<float>(context, dims_data1);
  this->h_eigenvals = new yoga_host_obj<float>(dims_data1,MA_PAGELOCK);

  dims_data2[1]     = dims_data1[1]; dims_data2[2] = dims_data1[1];
  this->h_mateigens = new yoga_host_obj<float>(dims_data2,MA_PAGELOCK);
  this->d_mateigens = new yoga_obj<float>(context, dims_data2);

  dims_data1[1]     = 2*nvalid;
  this->d_centroids = new yoga_obj<float>(context,dims_data1);
  dims_data1[1]     = nactu;
  this->d_com       = new yoga_obj<float>(context,dims_data1);
  this->d_err       = new yoga_obj<float>(context,dims_data1);
  this->d_gain      = new yoga_obj<float>(context,dims_data1);

  delete[] dims_data1;
  delete[] dims_data2;

  //yoga_checkCublasStatus(cublasCreate(&(this->cublas_handle)));
}


yoga_controler::~yoga_controler()
{
  current_context->set_activeDevice(device);
  delete this->streams;
  delete this->d_imat;
  delete this->d_cmat;
  delete this->d_gain;

  delete this->h_imat;
  delete this->h_cmat;
  delete this->d_eigenvals;
  delete this->h_eigenvals;
  delete this->h_mateigens;
  delete this->d_mateigens;
  delete this->h_mes2mod;
  delete this->d_mes2mod;
  delete this->h_mod2act;
  delete this->d_mod2act;

  delete this->d_centroids;
  if (this->delay > 0) delete this->d_cenbuff;
  delete this->d_com;
  delete this->d_err;

  //yoga_checkCublasStatus(cublasDestroy(this->cublas_handle));

  //delete this->current_context;
}

int yoga_controler::svdec_imat()
{
  /**/
  /////// HOST VERSION /////
  // transfering to host memory for svd
  this->h_imat->cpy_obj(this->d_imat, cudaMemcpyDeviceToHost);

  // doing svd
  yoga_cula_svd(this->h_imat,this->h_eigenvals, this->h_mod2act,this->h_mes2mod);

  this->h_mod2act->cpy_obj(this->d_mod2act, cudaMemcpyHostToDevice);
  this->h_mes2mod->cpy_obj(this->d_mes2mod, cudaMemcpyHostToDevice);

  /*
/////// DEV VERSION /////
// transfering to host memory for svd
this->h_imat->cpy_obj(this->d_imat, cudaMemcpyDeviceToHost);

// doing svd
yoga_cula_svd(this->d_imat,this->d_eigenvals, this->d_mod2act,this->d_mes2mod);
this->h_eigenvals->cpy_obj(this->d_eigenvals, cudaMemcpyDeviceToHost);
  */
  return EXIT_SUCCESS;
}

int yoga_controler::set_gain(float gain)
{
  this->gain = gain;
  return EXIT_SUCCESS;
}

int yoga_controler::load_mgain(float *mgain)
{
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}


int yoga_controler::set_delay(int delay)
{
  this->delay = delay;
  return EXIT_SUCCESS;
}


int yoga_controler::build_cmat(int nfilt, bool filt_tt)
{
  int nb_elem = this->h_eigenvals->getNbElem();

  memset(h_mateigens->getData(),0,sizeof(float)*nb_elem);

  // filtering modes
  for (int cc=0;cc<nb_elem-nfilt;cc++) {
    if ((cc < 2) && (filt_tt)) h_mateigens->getData()[cc+cc*nb_elem] = 0.0f;
    if (fabs(h_eigenvals->getData()[cc]) > 1.e-6)
      h_mateigens->getData()[cc+cc*nb_elem] = 1.0f/h_eigenvals->getData()[cc];
  }

  this->h_mateigens->cpy_obj(this->d_mateigens, cudaMemcpyHostToDevice);

  // computing cmat
  long *dims_data2 = new long[3]; dims_data2[0] = 2;
  dims_data2[1]   = nactu; dims_data2[2] = nactu;
  this->d_tmp     = new yoga_obj<float>(current_context,dims_data2);
  
  this->d_tmp->gemm('t','n',1.0f,this->d_mod2act,this->d_mod2act->getDims(1),this->d_mateigens,this->d_mateigens->getDims(1),0.0f,this->d_tmp->getDims(1));
  
  this->d_cmat->gemm('n','t',1.0f,this->d_tmp,this->d_tmp->getDims(1),this->d_mes2mod,this->d_mes2mod->getDims(1),0.0f,this->d_cmat->getDims(1));
  
  delete this->d_tmp;
  delete[] dims_data2;
  
  return EXIT_SUCCESS;
}

int yoga_controler::build_cmat(int nfilt)
{
  return this->build_cmat(nfilt,false);
}

int yoga_controler::frame_delay()
{
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  if (delay > 0) {
    for (int cc=0;cc<delay;cc++) 
      shift_buf(&((this->d_cenbuff->getData())[(cc)*2*this->nvalid]), 1, 2*this->nvalid ,this->device);

    cutilSafeCall( cudaMemcpy(&(this->d_cenbuff->getData()[(delay)*2*this->nvalid]),
			      this->d_centroids->getData(),sizeof(float)*2*this->nvalid,
			      cudaMemcpyDeviceToDevice));

    cutilSafeCall( cudaMemcpy(this->d_centroids->getData(),this->d_cenbuff->getData(),
			      sizeof(float)*2*this->nvalid,cudaMemcpyDeviceToDevice));
  }


  return EXIT_SUCCESS;
}

int yoga_controler::comp_com()
{
  if (this->nstreams > 1) {
    int nstreams = this->nstreams;
    float alpha = -1.0f;
    float beta  = 0.0f;
    //cout << this->d_cmat->getDims(1) << " x " << this->d_cmat->getDims(2) << endl;
    //cout << this->d_centroids->getNbElem()<< endl;
    //cout << this->d_err->getNbElem()<< endl;
    for (int i=0;i<nstreams;i++) {
      int istart1 = i * this->d_cmat->getDims(2) * this->d_cmat->getDims(1) / nstreams;
      int istart2 = i * this->d_cmat->getDims(1) / nstreams;

      //cout << istart1 << " " << istart2 << endl;

      cublasSetStream(cublas_handle, this->streams->get_stream(i));

      cublasOperation_t trans = yoga_char2cublasOperation('n');

      yoga_checkCublasStatus(cublasSgemv(cublas_handle,
					 trans,
					 this->d_cmat->getDims(1)/nstreams,
					 this->d_cmat->getDims(2),
					 &alpha,
					 &((this->d_cmat->getData())[istart1]),
					 this->d_cmat->getDims(1)/nstreams,
					 this->d_centroids->getData(),
					 1,
					 &beta,
					 &((this->d_err->getData())[istart2]),
					 1));
    }


    mult_int(this->d_com->getData(),this->d_err->getData(),this->d_gain->getData(),this->gain,this->nactu,this->device,this->streams);

    this->streams->wait_all_streams();

  } else {
    // compute error
    this->d_err->gemv('n',-1.0f,this->d_cmat,this->d_cmat->getDims(1),this->d_centroids,1,0.0f,1);
    
    // apply modal gain & loop gain
    mult_int(this->d_com->getData(),this->d_err->getData(),this->d_gain->getData(),this->gain,this->nactu,this->device);
  }

  return EXIT_SUCCESS;
}

