#include <sutra_controller_mv.h>
#include <string>

sutra_controller_mv::sutra_controller_mv(carma_context *context, long nvalid,
    long nactu, long delay) :
    sutra_controller(context, nvalid * 2, nactu) {

  this->nvalid = nvalid;
  this->nslope = 2 * nvalid;
  this->nactu = nactu;
  this->delay = delay;
  this->device = device;
  this->gain = 0.0f;

  this->nstreams = 1; //nvalid/10;
  while (nactu % this->nstreams != 0)
    nstreams--;
  cerr << "controller uses " << nstreams << " streams" << endl;
  streams = new carma_streams(nstreams);
  long dims_data2[3];
  dims_data2[0] = 2;
  dims_data2[1] = nslope;
  dims_data2[2] = nactu;
  this->d_imat = new carma_obj<float>(context, dims_data2);
  dims_data2[1] = nactu;
  dims_data2[2] = nslope;
  this->d_cmat = new carma_obj<float>(context, dims_data2);
  dims_data2[1] = dims_data2[2] = nactu;
  d_U = new carma_obj<float>(current_context, dims_data2);
  this->d_cenbuff = 0L;

  if (delay > 0) {
    dims_data2[1] = 2 * nvalid;
    dims_data2[2] = delay + 1;
    this->d_cenbuff = new carma_obj<float>(context, dims_data2);
  }

  long dims_data1[2];
  dims_data1[0] = 1;
  dims_data1[1] = nslope < nactu ? nslope : nactu;
  this->d_eigenvals = new carma_obj<float>(context, dims_data1);
  this->h_eigenvals = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);

  dims_data1[1] = nslope;
  this->d_centroids = new carma_obj<float>(context, dims_data1);
  dims_data1[1] = nactu;
  this->d_com = new carma_obj<float>(context, dims_data1);
  this->d_com1 = new carma_obj<float>(context, dims_data1);
  this->d_com2 = new carma_obj<float>(context, dims_data1);
  this->d_err = new carma_obj<float>(context, dims_data1);
  this->d_gain = new carma_obj<float>(context, dims_data1);

  // Florian features
  dims_data2[1] = nactu;
  dims_data2[2] = nactu;
  this->d_covmat = new carma_obj<float>(context, dims_data2);
  dims_data2[1] = 2 * nvalid;
  dims_data2[2] = nactu;
  this->d_KLbasis = new carma_obj<float>(context, dims_data2);

  cublas_handle = current_context->get_cublasHandle();
  //carma_checkCublasStatus(cublasCreate(&(this->cublas_handle)));
}

sutra_controller_mv::~sutra_controller_mv() {
  current_context->set_activeDevice(device);
  delete this->d_U;

  delete this->streams;
  delete this->d_imat;
  delete this->d_cmat;
  delete this->d_gain;

  delete this->d_eigenvals;
  delete this->h_eigenvals;

  if (this->delay > 0)
    delete this->d_cenbuff;
  delete this->d_com1;
  delete this->d_com2;
  delete this->d_err;
// Florian features
  delete this->d_covmat;
  delete this->d_KLbasis;

  //carma_checkCublasStatus(cublasDestroy(this->cublas_handle));

  //delete this->current_context;
}

string sutra_controller_mv::get_type() {
  return "mv";
}

int sutra_controller_mv::set_gain(float gain) {
  this->gain = gain;
  return EXIT_SUCCESS;
}

int sutra_controller_mv::load_mgain(float *mgain) {
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

// Florian features
int sutra_controller_mv::load_covmat(float *covmat) {
  this->d_covmat->host2device(covmat);
  return EXIT_SUCCESS;
}
int sutra_controller_mv::load_klbasis(float *klbasis) {
  this->d_KLbasis->host2device(klbasis);
  return EXIT_SUCCESS;
}
int sutra_controller_mv::set_delay(int delay) {
  this->delay = delay;
  return EXIT_SUCCESS;
}

// Florian features
int sutra_controller_mv::build_cmat(const char *dmtype) {
  float one = 1.;
  float zero = 0.;

  //if (strcmp(dmtype, "kl") == 0){
  fprintf(stderr, "[%s@%d] here ! \n", __FILE__, __LINE__);
  carma_obj<float> *d_tmp;
  carma_obj<float> *d_tmp2;
  // carma_obj<float> *d_tmp3;
  // carma_obj<float> *d_tmp4;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = nactu;
  dims_data2[2] = nactu;
  d_tmp = new carma_obj<float>(current_context, dims_data2);
  d_tmp2 = new carma_obj<float>(current_context, dims_data2);
  //dims_data2[1] = nactu;
  //dims_data2[2] = nslope;
  //d_tmp3 = new carma_obj<float>(current_context,dims_data2);
  //d_tmp4 = new carma_obj<float>(current_context,dims_data2);

  carma_gemm(cublas_handle, 't', 'n', nactu, nactu, nslope, one,
      d_imat->getData(), nslope, d_imat->getData(), nslope, zero,
      d_tmp->getData(), nactu);
  //add_md(d_tmp->getData(),d_tmp->getData(),d_covmat->getData(),nactu,this->device);

  carma_geam(cublas_handle, 'n', 'n', nactu, nactu, one, d_tmp->getData(),
      nactu, one, d_covmat->getData(), nactu, d_tmp2->getData(), nactu);
  carma_potri(d_tmp2);
  carma_gemm(cublas_handle, 'n', 't', nactu, nslope, nactu, one,
      d_tmp2->getData(), nactu, d_imat->getData(), nslope, zero,
      d_cmat->getData(), nactu);

  delete d_tmp;
  delete d_tmp2;
  //}
  /*
   else{

   //fprintf(stderr,"[%s@%d] here ! \n",__FILE__,__LINE__);
   carma_obj<float> *d_tmp;
   carma_obj<float> *d_tmp2;
   carma_obj<float> *d_tmp3;
   carma_obj<float> *d_tmp4;
   long *dims_data2 = new long[3];
   dims_data2[0] = 2;
   dims_data2[1] = nactu;
   dims_data2[2] = nactu;
   d_tmp = new carma_obj<float>(current_context, dims_data2);
   d_tmp2 = new carma_obj<float>(current_context,dims_data2);
   dims_data2[1] = nactu;
   dims_data2[2] = nslope;
   d_tmp3 = new carma_obj<float>(current_context,dims_data2);
   d_tmp4 = new carma_obj<float>(current_context,dims_data2);

   carma_gemm(cublas_handle,'t','n',nactu,nactu,nslope,one,d_KLbasis->getData(),nslope,d_KLbasis->getData(),nslope,zero,d_tmp->getData(),nactu);
   //carma_geam(cublas_handle,'n','n',nactu,nactu,one,d_tmp->getData(),nactu,one,d_covmat->getData(),nactu,d_tmp2->getData(),nactu);
   add_md(d_tmp->getData(),d_tmp->getData(),d_covmat->getData(),nactu,this->device);
   carma_potri(d_tmp);
   carma_gemm(cublas_handle,'n','t',nactu,nslope,nactu,one,d_tmp->getData(),nactu,d_KLbasis->getData(),nslope,zero,d_tmp4->getData(),nactu);

   carma_gemm(cublas_handle,'t','n',nactu,nactu,nslope,one,d_imat->getData(),nslope,d_imat->getData(),nslope,zero,d_tmp->getData(),nactu);
   carma_potri(d_tmp);
   carma_gemm(cublas_handle,'n','t',nactu,nslope,nactu,one,d_tmp->getData(),nactu,d_imat->getData(),nslope,zero,d_tmp3->getData(),nactu);
   carma_gemm(cublas_handle,'n','n',nactu,nactu,nslope,one,d_tmp3->getData(),nactu,d_KLbasis->getData(),nslope,zero,d_tmp2->getData(),nactu);

   carma_gemm(cublas_handle,'n','n',nactu,nslope,nactu,one,d_tmp2->getData(),nactu,d_tmp4->getData(),nactu,zero,d_cmat->getData(),nactu);

   delete d_tmp;
   delete d_tmp2;
   delete d_tmp3;
   delete d_tmp4;
   }
   */

  return EXIT_SUCCESS;
}

int sutra_controller_mv::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  if (delay > 0) {
    for (int cc = 0; cc < delay; cc++)
      shift_buf(&((this->d_cenbuff->getData())[cc * this->nslope]), 1,
          this->nslope, this->device);

    cutilSafeCall(
        cudaMemcpy(&(this->d_cenbuff->getData()[delay * this->nslope]),
            this->d_centroids->getData(), sizeof(float) * this->nslope,
            cudaMemcpyDeviceToDevice));

    cutilSafeCall(
        cudaMemcpy(this->d_centroids->getData(), this->d_cenbuff->getData(),
            sizeof(float) * this->nslope, cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

int sutra_controller_mv::comp_com() {

  long dims_data2[2];
  dims_data2[0] = 1;
  dims_data2[1] = nslope;
  carma_obj<float> d_tmp(current_context, dims_data2);
  carma_obj<float> d_tmp2(current_context, dims_data2);
  carma_obj<float> d_olmeas(current_context, dims_data2);
  dims_data2[1] = nactu;
  carma_obj<float> d_tmp3(current_context, dims_data2);
  carma_obj<float> d_tmp4(current_context, dims_data2);
  this->d_com2->copy(this->d_com1, 1, 1);
  this->d_com1->copy(this->d_com, 1, 1);
  // POLC equations
  carma_geam<float>(cublas_handle, 'n', 'n', nactu, 1, (float) delay, *d_com2,
      nactu, 1.0f - delay, *d_com1, nactu, d_tmp, nactu);
  //carma_gemv<float>(cublas_handle,'n',nslope,nactu,-1.0f,*d_imat,nslope,d_tmp,1,1.0f,*d_centroids,1);
  carma_gemv<float>(cublas_handle, 'n', nslope, nactu, 1.0f, *d_imat, nslope,
      d_tmp, 1, 0.0f, d_tmp2, 1);
  carma_geam<float>(cublas_handle, 'n', 'n', nslope, 1, 1.0f, *d_centroids,
      nslope, -1.0f, d_tmp2, nslope, d_olmeas, nslope);

  if (this->nstreams > 1) {
    int nstreams = this->nstreams;
    float alpha = -1.0f;
    float beta = 0.0f;
    //cout << this->d_cmat->getDims(1) << " x " << this->d_cmat->getDims(2) << endl;
    //cout << this->d_centroids->getNbElem()<< endl;
    //cout << this->d_err->getNbElem()<< endl;
    for (int i = 0; i < nstreams; i++) {
      int istart1 = i * this->d_cmat->getDims(2) * this->d_cmat->getDims(1)
          / nstreams;
      int istart2 = i * this->d_cmat->getDims(1) / nstreams;

      //cout << istart1 << " " << istart2 << endl;

      cublasSetStream(cublas_handle, this->streams->get_stream(i));

      cublasOperation_t trans = carma_char2cublasOperation('n');

      // carma_checkCublasStatus(
      //     cublasSgemv(cublas_handle, trans, this->d_cmat->getDims(1) / nstreams,
      //         this->d_cmat->getDims(2), &alpha,
      //         &((this->d_cmat->getData())[istart1]),
      //         this->d_cmat->getDims(1) / nstreams, this->d_centroids->getData(),
      //         1, &beta, &((this->d_err->getData())[istart2]), 1));
      carma_checkCublasStatus(
          cublasSgemv(cublas_handle, trans, this->d_cmat->getDims(1) / nstreams,
              this->d_cmat->getDims(2), &alpha,
              &((this->d_cmat->getData())[istart1]),
              this->d_cmat->getDims(1) / nstreams, d_olmeas, 1, &beta,
              &((this->d_err->getData())[istart2]), 1));
    }

    //mult_int(this->d_com->getData(), this->d_err->getData(),
    //   this->d_gain->getData(), this->gain, this->nactu, this->device,
    //   this->streams);

    this->streams->wait_all_streams();
    //this->d_com->copy(this->d_err,1,1);
  } else {
    // compute error
    //  this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->getDims(1),
    //      this->d_centroids, 1, 0.0f, 1);
    this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->getDims(1),
        &d_olmeas, 1, 0.0f, 1);
    // this->d_com->copy(this->d_err,1,1);
    // apply modal gain & loop gain
    //mult_int(this->d_com->getData(), this->d_err->getData(),
    //    this->d_gain->getData(), this->gain, this->nactu, this->device);
  }

  mult_int(d_tmp4, this->d_err->getData(), this->d_gain->getData(), this->gain,
      this->nactu, this->device);
  mult_int(d_tmp3, this->d_com1->getData(), this->d_gain->getData(),
      1 - (this->gain), this->nactu, this->device);
  carma_geam<float>(cublas_handle, 'n', 'n', nactu, 1, 1.0f, d_tmp4, nactu,
      1.0f, d_tmp3, nactu, *d_com, nactu);

  return EXIT_SUCCESS;
}
