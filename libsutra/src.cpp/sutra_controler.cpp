#include <sutra_controler.h>
#include <string>

sutra_controler::sutra_controler(carma_context *context, long nvalid,
    long nactu, long delay, int device, const char *typec) {
  this->d_imat = 0L;
  this->d_cmat = 0L;
  this->d_eigenvals = 0L;
  this->h_eigenvals = 0L;
  this->d_cenbuff = 0L;
  this->h_centroids = 0L;
  this->h_err = 0L;
  this->h_parcure = 0L;
  this->h_syscure = 0L;

  this->current_context = context;

  //long *dims_data2 = new long[3]; dims_data2[0] = 2; 
  //long *dims_data3 = new long[4]; dims_data3[0] = 3;

  this->typec = string(typec);
  this->nvalid = nvalid;
  this->nslope = 2 * nvalid;
  this->nactu = nactu;
  this->delay = delay;
  this->device = device;
  this->gain = 0.0f;

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;

  this->nstreams = 1; //nvalid/10;
  while (nactu % this->nstreams != 0)
    nstreams--;

  cerr << "controler uses " << nstreams << " streams" << endl;
  streams = new carma_streams(nstreams);

  if (this->typec == "ls") {
    dims_data2[1] = nslope;
    dims_data2[2] = nactu;
    this->d_imat = new carma_obj<float>(context, dims_data2);

    dims_data2[1] = nactu;
    dims_data2[2] = nslope;
    this->d_cmat = new carma_obj<float>(context, dims_data2);

    dims_data2[1] = dims_data2[2] = nactu;
    d_U = new carma_obj<float>(current_context, dims_data2);

    this->d_cenbuff = 0L;

    dims_data1[1] = nslope < nactu ? nslope : nactu;
    this->d_eigenvals = new carma_obj<float>(context, dims_data1);
    this->h_eigenvals = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);

  }

  dims_data1[1] = nslope;
  this->d_centroids = new carma_obj<float>(context, dims_data1);

  if (this->typec == "cured") {
    dims_data1[1] = nslope;
    this->h_centroids = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);
    dims_data1[1] = nactu;
    this->h_err = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);
  }

  if (delay > 0) {
    dims_data2[1] = nslope;
    dims_data2[2] = delay + 1;
    this->d_cenbuff = new carma_obj<float>(context, dims_data2);
  }

  dims_data1[1] = nactu;
  this->d_com = new carma_obj<float>(context, dims_data1);
  this->d_err = new carma_obj<float>(context, dims_data1);
  this->d_gain = new carma_obj<float>(context, dims_data1);

  delete[] dims_data1;
  delete[] dims_data2;

  cublas_handle = current_context->get_cublasHandle();
  //carma_checkCublasStatus(cublasCreate(&(this->cublas_handle)));
}

sutra_controler::~sutra_controler() {
  current_context->set_activeDevice(device);
  delete this->d_U;

  delete this->streams;
  if (this->d_imat != 0L)
    delete this->d_imat;
  if (this->d_cmat != 0L)
    delete this->d_cmat;
  if (this->d_eigenvals != 0L)
    delete this->d_eigenvals;
  if (this->h_eigenvals != 0L)
    delete this->h_eigenvals;
  delete this->d_gain;

  if (this->h_centroids != 0L)
    delete this->h_centroids;
  if (this->h_err != 0L)
    delete this->h_err;

  delete this->d_centroids;
  if (this->delay > 0)
    delete this->d_cenbuff;
  delete this->d_com;
  delete this->d_err;

  //carma_checkCublasStatus(cublasDestroy(this->cublas_handle));

  //delete this->current_context;
}

int sutra_controler::svdec_imat() {
  // doing U = Dt.D where D is i_mat
  float one = 1., zero = 0.;

  if (carma_syrk(cublas_handle, CUBLAS_FILL_MODE_LOWER, 't', nactu, nslope, one,
      d_imat->getData(), nslope, zero, d_U->getData(), nactu)) {
    return EXIT_FAILURE;
  }

  // we can skip this step syevd use only the lower part
  //fill_sym_matrix('U', d_U->getData(), nactu, nactu * nactu);

  // doing evd of U inplace
  if (carma_syevd<float>('V', d_U, *h_eigenvals) == EXIT_FAILURE) {
    //Case where MAGMA is not compiled

    //We fill the upper matrix part of the matrix
    fill_sym_matrix<float>('U', *d_U, nactu, nactu * nactu);

    carma_obj<float> d_tmp(d_U);
    carma_obj<float> d_tmp2(d_U);

    if (carma_cula_svd<float>(&d_tmp, d_eigenvals, d_U, &d_tmp2) == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
    d_eigenvals->device2host(*h_eigenvals);
    return EXIT_SUCCESS;
  };
  d_eigenvals->host2device(*h_eigenvals);

  return EXIT_SUCCESS;
}

int sutra_controler::set_gain(float gain) {
  this->gain = gain;
  return EXIT_SUCCESS;
}

int sutra_controler::load_mgain(float *mgain) {
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

int sutra_controler::set_delay(int delay) {
  this->delay = delay;
  return EXIT_SUCCESS;
}

int sutra_controler::build_cmat(int nfilt, bool filt_tt) {
  carma_obj<float> *d_eigenvals_inv;
  carma_host_obj<float> *h_eigenvals_inv;
  carma_obj<float> *d_tmp, *d_tmp2;

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;

  dims_data2[1] = dims_data2[2] = nactu;
  d_tmp = new carma_obj<float>(current_context, dims_data2);
  d_tmp2 = new carma_obj<float>(current_context, dims_data2);

  dims_data1[1] = nactu;
  d_eigenvals_inv = new carma_obj<float>(current_context, dims_data1);
  h_eigenvals_inv = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);

  float one = 1., zero = 0.;

  int nb_elem = this->h_eigenvals->getNbElem();
  memset(h_eigenvals_inv->getData(), 0, sizeof(float) * nb_elem);

  // filtering modes
  if (filt_tt)
    nb_elem -= 2;
  for (int cc = nfilt; cc < nb_elem; cc++) {
    float eigenval = h_eigenvals->getData()[cc];
    h_eigenvals_inv->getData()[cc] =
        (fabs(eigenval) > 1.e-6) ? 1.0f / eigenval : 0.f;
  }
  d_eigenvals_inv->host2device(h_eigenvals_inv->getData());

  carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, nactu, nactu, d_U->getData(),
      nactu, d_eigenvals_inv->getData(), one, d_tmp->getData(), nactu);
  carma_gemm(cublas_handle, 'n', 't', nactu, nactu, nactu, one,
      d_tmp->getData(), nactu, d_U->getData(), nactu, zero, d_tmp2->getData(),
      nactu);
  carma_gemm(cublas_handle, 'n', 't', nactu, nslope, nactu, one,
      d_tmp2->getData(), nactu, d_imat->getData(), nslope, zero,
      d_cmat->getData(), nactu);

  delete d_tmp;
  delete d_tmp2;
  delete d_eigenvals_inv;
  delete h_eigenvals_inv;

  return EXIT_SUCCESS;
}

int sutra_controler::build_cmat(int nfilt) {
  return this->build_cmat(nfilt, false);
}

int sutra_controler::frame_delay() {
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

int sutra_controler::comp_com() {
  if (typec == "ls") {
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

        carma_checkCublasStatus(
            cublasSgemv(cublas_handle, trans,
                this->d_cmat->getDims(1) / nstreams, this->d_cmat->getDims(2),
                &alpha, &((this->d_cmat->getData())[istart1]),
                this->d_cmat->getDims(1) / nstreams,
                this->d_centroids->getData(), 1, &beta,
                &((this->d_err->getData())[istart2]), 1));
      }

      mult_int(this->d_com->getData(), this->d_err->getData(),
          this->d_gain->getData(), this->gain, this->nactu, this->device,
          this->streams);

      this->streams->wait_all_streams();

    } else {
      // compute error
      this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->getDims(1),
          this->d_centroids, 1, 0.0f, 1);

      // apply modal gain & loop gain
      mult_int(this->d_com->getData(), this->d_err->getData(),
          this->d_gain->getData(), this->gain, this->nactu, this->device);
    }
  }

  if (this->typec == "cured") {
    h_centroids->cpy_obj(this->d_centroids, cudaMemcpyDeviceToHost);

    cured(this->h_syscure, this->h_parcure, this->h_centroids->getData(),
        this->h_err->getData(), this->gain);

    h_err->cpy_obj(this->d_com, cudaMemcpyHostToDevice);

    //mult_int(this->d_com->getData(),this->d_err->getData(),this->d_gain->getData(),this->gain,this->nactu,this->device);
  }

  return EXIT_SUCCESS;
}

int sutra_controler::init_cured(int nxsubs, int *isvalid) {
  this->h_syscure = cureSystem(nxsubs, this->nvalid, this->nactu, isvalid, 1);
  this->h_parcure = cureInit(this->h_syscure);
  return EXIT_SUCCESS;
}

