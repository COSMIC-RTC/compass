#include <sutra_controller_generic.h>

template <typename T, typename Tout>
sutra_controller_generic<T, Tout>::sutra_controller_generic(
    carma_context *context, long nvalid, long nslope, long nactu, float delay,
    sutra_dms *dms, int *idx_dms, int ndm, int nstates)
    : sutra_controller<T, Tout>(context, nvalid, nslope, nactu, delay, dms,
                                idx_dms, ndm) {
  this->command_law = "integrator";
  this->nstates = nstates;
  long dims_data1[2] = {1, nactu + nstates};
  this->d_gain = new carma_obj<T>(this->current_context, dims_data1);
  this->d_decayFactor = new carma_obj<T>(this->current_context, dims_data1);
  this->d_compbuff = new carma_obj<T>(this->current_context, dims_data1);
  if (this->d_com != nullptr) delete this->d_com;
  if (this->d_com1 != nullptr) delete this->d_com1;
  if (this->d_com2 != nullptr) delete this->d_com2;
  this->d_com = new carma_obj<T>(this->current_context, dims_data1);
  this->d_com1 = new carma_obj<T>(this->current_context, dims_data1);
  this->d_com2 = new carma_obj<T>(this->current_context, dims_data1);

  this->d_cmatPadded = nullptr;
  dims_data1[1] = this->nslope();
  this->d_olmeas = new carma_obj<T>(this->current_context, dims_data1);
  this->d_compbuff2 = new carma_obj<T>(this->current_context, dims_data1);

  long dims_data2[3] = {2, nactu + nstates, nslope};
  this->d_cmat = new carma_obj<T>(this->current_context, dims_data2);

  this->gain = 0.f;

  dims_data2[2] = nactu + nstates;
  this->d_matE = new carma_obj<T>(this->current_context, dims_data2);
  dims_data2[1] = this->nslope();
  this->d_imat = new carma_obj<T>(this->current_context, dims_data2);

  this->polc = false;

  if (std::is_same<T, half>::value) {
    int m = nactu + nstates;
    int n = nslope;
    while (m % 8 != 0) m++;
    while (n % 8 != 0) n++;
    dims_data2[1] = m;
    dims_data2[2] = n;
    this->d_cmatPadded = new carma_obj<T>(this->current_context, dims_data2);
    this->d_cmatPadded->reset();
  }
}

template <typename T, typename Tout>
sutra_controller_generic<T, Tout>::~sutra_controller_generic() {
  this->current_context->set_activeDevice(this->device, 1);
  delete this->d_gain;
  delete this->d_decayFactor;
  delete this->d_cmat;
  delete this->d_matE;
  delete this->d_compbuff;
}

template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::set_polc(bool p) {
  this->polc = p;
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
string sutra_controller_generic<T, Tout>::get_type() {
  this->current_context->set_activeDevice(this->device, 1);
  return "generic";
}

template <typename T, typename Tout>
string sutra_controller_generic<T, Tout>::get_commandlaw() {
  this->current_context->set_activeDevice(this->device, 1);
  return this->command_law;
}

template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::set_mgain(float *gain) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_gain->host2device(gain);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::set_decayFactor(float *decayFactor) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_decayFactor->host2device(decayFactor);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::set_matE(float *matE) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_matE->host2device(matE);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::set_imat(float *imat) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_imat->host2device(imat);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::set_cmat(float *cmat) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_cmat->host2device(cmat);
  this->fill_cmatPadded();
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::fill_cmatPadded() {
  return fill_cmatPadded_impl();
}

template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::set_commandlaw(string law) {
  this->current_context->set_activeDevice(this->device, 1);
  this->command_law = law;
  return EXIT_SUCCESS;
}
template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::comp_polc() {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_com2->copy(this->d_com1, 1, 1);
  this->d_com1->copy(this->d_com, 1, 1);
  // POLC equations
  this->d_compbuff->reset();
  this->d_compbuff->axpy(T(this->delay - 1), this->d_com2, 1, 1);
  this->d_compbuff->axpy(T(1 - (this->delay - 1)), this->d_com1, 1, 1);
  this->d_olmeas->copy(this->d_centroids, 1, 1);

  carma_gemv(this->cublas_handle(), 'n', this->nslope(), this->nactu(), T(1.0f),
             this->d_imat->getData(), this->nslope(),
             this->d_compbuff->getData(), 1, T(0.0f),
             this->d_compbuff2->getData(), 1);

  this->d_olmeas->axpy(T(-1.0f), this->d_compbuff2, 1, 1);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic<T, Tout>::comp_com() {
  this->current_context->set_activeDevice(this->device, 1);
  carma_obj<T> *centroids;
  T berta = T(1.0f);
  // cudaEvent_t startEv, stopEv;
  // carmaSafeCall(cudaEventCreate(&startEv));
  // carmaSafeCall(cudaEventCreate(&stopEv));
  // float gpuTime;
  int m, n;

  if (this->polc) {
    this->comp_polc();
    centroids = this->d_olmeas;
    berta = T(1.0f - this->gain);
  } else {
    centroids = this->d_centroids;
  }

  carma_obj<T> *cmat;
  if (std::is_same<T, half>::value) {
    cmat = this->d_cmatPadded;
    m = this->d_cmatPadded->getDims(1);
    n = this->d_cmatPadded->getDims(2);
  } else {
    cmat = this->d_cmat;
    m = this->d_cmat->getDims(1);
    n = this->d_cmat->getDims(2);
  }

  if (this->command_law == "integrator") {
    // cublasSetStream(this->cublas_handle(),
    //                 current_context->get_device(device)->get_stream());
    // carmaSafeCall(cudaEventRecord(startEv));
    carma_gemv(this->cublas_handle(), 'n', m, n, (T)(-1 * this->gain),
               cmat->getData(), m, centroids->getData(), 1, berta,
               this->d_com->getData(), 1);
    // carmaSafeCall(cudaEventRecord(stopEv));
    // carmaSafeCall(cudaEventSynchronize(stopEv));
    // carmaSafeCall(cudaEventElapsedTime(&gpuTime, startEv, stopEv));
    // DEBUG_TRACE("GEMV DONE IN %f ms", gpuTime);

  } else {
    // CMAT*s(k)
    carma_gemv(this->cublas_handle(), 'n', this->nactu(), this->nslope(),
               (T)1.0f, this->d_cmat->getData(), this->nactu() + this->nstates,
               centroids->getData(), 1, (T)0.0f, this->d_compbuff->getData(),
               1);
    // g*CMAT*s(k)
    mult_vect(this->d_compbuff->getData(), this->d_gain->getData(), (T)(-1.0f),
              this->nactu() + this->nstates, this->current_context->get_device(this->device));
    if (this->command_law == "modal_integrator") {
      // M2V * g * CMAT * s(k)
      carma_gemv(this->cublas_handle(), 'n', this->nactu() + this->nstates, this->nactu() + this->nstates,
                 (T)1.0f, this->d_matE->getData(), this->nactu() + this->nstates,
                 this->d_compbuff->getData(), 1, berta, this->d_com->getData(),
                 1);
    } else {  // 2matrices
      carma_gemv(this->cublas_handle(), 'n', this->nactu() + this->nstates, this->nactu() + this->nstates,
                 (T)1.0f, this->d_matE->getData(), this->nactu() + this->nstates,
                 this->d_com1->getData() + this->nstates, 1, (T)0.0f, this->d_com->getData(),
                 1);
      // v(k) = alpha*E*v(k-1)
      mult_vect(this->d_com->getData(), this->d_decayFactor->getData(), (T)1.0f,
                this->nactu() + this->nstates, this->current_context->get_device(this->device));
      // v(k) = alpha*E*v(k-1) + g*CMAT*s(k)
      this->d_com->axpy((T)1.0f, this->d_compbuff, 1, 1);
    }
  }
  //   cudaEventDestroy(startEv);
  // cudaEventDestroy(stopEv);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, half>::value, int>::type
sutra_controller_generic<T, Tout>::fill_cmatPadded_impl() {
  pad_cmat(this->d_cmat->getData(), this->d_cmat->getDims(1),
           this->d_cmat->getDims(2), this->d_cmatPadded->getData(),
           this->d_cmatPadded->getDims(1), this->d_cmatPadded->getDims(2),
           this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template class sutra_controller_generic<float, float>;
template class sutra_controller_generic<float, uint16_t>;
#ifdef CAN_DO_HALF
template class sutra_controller_generic<half, float>;
template class sutra_controller_generic<half, uint16_t>;
#endif
