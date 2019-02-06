#include <sutra_controller_generic.h>

template <typename T>
sutra_controller_generic<T>::sutra_controller_generic(carma_context *context,
                                                      long nvalid, long nslope,
                                                      long nactu, T delay,
                                                      sutra_dms *dms,
                                                      int *idx_dms, int ndm)
    : sutra_controller<T>(context, nvalid, nslope, nactu, delay, dms, idx_dms,
                          ndm) {
  this->command_law = "integrator";
  long dims_data1[2] = {1, nactu};
  this->d_gain = new carma_obj<T>(this->current_context, dims_data1);
  this->d_decayFactor = new carma_obj<T>(this->current_context, dims_data1);
  this->d_compbuff = new carma_obj<T>(this->current_context, dims_data1);
  if(this->d_com1 == nullptr) this->d_com1 = new carma_obj<T>(this->current_context, dims_data1);
  if(this->d_com2 == nullptr) this->d_com2 = new carma_obj<T>(this->current_context, dims_data1);
  
  dims_data1[1] = this->nslope();
  this->d_olmeas = new carma_obj<T>(this->current_context, dims_data1);
  this->d_compbuff2 = new carma_obj<T>(this->current_context, dims_data1);

  long dims_data2[3] = {2, nactu, nslope};
  this->d_cmat = new carma_obj<T>(this->current_context, dims_data2);
  this->gain = 0.3f;

  dims_data2[2] = nactu;
  this->d_matE = new carma_obj<T>(this->current_context, dims_data2);
  dims_data2[1] = this->nslope();
  this->d_imat = new carma_obj<T>(this->current_context, dims_data2);

  this->polc = false;

}

template <typename T>
sutra_controller_generic<T>::~sutra_controller_generic() {
  this->current_context->set_activeDevice(this->device, 1);
  delete this->d_gain;
  delete this->d_decayFactor;
  delete this->d_cmat;
  delete this->d_matE;
  delete this->d_compbuff;
}

template <typename T>
int sutra_controller_generic<T>::set_polc(bool p){
    this->polc=p;
  return EXIT_SUCCESS;
}

template <typename T>
string sutra_controller_generic<T>::get_type() {
  this->current_context->set_activeDevice(this->device, 1);
  return "generic";
}

template <typename T>
string sutra_controller_generic<T>::get_commandlaw() {
  this->current_context->set_activeDevice(this->device, 1);
  return this->command_law;
}

template <typename T>
int sutra_controller_generic<T>::set_mgain(float *gain) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_gain->host2device(gain);
  return EXIT_SUCCESS;
}

template <typename T>
int sutra_controller_generic<T>::set_gain(T gain) {
  this->current_context->set_activeDevice(this->device, 1);
  this->gain = gain;
  return EXIT_SUCCESS;
}

template <typename T>
int sutra_controller_generic<T>::set_decayFactor(float *decayFactor) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_decayFactor->host2device(decayFactor);
  return EXIT_SUCCESS;
}

template <typename T>
int sutra_controller_generic<T>::set_matE(float *matE) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_matE->host2device(matE);
  return EXIT_SUCCESS;
}

template <typename T>
int sutra_controller_generic<T>::set_imat(float *imat) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_imat->host2device(imat);
  return EXIT_SUCCESS;
}

template <typename T>
int sutra_controller_generic<T>::set_cmat(float *cmat) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_cmat->host2device(cmat);
  return EXIT_SUCCESS;
}

template <typename T>
int sutra_controller_generic<T>::set_commandlaw(string law) {
  this->current_context->set_activeDevice(this->device, 1);
  this->command_law = law;
  return EXIT_SUCCESS;
}
template <typename T>
int sutra_controller_generic<T>::comp_polc() {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_com2->copy(this->d_com1, 1, 1);
  this->d_com1->copy(this->d_com, 1, 1);
  // POLC equations
  this->d_compbuff->reset();
  this->d_compbuff->axpy(T(this->delay-1), this->d_com2, 1, 1);
  this->d_compbuff->axpy(T(1-(this->delay-1)), this->d_com1, 1, 1);
  this->d_olmeas->copy(this->d_centroids, 1, 1);

  carma_gemv(this->cublas_handle(), 'n', this->nslope(), this->nactu(), T(1.0f), this->d_imat->getData(),
                    this->nslope(), this->d_compbuff->getData(), 1, T(0.0f), this->d_compbuff2->getData(), 1);

  this->d_olmeas->axpy(T(-1.0f), this->d_compbuff2, 1, 1);

}

template <typename T>
int sutra_controller_generic<T>::comp_com() {
  this->current_context->set_activeDevice(this->device, 1);
  carma_obj<T> *centroids;
  T berta = T(1.0f);

  if(this->polc){
    this->comp_polc();
    centroids = this->d_olmeas;
    berta = T(1.0f - this->gain);
  }
  else{
    centroids = this->d_centroids;
  }

  if (this->command_law == "integrator") {
    // cublasSetStream(this->cublas_handle(),
    //                 current_context->get_device(device)->get_stream());
    carma_gemv(this->cublas_handle(), 'n', this->nactu(), this->nslope(),
               (T)(-1 * this->gain), this->d_cmat->getData(), this->nactu(),
               centroids->getData(), 1, berta, this->d_com->getData(),
               1);

  } else {
    // CMAT*s(k)
    carma_gemv(this->cublas_handle(), 'n', this->nactu(), this->nslope(),
               (T)1.0f, this->d_cmat->getData(), this->nactu(),
               centroids->getData(), 1, (T)0.0f,
               this->d_compbuff->getData(), 1);
    // g*CMAT*s(k)
    mult_vect(this->d_compbuff->getData(), this->d_gain->getData(), (T)(-1.0f),
              this->nactu(), this->current_context->get_device(this->device));

    carma_gemv(this->cublas_handle(), 'n', this->nactu(), this->nactu(),
               (T)1.0f, this->d_matE->getData(), this->nactu(),
               this->d_com1->getData(), 1, (T)0.0f, this->d_com->getData(), 1);
    // v(k) = alpha*E*v(k-1)
    mult_vect(this->d_com->getData(), this->d_decayFactor->getData(), (T)1.0f,
              this->nactu(), this->current_context->get_device(this->device));
    // v(k) = alpha*E*v(k-1) + g*CMAT*s(k)
    this->d_com->axpy((T)1.0f, this->d_compbuff, 1, 1);
  }
  return EXIT_SUCCESS;
}

template class sutra_controller_generic<float>;
#ifdef CAN_DO_HALF
template class sutra_controller_generic<half>;
#endif
