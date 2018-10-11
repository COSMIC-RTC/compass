#include <sutra_controller_generic.h>

sutra_controller_generic::sutra_controller_generic(carma_context *context,
                                                   long nvalid, long nslope,
                                                   long nactu, float delay,
                                                   sutra_dms *dms, int *idx_dms,
                                                   int ndm)
    : sutra_controller(context, nvalid, nslope, nactu, delay, dms, idx_dms,
                       ndm) {
  this->command_law = "integrator";
  long dims_data1[2] = {1, nactu};
  this->d_gain = new carma_obj<float>(current_context, dims_data1);
  this->d_decayFactor = new carma_obj<float>(current_context, dims_data1);
  this->d_compbuff = new carma_obj<float>(current_context, dims_data1);
  long dims_data2[3] = {2, nactu, nslope};
  this->d_cmat = new carma_obj<float>(current_context, dims_data2);

  dims_data2[2] = nactu;
  this->d_matE = new carma_obj<float>(current_context, dims_data2);
}

sutra_controller_generic::~sutra_controller_generic() {
  current_context->set_activeDevice(device, 1);
  delete this->d_gain;
  delete this->d_decayFactor;
  delete this->d_cmat;
  delete this->d_matE;
  delete this->d_compbuff;
}

string sutra_controller_generic::get_type() {
  current_context->set_activeDevice(device, 1);
  return "generic";
}

string sutra_controller_generic::get_commandlaw() {
  current_context->set_activeDevice(device, 1);
  return this->command_law;
}

int sutra_controller_generic::set_mgain(float *gain) {
  current_context->set_activeDevice(device, 1);
  this->d_gain->host2device(gain);
  return EXIT_SUCCESS;
}

int sutra_controller_generic::set_decayFactor(float *decayFactor) {
  current_context->set_activeDevice(device, 1);
  this->d_decayFactor->host2device(decayFactor);
  return EXIT_SUCCESS;
}

int sutra_controller_generic::set_matE(float *matE) {
  current_context->set_activeDevice(device, 1);
  this->d_matE->host2device(matE);
  return EXIT_SUCCESS;
}

int sutra_controller_generic::set_cmat(float *cmat) {
  current_context->set_activeDevice(device, 1);
  this->d_cmat->host2device(cmat);
  return EXIT_SUCCESS;
}

int sutra_controller_generic::set_commandlaw(string law) {
  current_context->set_activeDevice(device, 1);
  this->command_law = law;
  return EXIT_SUCCESS;
}

int sutra_controller_generic::comp_com() {
  current_context->set_activeDevice(device, 1);

  // CMAT*s(k)
  carma_gemv(cublas_handle(), 'n', nactu(), nslope(), 1.0f,
             this->d_cmat->getData(), nactu(), this->d_centroids->getData(), 1,
             0.0f, this->d_compbuff->getData(), 1);
  // g*CMAT*s(k)
  mult_vect(this->d_compbuff->getData(), this->d_gain->getData(), -1.0f,
            this->nactu(), this->current_context->get_device(device));

  if (this->command_law == "2matrices")  // v(k) = E.v(k-1)
    carma_gemv(cublas_handle(), 'n', nactu(), nactu(), 1.0f,
               this->d_matE->getData(), nactu(), this->d_com1->getData(), 1,
               0.0f, this->d_com->getData(), 1);
  else  // v(k) = v(k-1)
    this->d_com->copyFrom(this->d_com1->getData(), 1);
  // v(k) = alpha*E*v(k-1)
  mult_vect(this->d_com->getData(), this->d_decayFactor->getData(), 1.0f,
            this->nactu(), this->current_context->get_device(device));
  // v(k) = alpha*E*v(k-1) + g*CMAT*s(k)
  carma_axpy(cublas_handle(), nactu(), 1.0f, this->d_compbuff->getData(), 1,
             this->d_com->getData(), 1);

  return EXIT_SUCCESS;
}
