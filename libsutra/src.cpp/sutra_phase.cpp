#include <sutra_phase.h>

sutra_phase::sutra_phase(carma_context *current_context, long size) {
  this->current_context = current_context;
  this->screen_size = size;
  this->mat = 0;
  this->zernikes = 0;
  this->zerCoeff = 0;
  this->device = current_context->get_activeDevice();

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->screen_size;
  dims_data2[2] = this->screen_size;

  this->d_screen = new carma_obj<float>(current_context, dims_data2);
  delete[] dims_data2;
}

sutra_phase::~sutra_phase() {
  current_context->set_activeDevice(device, 1);
  delete this->d_screen;
}
