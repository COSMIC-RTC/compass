#include <sutra_phase.h>
#include <sutra_telescope.h>
#include <sutra_utils.h>

sutra_telescope::sutra_telescope(carma_context *current_context, long n_pup,
                                 long npos, float *pupil, long n_pup_m,
                                 float *pupil_m) {
  this->current_context = current_context;
  this->device = current_context->get_activeDevice();

  this->pup_size = n_pup;
  this->pup_size_m = n_pup_m;

  this->num_eleme_pup = npos;
  this->d_phase_ab_M1 = nullptr;
  this->d_phase_ab_M1_m = nullptr;

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->pup_size;
  dims_data2[2] = this->pup_size;

  this->d_pupil = new carma_obj<float>(this->current_context, dims_data2);
  this->d_pupil->host2device(pupil);

  long *dims_data3 = new long[3];
  dims_data3[0] = 2;
  dims_data3[1] = this->pup_size_m;
  dims_data3[2] = this->pup_size_m;

  this->d_pupil_m = new carma_obj<float>(this->current_context, dims_data3);
  this->d_pupil_m->host2device(pupil_m);

  delete[] dims_data2;
  delete[] dims_data3;
}

sutra_telescope::~sutra_telescope() {
  // delete this->current_context;
  current_context->set_activeDevice(device, 1);
  delete this->d_pupil;
  delete this->d_pupil_m;
  if (this->d_phase_ab_M1 != nullptr) delete this->d_phase_ab_M1;
  if (this->d_phase_ab_M1_m != nullptr) delete this->d_phase_ab_M1_m;
}

int sutra_telescope::set_phase_ab_M1(float *phase_ab_M1, int size) {
  current_context->set_activeDevice(device, 1);
  if (size == this->pup_size * this->pup_size) {
    long *dims_data2 = new long[3];
    dims_data2[0] = 2;
    dims_data2[1] = this->pup_size;
    dims_data2[2] = this->pup_size;
    if (this->d_phase_ab_M1 == nullptr)
      this->d_phase_ab_M1 =
          new carma_obj<float>(this->current_context, dims_data2, phase_ab_M1);
    else
      this->d_phase_ab_M1->host2device(phase_ab_M1);
  } else
    DEBUG_TRACE("Wrong dimensions");

  return EXIT_SUCCESS;
}

int sutra_telescope::set_phase_ab_M1_m(float *phase_ab_M1_m, int size) {
  current_context->set_activeDevice(device, 1);
  if (size == this->pup_size_m * this->pup_size_m) {
    long *dims_data2 = new long[3];
    dims_data2[0] = 2;
    dims_data2[1] = this->pup_size_m;
    dims_data2[2] = this->pup_size_m;
    if (this->d_phase_ab_M1_m == nullptr)
      this->d_phase_ab_M1_m = new carma_obj<float>(this->current_context,
                                                   dims_data2, phase_ab_M1_m);
    else
      this->d_phase_ab_M1_m->host2device(phase_ab_M1_m);
  } else
    DEBUG_TRACE("Wrong dimensions");

  return EXIT_SUCCESS;
}
