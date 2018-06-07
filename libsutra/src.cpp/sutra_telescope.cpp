#include <sutra_ao_utils.h>
#include <sutra_phase.h>
#include <sutra_telescope.h>

sutra_telescope::sutra_telescope(carma_context *current_context, long n_pup,
                                 long npos, float *pupil, float *phase_ab_M1,
                                 long n_pup_m, float *pupil_m,
                                 float *phase_ab_M1_m) {
  this->current_context = current_context;
  this->device = current_context->get_activeDevice();

  this->pup_size = n_pup;
  this->pup_size_m = n_pup_m;

  this->num_eleme_pup = npos;

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->pup_size;
  dims_data2[2] = this->pup_size;

  this->d_pupil = new carma_obj<float>(this->current_context, dims_data2);
  this->d_pupil->host2device(pupil);

  this->d_phase_ab_M1 = new carma_obj<float>(this->current_context, dims_data2);
  this->d_phase_ab_M1->host2device(phase_ab_M1);

  long *dims_data3 = new long[3];
  dims_data3[0] = 2;
  dims_data3[1] = this->pup_size_m;
  dims_data3[2] = this->pup_size_m;

  this->d_pupil_m = new carma_obj<float>(this->current_context, dims_data3);
  this->d_pupil_m->host2device(pupil_m);

  this->d_phase_ab_M1_m =
      new carma_obj<float>(this->current_context, dims_data3);
  this->d_phase_ab_M1_m->host2device(phase_ab_M1_m);

  delete[] dims_data2;
  delete[] dims_data3;
}

sutra_telescope::~sutra_telescope() {
  // delete this->current_context;
  current_context->set_activeDevice(device, 1);
  delete this->d_pupil;
  delete this->d_phase_ab_M1;
  delete this->d_pupil_m;
  delete this->d_phase_ab_M1_m;
}
