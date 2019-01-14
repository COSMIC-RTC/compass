#ifndef _sutra_controller_generic_H_
#define _sutra_controller_generic_H_

#include <sutra_acquisim.h>
#include <sutra_controller.h>

class sutra_controller_generic : public sutra_controller {
 public:
  float gain;
  carma_obj<float> *d_matE;
  carma_obj<float> *d_cmat;
  carma_obj<float> *d_gain;
  carma_obj<float> *d_decayFactor;
  carma_obj<float> *d_compbuff;

  string command_law;

 public:
  sutra_controller_generic(carma_context *context, long nvalid, long nslope,
                           long nactu, float delay, sutra_dms *dms,
                           int *idx_dms, int ndm);
  sutra_controller_generic(const sutra_controller_generic &controller);
  ~sutra_controller_generic();

  string get_type();
  string get_commandlaw();
  int set_decayFactor(float *decayFactor);
  int set_mgain(float *gain);
  int set_gain(float gain);
  int set_cmat(float *cmat);
  int set_matE(float *matE);
  int set_commandlaw(string law);
  int comp_com();
};

#endif  // _sutra_controller_generic_H_
