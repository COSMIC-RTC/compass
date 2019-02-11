#ifndef _sutra_controller_generic_H_
#define _sutra_controller_generic_H_

#include <sutra_acquisim.h>
#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class sutra_controller_generic : public sutra_controller<Tcomp, Tout> {
 public:
  Tcomp gain;
  carma_obj<Tcomp> *d_matE;
  carma_obj<Tcomp> *d_cmat;
  carma_obj<Tcomp> *d_gain;
  carma_obj<Tcomp> *d_decayFactor;
  carma_obj<Tcomp> *d_compbuff;  // Buffer for computations
  carma_obj<Tcomp> *d_compbuff2;
  carma_obj<Tcomp> *d_olmeas;  // Open-loop measurements for POLC
  carma_obj<Tcomp> *d_imat;

  bool polc;

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
  int set_polc(bool p);
  int set_imat(float *imat);
  int comp_polc();
  int comp_com();
};

#endif  // _sutra_controller_generic_H_
