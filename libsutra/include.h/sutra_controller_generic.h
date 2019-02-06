#ifndef _sutra_controller_generic_H_
#define _sutra_controller_generic_H_

#include <sutra_acquisim.h>
#include <sutra_controller.h>

template <typename T>
class sutra_controller_generic : public sutra_controller<T> {
 public:
  T gain;
  carma_obj<T> *d_matE;
  carma_obj<T> *d_cmat;
  carma_obj<T> *d_gain;
  carma_obj<T> *d_decayFactor;
  carma_obj<T> *d_com1;      // commands k-1 (for POLC)
  carma_obj<T> *d_com2;      // commands k-2 (for POLC)
  carma_obj<T> *d_compbuff;  // Buffer for computations
  carma_obj<T> *d_compbuff2;
  carma_obj<T> *d_olmeas;  // Open-loop measurements for POLC
  carma_obj<T> *d_imat;

  bool polc;

  string command_law;

 public:
  sutra_controller_generic(carma_context *context, long nvalid, long nslope,
                           long nactu, T delay, sutra_dms *dms, int *idx_dms,
                           int ndm);
  sutra_controller_generic(const sutra_controller_generic &controller);
  ~sutra_controller_generic();

  string get_type();
  string get_commandlaw();
  int set_decayFactor(float *decayFactor);
  int set_mgain(float *gain);
  int set_gain(T gain);
  int set_cmat(float *cmat);
  int set_matE(float *matE);
  int set_commandlaw(string law);
  int set_polc(bool p);
  int set_imat(float *imat);
  int comp_polc();
  int comp_com();
};

#endif  // _sutra_controller_generic_H_
