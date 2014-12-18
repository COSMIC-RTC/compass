#ifndef _sutra_controller_ls_H_
#define _sutra_controller_ls_H_

#include <sutra_controller.h>

class sutra_controller_ls: public sutra_controller {
public:
  int delay;
  float gain;

  carma_obj<float> *d_imat;
  carma_obj<float> *d_cmat;
  carma_obj<float> *d_gain;

  // svd computations
  carma_obj<float> *d_eigenvals;
  carma_host_obj<float> *h_eigenvals;
  carma_obj<float> *d_U;

  // loop components
  carma_obj<float> *d_cenbuff; // centroids circular buffer
  carma_obj<float> *d_err; // current error

  // Modal optimization components
  int is_modopti; // Flag for using modal optimization
  int nrec;
  int nmodes;
  carma_obj<float> *d_M2V; // Modes to Volts matrix
  carma_obj<float> *d_S2M; // Slopes to Modes matrix
  carma_obj<float> *d_slpol; // Open-loop measurements buffer, recorded and loaded from Yorick

public:
  sutra_controller_ls(carma_context *context, long nvalid, long nactu,
      long delay);
  sutra_controller_ls(const sutra_controller_ls& controller);
  ~sutra_controller_ls();

  string
  get_type();

  int
  svdec_imat();
  int
  build_cmat(int nfilt, bool filt_tt);
  int
  build_cmat(int nfilt);
  int
  frame_delay();
  int
  comp_com();
  int
  set_gain(float gain);
  int
  load_mgain(float *mgain);
  int
  set_delay(int delay);
  int
  init_modalOpti(int nmodes, int nrec);
  int
  loadOpenLoopSlp(float *ol_slopes);

};

#endif // _sutra_controller_ls_H_
