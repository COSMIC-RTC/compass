#ifndef _sutra_controller_ls_H_
#define _sutra_controller_ls_H_

#include <sutra_controller.h>

template <typename T>
class sutra_controller_ls : public sutra_controller<T> {
 public:
  T gain;

  carma_obj<T> *d_imat;
  carma_obj<T> *d_cmat;
  carma_obj<T> *d_gain;

  // svd computations
  carma_obj<T> *d_eigenvals;
  carma_host_obj<T> *h_eigenvals;
  carma_obj<T> *d_U;

  // loop components
  carma_obj<T> *d_cenbuff;  // centroids circular buffer
  carma_obj<T> *d_err;      // current error

  // Modal optimization components
  int is_modopti;             // Flag for using modal optimization
  int nrec;                   // Number of recorded open slopes measurements
  int nmodes;                 // Number of modes
  T gmin;                     // Gain min
  T gmax;                     // Gain max
  int ngain;                  // Number of tested gains between gmin and gmax
  T Fs;                       // Sampling frequency
  int cpt_rec;                // Counter for modal gains refresh
  carma_obj<T> *d_M2V;        // Modes to Volts matrix
  carma_obj<T> *d_S2M;        // Slopes to Modes matrix
  carma_obj<T> *d_slpol;      // Open-loop measurements buffer, recorded and
                              // loaded from Yorick
  carma_obj<T> *d_Hcor;       // Transfer function
  carma_obj<T> *d_com1;       // Command k-1 for POLC
  carma_obj<T> *d_com2;       // Command k-2 for POLC
  carma_obj<T> *d_compbuff;   // Buffer for POLC computation
  carma_obj<T> *d_compbuff2;  // Buffer for POLC computation

 public:
  sutra_controller_ls(carma_context *context, long nvalid, long nslope,
                      long nactu, T delay, sutra_dms *dms, int *idx_dms,
                      int ndm);
  sutra_controller_ls(const sutra_controller_ls &controller);
  ~sutra_controller_ls();

  string get_type();

  int svdec_imat();
  int build_cmat(int nfilt, bool filt_tt);
  int build_cmat(int nfilt);
  int build_cmat_modopti();
  int frame_delay();
  int comp_com();
  int set_gain(T gain);
  int set_mgain(T *mgain);
  int set_cmat(T *cmat);
  int set_imat(T *imat);
  int set_delay(T delay);
  int init_modalOpti(int nmodes, int nrec, T *M2V, T gmin, T gmax, int ngain,
                     T Fs);
  int loadOpenLoopSlp(T *ol_slopes);
  int modalControlOptimization();
  int compute_Hcor();
};

#endif  // _sutra_controller_ls_H_
