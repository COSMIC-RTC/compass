#ifndef _sutra_controller_ls_H_
#define _sutra_controller_ls_H_

#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class sutra_controller_ls : public sutra_controller<Tcomp, Tout> {
 public:
  carma_obj<Tcomp> *d_imat;
  carma_obj<Tcomp> *d_cmat;
  carma_obj<Tcomp> *d_gain;

  // svd computations
  carma_obj<Tcomp> *d_eigenvals;
  carma_host_obj<Tcomp> *h_eigenvals;
  carma_obj<Tcomp> *d_U;

  // loop components
  carma_obj<Tcomp> *d_cenbuff;  // centroids circular buffer
  carma_obj<Tcomp> *d_err;      // current error

  // Modal optimization components
  int is_modopti;                // Flag for using modal optimization
  int nrec;                      // Number of recorded open slopes measurements
  int nmodes;                    // Number of modes
  Tcomp gmin;                    // Gain min
  Tcomp gmax;                    // Gain max
  int ngain;                     // Number of tested gains between gmin and gmax
  Tcomp Fs;                      // Sampling frequency
  int cpt_rec;                   // Counter for modal gains refresh
  carma_obj<Tcomp> *d_M2V;       // Modes to Volts matrix
  carma_obj<Tcomp> *d_S2M;       // Slopes to Modes matrix
  carma_obj<Tcomp> *d_slpol;     // Open-loop measurements buffer, recorded and
                                 // loaded from Yorick
  carma_obj<Tcomp> *d_Hcor;      // Transfer function
  carma_obj<Tcomp> *d_compbuff;  // Buffer for POLC computation
  carma_obj<Tcomp> *d_compbuff2;  // Buffer for POLC computation

 public:
  sutra_controller_ls(carma_context *context, long nvalid, long nslope,
                      long nactu, float delay, sutra_dms *dms, int *idx_dms,
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
  int set_mgain(Tcomp *mgain);
  int set_cmat(Tcomp *cmat);
  int set_imat(Tcomp *imat);
  int init_modalOpti(int nmodes, int nrec, Tcomp *M2V, Tcomp gmin, Tcomp gmax,
                     int ngain, Tcomp Fs);
  int loadOpenLoopSlp(Tcomp *ol_slopes);
  int modalControlOptimization();
  int compute_Hcor();
};

#endif  // _sutra_controller_ls_H_
