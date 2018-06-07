#ifndef _sutra_controller_ls_H_
#define _sutra_controller_ls_H_

#include <sutra_controller.h>

class sutra_controller_ls : public sutra_controller {
 public:
  float gain;

  carma_obj<float> *d_imat;
  carma_obj<float> *d_cmat;
  carma_obj<float> *d_gain;

  // svd computations
  carma_obj<float> *d_eigenvals;
  carma_host_obj<float> *h_eigenvals;
  carma_obj<float> *d_U;

  // loop components
  carma_obj<float> *d_cenbuff;  // centroids circular buffer
  carma_obj<float> *d_err;      // current error

  // Modal optimization components
  int is_modopti;                // Flag for using modal optimization
  int nrec;                      // Number of recorded open slopes measurements
  int nmodes;                    // Number of modes
  float gmin;                    // Gain min
  float gmax;                    // Gain max
  int ngain;                     // Number of tested gains between gmin and gmax
  float Fs;                      // Sampling frequency
  int cpt_rec;                   // Counter for modal gains refresh
  carma_obj<float> *d_M2V;       // Modes to Volts matrix
  carma_obj<float> *d_S2M;       // Slopes to Modes matrix
  carma_obj<float> *d_slpol;     // Open-loop measurements buffer, recorded and
                                 // loaded from Yorick
  carma_obj<float> *d_Hcor;      // Transfer function
  carma_obj<float> *d_com1;      // Command k-1 for POLC
  carma_obj<float> *d_com2;      // Command k-2 for POLC
  carma_obj<float> *d_compbuff;  // Buffer for POLC computation
  carma_obj<float> *d_compbuff2;  // Buffer for POLC computation

 public:
  sutra_controller_ls(carma_context *context, long nvalid, long nactu,
                      float delay, sutra_dms *dms, char **type, float *alt,
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
  int set_gain(float gain);
  int set_mgain(float *mgain);
  int set_cmat(float *cmat);
  int set_delay(float delay);
  int init_modalOpti(int nmodes, int nrec, float *M2V, float gmin, float gmax,
                     int ngain, float Fs);
  int loadOpenLoopSlp(float *ol_slopes);
  int modalControlOptimization();
  int compute_Hcor();
};

#endif  // _sutra_controller_ls_H_
