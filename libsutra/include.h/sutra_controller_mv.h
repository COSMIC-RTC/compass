#ifndef _sutra_controller_mv_H_
#define _sutra_controller_mv_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_centroider.h>
#include <sutra_controller.h>
#include <sutra_dm.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>

class sutra_controller_mv : public sutra_controller {
 public:
  float gain;

  carma_obj<float> *d_imat;
  carma_obj<float> *d_cmat;
  carma_obj<float> *d_gain;

  // Cphim & Cmm features
  carma_obj<float> *d_covmat;
  carma_obj<float> *d_KLbasis;
  carma_obj<float> *d_noisemat;
  carma_obj<float> *d_Cmm;
  carma_obj<float> *d_Cphim;
  // svd computations
  carma_host_obj<float> *h_Cmmeigenvals;
  carma_host_obj<float> *h_eigenvals;
  // carma_obj<float> *d_U;

  // loop components
  carma_obj<float> *d_cenbuff;   // centroids circular buffer
  carma_obj<float> *d_com1;      // commands k-1 (for POLC)
  carma_obj<float> *d_com2;      // commands k-2 (for POLC)
  carma_obj<float> *d_compbuff;  // Buffer for computations
  carma_obj<float> *d_compbuff2;
  carma_obj<float> *d_olmeas;  // Open-loop measurements for POLC
  carma_obj<float> *d_err;     // current error

  cublasHandle_t cublas_handle;

 public:
  sutra_controller_mv(carma_context *context, long nvalid, long nactu,
                      float delay, sutra_dms *dms, int *idx_dms, int ndm);
  sutra_controller_mv(const sutra_controller_mv &controller);
  ~sutra_controller_mv();

  string get_type();

  int svdec_imat();
  int build_cmat(const char *dmtype, char *method);
  int build_cmat(float cond);
  int frame_delay();
  int comp_com();
  int set_gain(float gain);
  int set_mgain(float *mgain);
  int set_delay(float delay);
  int set_cmat(float *cmat);
  // Florian features
  int load_noisemat(float *noise);
  int do_covmat(sutra_dm *ydm, char *method, int *indx_pup, long dim,
                float *xpos, float *ypos, long Nkl, float norm, float ampli);
  int do_geomat(carma_obj<float> *d_geocov, carma_obj<float> *d_IF, long n_pts,
                float ampli);
  int piston_filt(carma_obj<float> *d_statcov);
  int piston_filt_cphim(carma_obj<float> *d_cphim, float *F);
  int filter_cphim(float *F, float *Nact);
  int filter_cmat(float cond);
  int invgen(carma_obj<float> *d_mat, float cond, int job);
  int invgen(carma_obj<float> *d_mat, carma_host_obj<float> *h_eigen,
             float cond);
  int invgen_cpu(carma_obj<float> *d_mat, carma_host_obj<float> *h_eigen,
                 float cond);
  // int
  // do_statmat(float *statcov,long dim, float *xpos, float *ypos, float norm,
  // carma_device *device);
  int DDiago(carma_obj<float> *d_statcov, carma_obj<float> *d_geocov);
  int load_covmat(float *covmat);
  int load_klbasis(float *klbasis);
  int compute_Cmm(sutra_atmos *atmos, sutra_sensors *sensors, double *L0,
                  double *cn2, double *alphaX, double *alphaY, double diamTel,
                  double cobs);
  int compute_Cphim(sutra_atmos *atmos, sutra_sensors *sensors, sutra_dms *dms,
                    double *L0, double *cn2, double *alphaX, double *alphaY,
                    double *X, double *Y, double *xactu, double *yactu,
                    double diamTel, double *k2, long *NlayerDm,
                    long *indLayerDm, double FoV, double *pitch,
                    double *alt_dm);
};

#endif  // _sutra_controller_mv_H_
