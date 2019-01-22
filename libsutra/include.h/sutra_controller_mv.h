#ifndef _sutra_controller_mv_H_
#define _sutra_controller_mv_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_centroider.h>
#include <sutra_controller.h>
#include <sutra_dm.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>

template <typename T>
class sutra_controller_mv : public sutra_controller<T> {
 public:
  T gain;

  carma_obj<T> *d_imat;
  carma_obj<T> *d_cmat;
  carma_obj<T> *d_gain;

  // Cphim & Cmm features
  carma_obj<T> *d_covmat;
  carma_obj<T> *d_KLbasis;
  carma_obj<T> *d_noisemat;
  carma_obj<T> *d_Cmm;
  carma_obj<T> *d_Cphim;
  // svd computations
  carma_host_obj<T> *h_Cmmeigenvals;
  carma_host_obj<T> *h_eigenvals;
  // carma_obj<T> *d_U;

  // loop components
  carma_obj<T> *d_cenbuff;   // centroids circular buffer
  carma_obj<T> *d_com1;      // commands k-1 (for POLC)
  carma_obj<T> *d_com2;      // commands k-2 (for POLC)
  carma_obj<T> *d_compbuff;  // Buffer for computations
  carma_obj<T> *d_compbuff2;
  carma_obj<T> *d_olmeas;  // Open-loop measurements for POLC
  carma_obj<T> *d_err;     // current error

  cublasHandle_t cublas_handle;

 public:
  sutra_controller_mv(carma_context *context, long nvalid, long nslope,
                      long nactu, T delay, sutra_dms *dms, int *idx_dms,
                      int ndm);
  sutra_controller_mv(const sutra_controller_mv &controller);
  ~sutra_controller_mv();

  string get_type();

  int svdec_imat();
  int build_cmat(const char *dmtype, char *method);
  int build_cmat(T cond);
  int frame_delay();
  int comp_com();
  int set_gain(T gain);
  int set_mgain(T *mgain);
  int set_delay(T delay);
  int set_cmat(T *cmat);
  int set_imat(T *imat);
  // Florian features
  int load_noisemat(T *noise);
  int do_covmat(sutra_dm *ydm, char *method, int *indx_pup, long dim, T *xpos,
                T *ypos, long Nkl, T norm, T ampli);
  int do_geomat(carma_obj<T> *d_geocov, carma_obj<T> *d_IF, long n_pts,
                T ampli);
  int piston_filt(carma_obj<T> *d_statcov);
  int piston_filt_cphim(carma_obj<T> *d_cphim, T *F);
  int filter_cphim(T *F, T *Nact);
  int filter_cmat(T cond);
  int invgen(carma_obj<T> *d_mat, T cond, int job);
  int invgen(carma_obj<T> *d_mat, carma_host_obj<T> *h_eigen, T cond);
  int invgen_cpu(carma_obj<T> *d_mat, carma_host_obj<T> *h_eigen, T cond);
  // int
  // do_statmat(T *statcov,long dim, T *xpos, T *ypos, T norm,
  // carma_device *device);
  int DDiago(carma_obj<T> *d_statcov, carma_obj<T> *d_geocov);
  int load_covmat(T *covmat);
  int load_klbasis(T *klbasis);
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
