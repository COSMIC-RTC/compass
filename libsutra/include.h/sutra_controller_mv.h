#ifndef _sutra_controller_mv_H_
#define _sutra_controller_mv_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_centroider.h>
#include <sutra_controller.h>
#include <sutra_dm.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>

template <typename Tcomp, typename Tout>
class sutra_controller_mv : public sutra_controller<Tcomp, Tout> {
 public:
  carma_obj<Tcomp> *d_imat;
  carma_obj<Tcomp> *d_cmat;
  carma_obj<Tcomp> *d_gain;

  // Cphim & Cmm features
  carma_obj<Tcomp> *d_covmat;
  carma_obj<Tcomp> *d_KLbasis;
  carma_obj<Tcomp> *d_noisemat;
  carma_obj<Tcomp> *d_Cmm;
  carma_obj<Tcomp> *d_Cphim;
  // svd computations
  carma_host_obj<Tcomp> *h_Cmmeigenvals;
  carma_host_obj<Tcomp> *h_eigenvals;
  // carma_obj<Tcomp> *d_U;

  // loop components
  carma_obj<Tcomp> *d_cenbuff;   // centroids circular buffer
  carma_obj<Tcomp> *d_compbuff;  // Buffer for computations
  carma_obj<Tcomp> *d_compbuff2;
  carma_obj<Tcomp> *d_olmeas;  // Open-loop measurements for POLC
  carma_obj<Tcomp> *d_err;     // current error

  cublasHandle_t cublas_handle;

 public:
  sutra_controller_mv(carma_context *context, long nvalid, long nslope,
                      long nactu, float delay, sutra_dms *dms, int *idx_dms,
                      int ndm);
  sutra_controller_mv(const sutra_controller_mv &controller);
  ~sutra_controller_mv();

  string get_type();

  int svdec_imat();
  int build_cmat(const char *dmtype, char *method);
  int build_cmat(Tcomp cond);
  int frame_delay();
  int comp_com();
  int set_mgain(Tcomp *mgain);
  int set_cmat(Tcomp *cmat);
  int set_imat(Tcomp *imat);
  // Florian features
  int load_noisemat(Tcomp *noise);
  int do_covmat(sutra_dm *ydm, char *method, int *indx_pup, long dim,
                Tcomp *xpos, Tcomp *ypos, long Nkl, Tcomp norm, Tcomp ampli);
  int do_geomat(carma_obj<Tcomp> *d_geocov, carma_obj<Tcomp> *d_IF, long n_pts,
                Tcomp ampli);
  int piston_filt(carma_obj<Tcomp> *d_statcov);
  int piston_filt_cphim(carma_obj<Tcomp> *d_cphim, Tcomp *F);
  int filter_cphim(Tcomp *F, Tcomp *Nact);
  int filter_cmat(Tcomp cond);
  int invgen(carma_obj<Tcomp> *d_mat, Tcomp cond, int job);
  int invgen(carma_obj<Tcomp> *d_mat, carma_host_obj<Tcomp> *h_eigen,
             Tcomp cond);
  int invgen_cpu(carma_obj<Tcomp> *d_mat, carma_host_obj<Tcomp> *h_eigen,
                 Tcomp cond);
  // int
  // do_statmat(T *statcov,long dim, T *xpos, T *ypos, T norm,
  // carma_device *device);
  int DDiago(carma_obj<Tcomp> *d_statcov, carma_obj<Tcomp> *d_geocov);
  int load_covmat(Tcomp *covmat);
  int load_klbasis(Tcomp *klbasis);
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
