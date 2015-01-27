#ifndef _sutra_controller_mv_H_
#define _sutra_controller_mv_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_wfs.h>
#include <sutra_dm.h>
#include <sutra_centroider.h>
#include <sutra_controller.h>
#include <sutra_ao_utils.h>

using namespace std;

class sutra_controller_mv: public sutra_controller {
public:
  int delay;
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
  carma_obj<float> *d_eigenvals;
  carma_host_obj<float> *h_eigenvals;
  carma_obj<float> *d_U;

  // loop components
  carma_obj<float> *d_cenbuff; // centroids circular buffer
  carma_obj<float> *d_com1; // commands k-1 (for POLC)
  carma_obj<float> *d_com2; // commands k-2 (for POLC)
  carma_obj<float> *d_compbuff; // Buffer for computations
  carma_obj<float> *d_compbuff2;
  carma_obj<float> *d_olmeas; // Open-loop measurements for POLC
  carma_obj<float> *d_err; // current error

  cublasHandle_t cublas_handle;

public:
  sutra_controller_mv(carma_context *context, long nvalid, long nactu,
      long delay);
  sutra_controller_mv(const sutra_controller_mv& controller);
  ~sutra_controller_mv();

  string
  get_type();

  int
  svdec_imat();
  int
  build_cmat(const char *dmtype, char *method);
  int
  build_cmat(float *Dm, float *Dtt, float cond);
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
  // Florian features
  int
  load_noisemat(float *noise);
  int
  do_covmat(sutra_dm *ydm,char *method, int *indx_pup, long dim, float *xpos, float *ypos, long Nkl, float norm, float ampli);
  int
  do_geomat(carma_obj<float> *d_geocov, carma_obj<float> *d_IF, long n_pts, float ampli);
  int
  piston_filt(carma_obj<float> *d_statcov);
  int
  piston_filt_cphim(carma_obj<float> *d_cphim);
  int
  invgen(carma_obj<float> *d_mat, float cond, int job);
 // int
 // do_statmat(float *statcov,long dim, float *xpos, float *ypos, float norm, carma_device *device);
  int
  DDiago(carma_obj<float> *d_statcov, carma_obj<float> *d_geocov);
  int
  load_covmat(float *covmat);
  int
  load_klbasis(float *klbasis);
  int
  compute_Cmm(sutra_atmos *atmos, sutra_sensors *sensors, float *L0, float *cn2, float *alphaX, float *alphaY, float diamTel, float cobs);
  int
  compute_Cphim(sutra_atmos *atmos, sutra_sensors *sensors, sutra_dms *dms, float *L0, float *cn2, float *alphaX, float *alphaY, float *X, float *Y, float *xactu, float *yactu, float diamTel, float k2, float *Nact);
};

#endif // _sutra_controller_mv_H_
