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

class sutra_controller_mv : public sutra_controller {
public:
  int nvalid;
  int nslope;
  int nactu;
  int delay;
  float gain;

  carma_obj<float> *d_imat;
  carma_obj<float> *d_cmat;
  carma_obj<float> *d_gain;

  // Florian features
  carma_obj<float> *d_covmat;
  carma_obj<float> *d_KLbasis;

  // svd computations
  carma_obj<float> *d_eigenvals;
  carma_host_obj<float> *h_eigenvals;
  carma_obj<float> *d_U;

  // loop components
  carma_obj<float> *d_cenbuff; // centroids circular buffer
  carma_obj<float> *d_com1; // commands k-1 (for POLC)
  carma_obj<float> *d_com2; // commands k-2 (for POLC)
  //carma_obj<float>          *d_combuff;   // command circular buffer 
  carma_obj<float> *d_err; // current error
  //carma_obj<float>          *d_err;       // error circular buffer  

  carma_streams *streams;
  int nstreams;
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
  build_cmat(const char *dmtype);
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
  load_covmat(float *covmat);
  int
  load_klbasis(float *klbasis);
  int
  comp_mvcom();

};

#endif // _sutra_controller_mv_H_
