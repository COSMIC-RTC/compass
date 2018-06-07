// sutra_controller_kalman.h

#ifndef __LAM__SUTRA_CONTROLLER_KALMAN_H__
#define __LAM__SUTRA_CONTROLLER_KALMAN_H__

#include "sutra_controller.h"

// include temporaires pour tests avec lesctures des matrices
// a partir des fichiers .dat
#include <fstream>
#include <iterator>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#define __SP setprecision(20) <<

#ifdef COMPILATION_LAM
#include "kp_KFPP.h"
#include "kp_vector.h"
#endif

class kp_kalman_core_sparse;
class kp_kalman_core_full;

class sutra_controller_kalman : public sutra_controller {
 public:
  sutra_controller_kalman(carma_context* context, int nvalid_, int nactu_,
                          sutra_dms* dms, char** type, float* alt, int ndm);

  ~sutra_controller_kalman();

  void init_kalman(carma_host_obj<float>& D_Mo, carma_host_obj<float>& N_Act,
                   carma_host_obj<float>& PROJ, bool is_zonal, bool is_sparse,
                   bool is_GPU);

  double gettime();
  double gettime_op1();
  double gettime_op2();
  double gettime_op3();

  void calculate_gain(float bruit, carma_host_obj<float>& SigmaV,
                      carma_host_obj<float>& atur, carma_host_obj<float>& btur);

  virtual string get_type() {
    if (isInit) {
      if (isGPU)
        return "kalman_GPU";
      else
        return "kalman_CPU";
    } else
      return "kalman_uninitialized";
  };

  virtual int comp_com();

  int set_gain(float k_W);

 private:
  cusparseHandle_t cusparseHandle;
  kp_kalman_core_sparse* core_sparse;
  kp_kalman_core_full* core_full;
  bool isGPU;
  bool isSparse;
  bool isZonal;
  bool isInit;
  bool isGainSet;
  float gain;
};
#endif  //__LAM__SUTRA_CONTROLLER_KALMAN_H__
