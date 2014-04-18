//sutra_controller_kalman.h

#ifndef __LAM__SUTRA_CONTROLLER_KALMAN_H__
#define __LAM__SUTRA_CONTROLLER_KALMAN_H__

#include "sutra_controller.h"

class sutra_controller_kalman: public sutra_controller {
public:
  sutra_controller_kalman(carma_context* context, carma_obj<double>* D_Mo,
      carma_obj<double>* N_Act, carma_obj<double>* PROJ, bool is_zonal);

  ~sutra_controller_kalman();

  void calculate_gain(double bruit, double k_W, carma_obj<double>* SigmaV,
      carma_obj<double>* atur, carma_obj<double>* btur);

  virtual string get_type() {
    return "kalman";
  }
  ;

  virtual int comp_com();
};

carma_obj<double>* calculate_D_Mo(carma_context* context, int n_slopes);
carma_obj<double>* calculate_N_Act(carma_context* context, int n_actus);
carma_obj<double>* calculate_PROJ();
carma_obj<double>* calculate_btur();
carma_obj<double>* calculate_SigmaV();
carma_obj<double>* calculate_atur();

#endif //__LAM__SUTRA_CONTROLLER_KALMAN_H__
