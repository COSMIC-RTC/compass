//sutra_controller_kalman.cpp

#include "sutra_controller_kalman.h"

carma_obj<float>* calculate_D_Mo(carma_context* context, int n_slopes) {
  long dims[] = { 1, n_slopes };
  carma_obj<float>* D_Mo = new carma_obj<float>(context, dims);
  return D_Mo;
}
carma_obj<float>* calculate_N_Act(carma_context* context, int n_actus) {
  long dims[] = { 1, n_actus };
  carma_obj<float>* N_Act = new carma_obj<float>(context, dims);
  return N_Act;
}
carma_obj<float>* calculate_btur() {
  return NULL;
}

sutra_controller_kalman::sutra_controller_kalman(carma_context* context,
    carma_obj<float>* cD_Mo, carma_obj<float>* cN_Act, carma_obj<float>* cPROJ,
    bool is_zonal) :
    sutra_controller(context, cD_Mo->getDims(1), cN_Act->getDims(1)) {
}

sutra_controller_kalman::~sutra_controller_kalman() {
}

void sutra_controller_kalman::calculate_gain(double bruit, double k_W,
    carma_obj<float>* cSigmaV, carma_obj<float>* catur,
    carma_obj<float>* cbtur) {
}
//                                                                                     
int sutra_controller_kalman::comp_com() {
  return -378;
}
//                                                                                     
