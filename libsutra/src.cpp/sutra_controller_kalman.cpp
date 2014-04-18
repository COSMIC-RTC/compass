//sutra_controller_kalman.cpp

#include "sutra_controller_kalman.h"

carma_obj<double>* calculate_D_Mo(carma_context* context, int n_slopes) {
  long dims[] = { 1, n_slopes };
  carma_obj<double>* D_Mo = new carma_obj<double>(context, dims);
  return D_Mo;
}
carma_obj<double>* calculate_N_Act(carma_context* context, int n_actus) {
  long dims[] = { 1, n_actus };
  carma_obj<double>* N_Act = new carma_obj<double>(context, dims);
  return N_Act;
}
carma_obj<double>* calculate_PROJ() {
  return NULL;
}
carma_obj<double>* calculate_btur() {
  return NULL;
}
carma_obj<double>* calculate_SigmaV() {
  return NULL;
}
carma_obj<double>* calculate_atur() {
  return NULL;
}

sutra_controller_kalman::sutra_controller_kalman(carma_context* context,
    carma_obj<double>* cD_Mo, carma_obj<double>* cN_Act, carma_obj<double>* cPROJ,
    bool is_zonal) :
    sutra_controller(context, cD_Mo->getDims(1), cN_Act->getDims(1)) {
}

sutra_controller_kalman::~sutra_controller_kalman() {
}

void sutra_controller_kalman::calculate_gain(double bruit, double k_W,
    carma_obj<double>* cSigmaV, carma_obj<double>* catur,
    carma_obj<double>* cbtur) {
}
//                                                                                     
int sutra_controller_kalman::comp_com() {
  return -378;
}
//                                                                                     
