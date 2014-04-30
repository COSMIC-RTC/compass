//sutra_controller_kalman.cpp

#include "sutra_controller_kalman.h"

carma_obj<float>* calculate_D_Mo(carma_context* context, int n_slopes) {
  long dims[] = { 2, n_slopes, 877 }; //877 pour cas particulier en 16m
  ifstream vec_stream("/home/tgautrais/Documents/v2_kalman_kp/data/D_Mo_16P32_AR1Za980.dat");
  istream_iterator<float> start(vec_stream), eos;
  vector<float> vec(start, eos); 
  carma_obj<float>* D_Mo = new carma_obj<float>(context, dims, &vec[0]);
  return D_Mo;
}
carma_obj<float>* calculate_N_Act(carma_context* context, int n_actus) {
  long dims[] = { 2, n_actus, n_actus };

  ifstream vec_stream("/home/tgautrais/Documents/v2_kalman_kp/data/N_Act_16P32_AR1Za980.dat");
  istream_iterator<float> start(vec_stream), eos;
  vector<float> vec(start, eos); 
  carma_obj<float>* N_Act = new carma_obj<float>(context, dims, &vec[0]);  
  return N_Act;
}
carma_obj<float>* calculate_PROJ(carma_context* context, int n_actus) {
  long dims[] = { 2, n_actus, n_actus };
  ifstream vec_stream("/home/tgautrais/Documents/v2_kalman_kp/data/PROJ_16P32_AR1Za980.dat");
  istream_iterator<float> start(vec_stream), eos;
  vector<float> vec(start, eos); 
  carma_obj<float>* PROJ = new carma_obj<float>(context, dims, &vec[0]);
  return PROJ;
}
carma_obj<float>* calculate_btur(carma_context* context, int n_actus) {
  long dims[] = { 1, n_actus };
  float* zeros_tmp = new float[n_actus];
  for (int i=0 ; i<n_actus ; i++) zeros_tmp[i] = 0;
  carma_obj<float>* btur = new carma_obj<float>(context, dims, zeros_tmp);
  delete [] zeros_tmp ; zeros_tmp=NULL;
  return btur;
}
carma_obj<float>* calculate_SigmaV(carma_context* context, int n_actus) {
  long dims[] = { 2, n_actus, n_actus };
  ifstream vec_stream("/home/tgautrais/Documents/v2_kalman_kp/data/SigmaV_16P32_AR1Za980.dat");
  istream_iterator<float> start(vec_stream), eos;
  vector<float> vec(start, eos); 
  carma_obj<float>* SigmaV = new carma_obj<float>(context, dims, &vec[0]);

  return SigmaV;
}
carma_obj<float>* calculate_atur(carma_context* context, int n_actus) {
  long dims[] = { 1, n_actus };
  ifstream vec_stream("/home/tgautrais/Documents/v2_kalman_kp/data/atur_16P32_AR1Za980.dat");
  istream_iterator<float> start(vec_stream), eos;
  vector<float> vec(start, eos); 
  carma_obj<float>* atur = new carma_obj<float>(context, dims, &vec[0]);
  return atur;
}
sutra_controller_kalman::sutra_controller_kalman(carma_context* context,
    carma_obj<float>& cD_Mo, carma_obj<float>& cN_Act, carma_obj<float>& cPROJ,
    bool is_zonal) :
    sutra_controller(context, cD_Mo.getDims(1), cN_Act.getDims(1)) {
}

sutra_controller_kalman::~sutra_controller_kalman() {
}

void sutra_controller_kalman::calculate_gain(double bruit, double k_W,
    carma_obj<float>& cSigmaV, carma_obj<float>& catur,
    carma_obj<float>& cbtur) {
}
//                                                                                     
int sutra_controller_kalman::comp_com() {
  return -378;
}
//                                                                                     
