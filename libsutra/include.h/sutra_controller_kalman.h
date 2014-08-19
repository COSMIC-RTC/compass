//sutra_controller_kalman.h

#ifndef __LAM__SUTRA_CONTROLLER_KALMAN_H__
#define __LAM__SUTRA_CONTROLLER_KALMAN_H__

#include "sutra_controller.h"

//include temporaires pour tests avec lesctures des matrices
//a partir des fichiers .dat
#include <iterator>
#include <fstream>


#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#define __SP setprecision(20)<<

#ifdef COMPILATION_LAM
#include "/home/tgautrais/compass_debug/lam/kalman_CPU_GPU/kp_real.h"
#endif

class kp_kalman_core_sparse;
class kp_kalman_core_full;


class sutra_controller_kalman: public sutra_controller {
public:
   sutra_controller_kalman(carma_context* context, int nslope_, int nactu_);

   ~sutra_controller_kalman();
   
   void init_kalman(carma_host_obj<float>& D_Mo, carma_host_obj<float>& N_Act, carma_host_obj<float>& PROJ, bool is_zonal, bool is_sparse, bool is_GPU);

   void calculate_gain(double bruit, carma_host_obj<float>& SigmaV,
		     carma_host_obj<float>& atur, carma_host_obj<float>& btur);
   
   virtual string get_type() {if(isGPU) return "kalman_GPU";else return "kalman_CPU";};
   
   virtual int comp_com();
   //int frame_delay();

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
   //int delay;
   //carma_obj<float> *d_cenbuff; // centroids circular buffer

 public:
#ifdef COMPILATION_LAM
   real** Yk;
#endif
   int ind_Yk;
   bool matrices_matlab;
   bool pentes_matlab;
   bool sigmaVmatlab;

};

carma_obj<float>* calculate_D_Mo(carma_context* context, int nslopes, int n_actu_zern, bool is_zonal);
carma_obj<float>* calculate_N_Act(carma_context* context, int n_actu_zern, int n_actus, bool is_zonal);
carma_obj<float>* calculate_PROJ(carma_context* context, int n_actus, int n_actu_zern, bool is_zonal);
carma_obj<float>* calculate_btur(carma_context* context, int n_actu_zern, bool is_zonal);
carma_obj<float>* calculate_SigmaV(carma_context* context, int n_actu_zern, bool is_zonal);
carma_obj<float>* calculate_atur(carma_context* context, int n_actu_zern, bool is_zonal);

#endif //__LAM__SUTRA_CONTROLLER_KALMAN_H__
