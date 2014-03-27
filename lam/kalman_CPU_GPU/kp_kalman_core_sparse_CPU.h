//kp_kalman_core_sparse_CPU.h
//Kalman controller (core)


#ifndef __LAM__KP_KALMAN_CORE_SPARSE_CPU__H__
#define __LAM__KP_KALMAN_CORE_SPARSE_CPU__H__

#include "kp_kalman_core_sparse.h"
#include <math.h>
#include "kp_randn.h"
#include "kp_tlib.h"
#include "kp_timer.h"


class kp_kalman_core_sparse_CPU : public kp_kalman_core_sparse
{
 public:

   
   kp_kalman_core_sparse_CPU(const kp_smatrix& D_Mo_,
			       const kp_smatrix& N_Act_,
			       const kp_smatrix& PROJ_,
			       bool isZonal_);

   virtual void calculate_gain(real bruit_pix,
			       real k_W,
			       const  kp_matrix&  SigmaV,
			       const  kp_vector& atur_,
			       const  kp_vector& btur_);
  

   virtual void next_step(const kp_vector& Y_k, kp_vector& U_k);

   virtual ~kp_kalman_core_sparse_CPU(){}


private	:
  kp_vector U_km2 ;
  kp_vector U_km1; 
  kp_vector X_kskm1;

  kp_smatrix D_Mo, N_Act, PROJ;

  kp_vector atur;
  kp_vector btur;
  kp_matrix H_inf;

  //variables de next_step
  kp_vector Nact_Ukm2;
  kp_vector tmp_vec1;
  kp_vector Y_kskm1;
  kp_vector innovation;
  kp_vector X_kskm1_tmp; 
  kp_vector A1_00_Xkdebut;
  kp_vector A1_01_Xkfin;
  kp_vector X_kp1sk_debut;
  kp_vector X_kp1sk_tmp;	
  kp_vector X_kp1sk;

//public:
  //kp_timer temps_op1, temps_op2, temps_op3 ;

     
};


#endif
