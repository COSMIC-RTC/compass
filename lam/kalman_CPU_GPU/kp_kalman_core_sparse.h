//kp_kalman_core_sparse.h
//Kalman controller (core)


#ifndef __LAM__KP_KALMAN_CORE_SPARSE__H__
#define __LAM__KP_KALMAN_CORE_SPARSE__H__

#include "kp_smatrix.h"
#include <string>


class kp_kalman_core_sparse
{
 public:

   //make initialization and calculate gain
   kp_kalman_core_sparse(const kp_smatrix& D_Mo_,
		  const kp_smatrix& N_Act_,
		  const kp_smatrix& PROJ_,
		  bool isZonal_);
   
   virtual void calculate_gain(real bruit,
			       real k_W,
			       const  kp_matrix&  SigmaV,
			       const  kp_vector& atur_,
			       const  kp_vector& btur_) = 0;
  

   virtual void next_step(const kp_vector& Y_k, kp_vector& U_k) = 0;

   virtual ~kp_kalman_core_sparse(){}

 protected:
   //kp_smatrix D_Mo, N_Act, PROJ;
   //kp_vector atur;
   //kp_vector btur;
   int ordreAR;
   int nb_p, nb_act, nb_z, nb_az, nb_n;
   bool gainComputed;
   bool isZonal;
};


#endif
