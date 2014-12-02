//kp_kalman_core_sparse.h
//Kalman controller (core)


#ifndef __LAM__KP_KALMAN_CORE_SPARSE__H__
#define __LAM__KP_KALMAN_CORE_SPARSE__H__

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include "kp_KFPP.h"
#include "kp_cu_timer.h"
#include <string>

#include "kp_vector.h"
#include "kp_matrix.h"
#include "kp_smatrix.h"

using namespace std;

class kp_kalman_core_sparse
{
 public:

   //make initialization and calculate gain
   kp_kalman_core_sparse(const kp_smatrix<KFPP>& D_Mo_,
		  const kp_smatrix<KFPP>& N_Act_,
		  const kp_smatrix<KFPP>& PROJ_,
		  bool isZonal_);
   
   virtual void calculate_gain(float bruit,
			       float k_W,
			       const  kp_matrix<double>&  SigmaV,
			       const  kp_vector<double>& atur_,
			       const  kp_vector<double>& btur_) = 0;
  

   virtual void next_step(const kp_vector<KFPP>& Y_k, kp_vector<KFPP>& U_k) = 0;

   virtual ~kp_kalman_core_sparse(){}

 protected:
 //public:
   //kp_smatrix D_Mo, N_Act, PROJ;
   //kp_vector atur;
   //kp_vector btur;
   int ordreAR;
   int nb_p, nb_act, nb_z, nb_az, nb_n;
 protected :
   bool gainComputed;
   bool isZonal;


    public:
  kp_cu_timer temps_boucle;//, temps_op2, temps_op3 ;
  kp_cu_timer temps_boucle_op1;//, temps_op2, temps_op3 ;
  kp_cu_timer temps_boucle_op2;//, temps_op2, temps_op3 ;
  kp_cu_timer temps_boucle_op3;//, temps_op2, temps_op3 ;
};


#endif
