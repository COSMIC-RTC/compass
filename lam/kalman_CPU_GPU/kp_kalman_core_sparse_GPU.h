//kp_kalman_core_sparse_GPU.h
//Kalman controller (core)


#ifndef __LAM__KP_KALMAN_CORE_SPARSE_GPU__H__
#define __LAM__KP_KALMAN_CORE_SPARSE_GPU__H__

#include "kp_kalman_core_sparse.h"
#include <math.h>
#include "kp_randn.h"
#include "kp_timer.h"
#include "kp_cu_timer.h"
#include "kp_cu_smatrix.h"
#include "kp_cu_matrix.h"
#include "kp_tlib.h"
#include "cublas_v2.h"
#include "cusparse_v2.h"
#include "magma.h"
#include "kp_cu_reduce.h"
#include "kp_cu2kp.h"

class kp_cu_matrix;

class kp_kalman_core_sparse_GPU : public kp_kalman_core_sparse
{
 public:

   
   kp_kalman_core_sparse_GPU(const kp_smatrix& D_Mo_,
			       const kp_smatrix& N_Act_,
			       const kp_smatrix& PROJ_,
			       bool isZonal_,
			       cublasHandle_t cublasHandle_,
			       cusparseHandle_t cusparseHandle_);
   

   virtual void calculate_gain(real bruit_pix,
			       real k_W,
			       const  kp_matrix&  SigmaV,
			       const  kp_vector& atur_,
			       const  kp_vector& btur_);
   
   virtual void next_step(const kp_vector& Y_k, kp_vector& U_k);
 


   virtual ~kp_kalman_core_sparse_GPU();
   //void init_for_next_step();

private :
  kp_cu_vector cu_U_km2 ;
  kp_cu_vector cu_U_km1; 
  kp_cu_vector cu_X_kskm1;

  kp_cu_smatrix cu_D_Mo;
  kp_cu_smatrix cu_N_Act;
  kp_cu_smatrix cu_PROJ;
  
  kp_cu_vector cu_atur;
  kp_cu_vector cu_btur;

  kp_cu_matrix cu_H_inf;

  cublasHandle_t cublasHandle;
  cusparseHandle_t cusparseHandle;



//variables de next_step
	kp_cu_vector cu_Y_k;	
	kp_cu_vector cu_U_k;
	kp_cu_vector cu_Nact_Ukm2;
	kp_cu_vector cu_tmp_vec1; 
	kp_cu_vector cu_innovation;
	kp_cu_vector cu_X_kskm1_tmp;		
	kp_cu_vector cu_X_kp1sk;
	kp_cu_vector cu_X_kp1sk_tmp;
	kp_cu_vector cu_Y_kskm1; 
	kp_cu_vector cu_A1_00_Xkdebut;
	kp_cu_vector cu_A1_01_Xkfin;
	kp_cu_vector cu_X_kp1sk_debut;

    //public:
  //kp_cu_timer temps_op1, temps_op2, temps_op3 ;


};


#endif
