//sutra_controller_kalman.cpp

#include "sutra_controller_kalman.h"
//#include "magma.h"

#ifdef COMPILATION_LAM

#include "kp_kalman_core_sparse_GPU.h"
#include "kp_kalman_core_full_GPU.h"
#include "kp_kalman_core_sparse_CPU.h"
#include "kp_kalman_core_full_CPU.h"
#include "kp_smatrix.h"
#include "kp_carma_tools.h"
#include <fstream>

#include "kp_vector.h"

#include "kp_randn.h"
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <cuda.h>

sutra_controller_kalman::sutra_controller_kalman(carma_context* context_, int nvalid_, int nactu_) : sutra_controller(context_, nvalid_ * 2, nactu_) {
   core_sparse = NULL;
   core_full = NULL;
   cusparseHandle = NULL;
   isGPU = true;
   isZonal = true;
   isSparse = true;
   isInit = false;
   isGainSet = false;
}



void sutra_controller_kalman::init_kalman(carma_host_obj<float>& chD_Mo, carma_host_obj<float>& chN_Act, carma_host_obj<float>& chPROJ, bool is_zonal, bool is_sparse, bool is_GPU) {
  current_context->set_activeDevice(device);
   if (isInit)
   {
           cerr << "Error |sutra_controller_kalman::init_kalman | Kalman controller has already been initialized"<<endl;
           exit(EXIT_FAILURE);
   }

   core_sparse = NULL;
   core_full = NULL;
   isGPU = is_GPU;
   isZonal = is_zonal;
   isSparse = is_sparse;

   //convert from carma_host_obj to kp_matrix
   kp_matrix<KFPP> kD_Mo, kN_Act, kPROJ;
   kp_carma_host_obj_to_kp_matrix(chD_Mo,  kD_Mo);
   kp_carma_host_obj_to_kp_matrix(chN_Act, kN_Act);
   kp_carma_host_obj_to_kp_matrix(chPROJ,  kPROJ);


   if (is_sparse)
   {
      //convert from kp_matrix to kp_smatrix (full -> sparse)
      kp_smatrix<KFPP> sD_Mo, sN_Act, sPROJ;
      sD_Mo.init_from_matrix(kD_Mo,pow(10.0,-12.0));sD_Mo.resize2colMajor();
      sN_Act.init_from_matrix(kN_Act,pow(10.0,-12.0));sN_Act.resize2colMajor();
      sPROJ.init_from_matrix(kPROJ,pow(10.0,-12.0));sPROJ.resize2colMajor();
      
      if (is_GPU)
      {
         cusparseStatus_t cusparseStat = cusparseCreate(&cusparseHandle);
         if (cusparseStat != CUSPARSE_STATUS_SUCCESS)
         { 
            cerr<<"Error | sutra_controller_kalman::sutra_controller_kalman  | cusparseCreate failed "<<endl;
            exit(EXIT_FAILURE);
         }

         core_sparse = new kp_kalman_core_sparse_GPU(sD_Mo, sN_Act, sPROJ,
					       is_zonal, 
					       current_context->get_cublasHandle(),
					       cusparseHandle);
      }
      else 
         core_sparse = new kp_kalman_core_sparse_CPU(sD_Mo, sN_Act, sPROJ, is_zonal);

   }  
   else if (is_GPU)
      core_full = new kp_kalman_core_full_GPU(kD_Mo, kN_Act, kPROJ, is_zonal, current_context->get_cublasHandle());

   else
      core_full = new kp_kalman_core_full_CPU(kD_Mo, kN_Act, kPROJ, is_zonal);

   isInit = true;
}

sutra_controller_kalman::~sutra_controller_kalman() {
  current_context->set_activeDevice(device);
   if(core_full)
   {
      delete core_full;
      core_full=NULL;
   }
   if(core_sparse)
   {
      delete core_sparse;
      core_sparse=NULL;
   }

   if (cusparseHandle)
   {
      cusparseDestroy(cusparseHandle);
      cusparseHandle = NULL;
   }
}

void sutra_controller_kalman::calculate_gain(float bruit,
    carma_host_obj<float>& chSigmaV, carma_host_obj<float>& chatur,
    carma_host_obj<float>& chbtur) {
  current_context->set_activeDevice(device);

   if (!isInit)
   {
      cerr << "Error | sutra_controller_kalman::calculate_gain | Kalman controller has not been initialiez"<<endl;
      exit(EXIT_FAILURE);
   }
   //convert carma_obj to kp_matrix
   kp_matrix<double> kSigmaV;
   kp_carma_host_obj_to_kp_matrix(chSigmaV, kSigmaV);
   
   //convert carma_obj to kp_vector
   kp_vector<double> katur, kbtur;
   kp_carma_host_obj_to_kp_vector(chatur, katur);
   kp_carma_host_obj_to_kp_vector(chbtur, kbtur);

   //gain (attribut de la classe) correspond a k_W
   if (!isGainSet)
   {
      cerr << "k_W (gain parameter) has not been set" << endl;
      exit(EXIT_FAILURE);
   }
   if (core_sparse)
   {
           core_sparse->calculate_gain(bruit, gain, kSigmaV, katur, kbtur);
   }
   else if (core_full)
      core_full->calculate_gain(bruit, gain, kSigmaV, katur, kbtur);

}

double sutra_controller_kalman::gettime(){
  current_context->set_activeDevice(device);
  if(core_sparse)
      return core_sparse->temps_boucle.rez();
   else if(core_full)
      return core_full->temps_boucle.rez();
   else return -1;
}
double sutra_controller_kalman::gettime_op1(){
  current_context->set_activeDevice(device);
  if(core_sparse)
      return core_sparse->temps_boucle_op1.rez();
   else if(core_full)
      return core_full->temps_boucle_op1.rez();
   else return -1;
}
double sutra_controller_kalman::gettime_op2(){
  current_context->set_activeDevice(device);
 if(core_sparse)
      return core_sparse->temps_boucle_op2.rez();
   else if(core_full)
      return core_full->temps_boucle_op2.rez();
   else return -1;
}
double sutra_controller_kalman::gettime_op3(){
  current_context->set_activeDevice(device);
   if(core_sparse)
      return core_sparse->temps_boucle_op3.rez();
   else if(core_full)
      return core_full->temps_boucle_op3.rez();
   else return -1;
}

int sutra_controller_kalman::comp_com() {
   
  current_context->set_activeDevice(device);
   kp_vector<KFPP> Y_k,Y_k_tmp, U_k;
   kp_carma_obj_to_kp_vector(*d_centroids, Y_k); 
   if (core_sparse)
      core_sparse->next_step(Y_k, U_k);
   else if (core_full)
      core_full->next_step(Y_k, U_k);
   kp_kp_vector_to_carma_obj(U_k, *d_com);
  
   return -378;
}

#else
sutra_controller_kalman::sutra_controller_kalman(carma_context* context_,
		int nvalid_, int nactu_, sutra_dms *dms, char **type, float *alt,
		int ndm) :
		sutra_controller(context_, nvalid_ * 2, nactu_,0.0f, dms, type, alt, ndm) {
  current_context->set_activeDevice(device,1);
   core_sparse = NULL;
   core_full = NULL;
   cusparseHandle = NULL;
   isGPU = true;
   isZonal = true;
   isSparse = true;
   isInit = false;
}

void sutra_controller_kalman::init_kalman(carma_host_obj<float>& chD_Mo, 
		carma_host_obj<float>& chN_Act, carma_host_obj<float>& chPROJ, 
		bool is_zonal, bool is_sparse, bool is_GPU) {
}
sutra_controller_kalman::~sutra_controller_kalman() {
}

void sutra_controller_kalman::calculate_gain(float bruit,
    carma_host_obj<float>& chSigmaV, carma_host_obj<float>& chatur,
    carma_host_obj<float>& chbtur) {
}
//
int sutra_controller_kalman::comp_com() {
  return -378;
}

//int sutra_controller_kalman::frame_delay() {}
#endif

int sutra_controller_kalman::set_gain(float k_W) {
  this->gain = k_W;
   isGainSet = true;
  return EXIT_SUCCESS;
}


//                                                                                     
