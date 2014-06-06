//sutra_controller_kalman.cpp

#include "sutra_controller_kalman.h"

#ifdef COMPILATION_LAM
#include "kp_kalman_core_sparse_GPU.h"
#include "kp_kalman_core_full_GPU.h"
#include "kp_kalman_core_sparse_CPU.h"
#include "kp_kalman_core_full_CPU.h"
#include "kp_smatrix.h"
#include "kp_carma_tools.h"

#include <cuda.h>
carma_obj<float>* calculate_D_Mo(carma_context* context, int n_slopes, int n_actu_zern, bool is_zonal) {
  /*long dims[] = { 2, n_slopes, n_actu_zern }; //877 pour cas particulier en 16m
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/D_Mo_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/D_Mo_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos); */
  carma_obj<float>* D_Mo = NULL;//new carma_obj<float>(context, dims, &vec[0]);
  //delete vec_stream;
  return D_Mo;
}
carma_obj<float>* calculate_N_Act(carma_context* context, int n_actu_zern, int n_actus, bool is_zonal) {
  /*long dims[] = { 2, n_actu_zern, n_actus };
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/N_Act_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/N_Act_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos); */
  carma_obj<float>* N_Act = NULL;//new carma_obj<float>(context, dims, &vec[0]);  
  //delete vec_stream;
  return N_Act;
}
carma_obj<float>* calculate_PROJ(carma_context* context, int n_actus, int n_actu_zern, bool is_zonal) {
  /*long dims[] = { 2, n_actus, n_actu_zern };
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/PROJ_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/PROJ_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos); */
  carma_obj<float>* PROJ = NULL;//new carma_obj<float>(context, dims, &vec[0]);
  //delete vec_stream;
  return PROJ;
}
carma_obj<float>* calculate_btur(carma_context* context, int n_actu_zern, bool is_zonal) {
  /*long dims[] = { 1, n_actu_zern };
  float* zeros_tmp = new float[n_actu_zern];
  for (int i=0 ; i<n_actu_zern ; i++) zeros_tmp[i] = 0;*/
  carma_obj<float>* btur = NULL;//new carma_obj<float>(context, dims, zeros_tmp);
  //delete [] zeros_tmp ; zeros_tmp=NULL;
  return btur;
}
carma_obj<float>* calculate_SigmaV(carma_context* context, int n_actu_zern, bool is_zonal) {
  /*long dims[] = { 2, n_actu_zern, n_actu_zern };
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/SigmaV_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/SigmaV_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos);*/
 //for (int i=0 ; i<vec.size() ; i++) vec[i] = vec[i]/100; 
  carma_obj<float>* SigmaV = NULL;//new carma_obj<float>(context, dims, &vec[0]);
  //delete vec_stream;
  return SigmaV;
}
carma_obj<float>* calculate_atur(carma_context* context, int n_actu_zern, bool is_zonal) {
  /*long dims[] = { 1, n_actu_zern };
  ifstream* vec_stream;
  if (is_zonal)
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/atur_08P16_AR1Za980.dat");
  else
     vec_stream = new ifstream("/home/tgautrais/Documents/v2_kalman_kp/data/atur_08P16_AR1Mo30.dat");
  istream_iterator<float> start(*vec_stream), eos;
  vector<float> vec(start, eos); */
  carma_obj<float>* atur = NULL;//new carma_obj<float>(context, dims, &vec[0]);
  //delete vec_stream;
  return atur;
}


sutra_controller_kalman::sutra_controller_kalman(carma_context* context_, int nslope_, int nactu_) : sutra_controller(context_, nslope_, nactu_) {
   core_sparse = NULL;
   core_full = NULL;
   cusparseHandle = NULL;
   isGPU = true;
   isZonal = true;
   isSparse = true;
   isInit = false;
}



void sutra_controller_kalman::init_kalman(carma_host_obj<float>& chD_Mo, carma_host_obj<float>& chN_Act, carma_host_obj<float>& chPROJ, bool is_zonal, bool is_sparse, bool is_GPU) {
   
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

   //convert from carma_obj to kp_matrix
   kp_matrix kD_Mo, kN_Act, kPROJ;
   kp_carma_host_obj_to_kp_matrix(chD_Mo,  kD_Mo);
   kp_carma_host_obj_to_kp_matrix(chN_Act, kN_Act);
   kp_carma_host_obj_to_kp_matrix(chPROJ,  kPROJ);

   cudaSetDevice(2);

   if (is_sparse)
   {
      //convert from kp_matrix to kp_smatrix (full -> sparse)
      kp_smatrix sD_Mo, sN_Act, sPROJ;
      sD_Mo.init_from_matrix(kD_Mo);sD_Mo.resize2rowMajor();
      sN_Act.init_from_matrix(kN_Act);sN_Act.resize2rowMajor();
      sPROJ.init_from_matrix(kPROJ);sPROJ.resize2rowMajor();
      
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
         core_sparse = new kp_kalman_core_sparse_CPU(sD_Mo, sN_Act, sPROJ,
					       is_zonal);

   }  
   else if (is_GPU)
      core_full = new kp_kalman_core_full_GPU(kD_Mo, kN_Act, kPROJ, is_zonal, current_context->get_cublasHandle());
   else
      core_full = new kp_kalman_core_full_CPU(kD_Mo, kN_Act, kPROJ, is_zonal);

   cudaSetDevice(0);
   
   isInit = true;
}

sutra_controller_kalman::~sutra_controller_kalman() {
   if (cusparseHandle)
   {
      cusparseDestroy(cusparseHandle);
      cusparseHandle = NULL;
   }
}

void sutra_controller_kalman::calculate_gain(double bruit, double k_W,
    carma_host_obj<float>& chSigmaV, carma_host_obj<float>& chatur,
    carma_host_obj<float>& chbtur) {
	
   if (!isInit)
   {
      cerr << "Error | sutra_controller_kalman::calculate_gain | Kalman controller has not been initialiez"<<endl;
      exit(EXIT_FAILURE);
   }
   //convert carma_obj to kp_matrix
   kp_matrix kSigmaV;
   kp_carma_host_obj_to_kp_matrix(chSigmaV, kSigmaV);
   //convert carma_obj to kp_vector
   kp_vector katur, kbtur;
   kp_carma_host_obj_to_kp_vector(chatur, katur);
   kp_carma_host_obj_to_kp_vector(chbtur, kbtur);

   cudaSetDevice(2);
   if (core_sparse)
      core_sparse->calculate_gain(bruit, k_W, kSigmaV, katur, kbtur);
   else if (core_full)
      core_full->calculate_gain(bruit, k_W, kSigmaV, katur, kbtur);
   cudaSetDevice(0);
}
//                                                                                     
int sutra_controller_kalman::comp_com() {
   kp_vector Y_k, U_k;
   
   kp_carma_obj_to_kp_vector(*d_centroids, Y_k);
   cudaSetDevice(2);
   if (core_sparse)
   {
      core_sparse->next_step(Y_k, U_k);
   }
   else if (core_full)
   {
      core_full->next_step(Y_k, U_k);
   }
   cudaSetDevice(0);
   kp_kp_vector_to_carma_obj(U_k, *d_com);
  return -378;
}

#else

carma_obj<float>* calculate_D_Mo(carma_context* context, int n_slopes, int n_actu_zern, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_N_Act(carma_context* context, int n_actu_zern, int n_actus, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_PROJ(carma_context* context, int n_actus, int n_actu_zern, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_btur(carma_context* context, int n_actu_zern, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_SigmaV(carma_context* context, int n_actu_zern, bool is_zonal) {
  return NULL;
}
carma_obj<float>* calculate_atur(carma_context* context, int n_actu_zern, bool is_zonal) {
  return NULL;
}

sutra_controller_kalman::sutra_controller_kalman(carma_context* context_, int nslope_, int nactu_) : sutra_controller(context_, nslope_, nactu_) {
   core_sparse = NULL;
   core_full = NULL;
   cusparseHandle = NULL;
   isGPU = true;
   isZonal = true;
   isSparse = true;
   isInit = false;
}

sutra_controller_kalman::~sutra_controller_kalman() {
}

void sutra_controller_kalman::calculate_gain(double bruit, double k_W,
    carma_host_obj<float>& chSigmaV, carma_host_obj<float>& chatur,
    carma_host_obj<float>& chbtur) {
}
//
int sutra_controller_kalman::comp_com() {
  return -378;
}

#endif


//                                                                                     
