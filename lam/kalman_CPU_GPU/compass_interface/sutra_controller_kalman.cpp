//sutra_controller_kalman.cpp

#include "sutra_controller_kalman.h"
#include "kp_carma_tools.h"


//nslope = dim1(D_Mo)
//nactu  = dim1(N_Act)
sutra_controller_kalman::sutra_controller_kalman(carma_context* context,
						 carma_obj<float>& cD_Mo,
						 carma_obj<float>& cN_Act,
						 carma_obj<float>& cPROJ,
						 bool is_zonal)
:sutra_controller(context, kp_carma_obj_get_dim1(cD_Mo), kp_carma_obj_get_dim1(cN_Act))
{
   //convert from carma_obj to kp_matrix
   kp_matrix kD_Mo, kN_Act, kPROJ;
   kp_carma_obj_to_kp_matrix(cD_Mo,  kD_Mo);
   kp_carma_obj_to_kp_matrix(cN_Act, kN_Act);
   kp_carma_obj_to_kp_matrix(cPROJ,  kPROJ);

   //convert from kp_matrix to kp_smatrix (full -> sparce)
   kp_smatrix sD_Mo, sN_Act, sPROJ;
   sD_Mo.init_from_matrix(kD_Mo);
   sN_Act.init_from_matrix(kN_Act);
   sPROJ.init_from_matrix(kPROJ);
   
   
   cusparseStatus_t cusparseStat = cusparseCreate(&cusparseHandle);
   if (cusparseStat != CUSPARSE_STATUS_SUCCESS)
     {
	cerr<<"Error | sutra_controller_kalman::sutra_controller_kalman  | cusparseCreate failed "<<endl;
	exit(EXIT_FAILURE);
     }
   
   core_sparce = new kp_kalman_core_sparse_GPU(sD_Mo, 
					       sN_Act, 
					       sPROJ,
					       is_zonal, 
					       context->get_cublasHandle(), 
					       cusparseHandle);
}
//                                                                                     
void sutra_controller_kalman::calculate_gain(double bruit,
					     double k_W,
					     carma_obj<float>& cSigmaV,
					     carma_obj<float>& catur,
					     carma_obj<float>& cbtur)
{
   //convert carma_obj to kp_matrix
   kp_matrix kSigmaV;
   kp_carma_obj_to_kp_matrix(cSigmaV, kSigmaV);
   
   //convert carma_obj to kp_vector
   kp_vector katur, kbtur;
   kp_carma_obj_to_kp_vector(catur, katur);
   kp_carma_obj_to_kp_vector(cbtur, kbtur);
   
   core_sparce->calculate_gain(bruit, k_W, kSigmaV, katur, kbtur);
}
//                                                                                     
int sutra_controller_kalman::comp_com()
{
   kp_vector Y_k, U_k;
   kp_carma_obj_to_kp_vector(*d_centroids, Y_k);
   
   core_sparce->next_step(Y_k, U_k);
   
   kp_kp_vector_to_carma_obj(U_k, *d_com);
   //we don't know meaning of this int
   return -378;
}
//                                                                                     
