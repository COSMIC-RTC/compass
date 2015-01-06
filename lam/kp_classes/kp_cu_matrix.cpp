#include "kp_cu_matrix.h"



template <>
void kp_cu_getrf_core<double>(int& m, int& n,
		double* dA, int& ldda,
		int *ipiv, int *info)
{
	magma_int_t status;

	status = magma_dgetrf_gpu(m, n, dA, ldda, ipiv, info);

	if (status != 0)
	{	
		cerr<<"Erreur | kp_cu_matrix::inverse | Erreur lors de l'execution de magma_dgetrf_gpu"<<endl;
		//exit(EXIT_FAILURE);
	}

}
template <>
void kp_cu_getrf_core<float>(int& m, int& n,
		float* dA, int& ldda,
		int *ipiv, int *info )
{
	magma_int_t status;

	status = magma_sgetrf_gpu(m, n, dA, ldda, ipiv, info);

	if (status != 0)
	{	
		cerr<<"Erreur | kp_cu_matrix::inverse | Erreur lors de l'execution de magma_sgetrf_gpu"<<endl;
		//exit(EXIT_FAILURE);
	}
}


template <>
void kp_cu_getri_core<double>(int& n, double* dA,
		int& ldda, int* ipiv, double* dwork,
		int& lwork, int* info)
{
  int status = magma_dgetri_gpu(n, dA, ldda, ipiv, dwork, lwork, info);

   if (status != 0)
   {
      cerr<<"Erreur | kp_cu_matrix::inverse | Erreur lors de l'execution de magma_sgetri_gpu"<<endl;
      //exit(EXIT_FAILURE);
   }
}
template <>
void kp_cu_getri_core<float>(int& n, float* dA,
		int& ldda, int* ipiv, float* dwork,
		int& lwork, int* info)
{
   int status = magma_sgetri_gpu(n, dA, ldda, ipiv, dwork, lwork, info);

   if (status != 0)
   {
      cerr<<"Erreur | kp_cu_matrix::inverse | Erreur lors de l'execution de magma_sgetri_gpu"<<endl;
      //exit(EXIT_FAILURE);
   }
}


template <>
int kp_cu_get_getri_nb_core<double>(int& dim)
{
	return magma_get_dgetri_nb(dim);
}
template <>
int kp_cu_get_getri_nb_core<float>(int& dim)
{
	return magma_get_sgetri_nb(dim);
}


