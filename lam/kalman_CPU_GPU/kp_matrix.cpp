#include "kp_matrix.h"

template<>                                                                                   
void kp_getrf_core<double>(int* dim1, int* dim2, double* data, int* ipiv, int* info)
{
   dgetrf(dim1, dim2, data, dim1, ipiv, info );
}
template<>                                                                                   
void kp_getrf_core<float>(int* dim1, int* dim2, float* data, int* ipiv, int* info)
{
   sgetrf(dim1, dim2, data, dim1, ipiv, info );
}

                                                    
template<>                                                                                   
void kp_getri_core<double>(int* dim1, double* data, int* ipiv, double* wkopt, int* lwork, int* info)
{
   dgetri(dim1, data, dim1, ipiv, wkopt, lwork, info );
}
template<>                                                                                   
void kp_getri_core<float>(int* dim1, float* data, int* ipiv, float* wkopt, int* lwork, int* info)
{
   sgetri(dim1, data, dim1, ipiv, wkopt, lwork, info );
}


template<> 
void kp_gemv_core<double>(char* op_A, const int* A_dim1, const int* A_dim2, double* alpha, double* A_data, double* x_data, int* one, double* beta, double* y_data )
{
   dgemv(op_A, A_dim1, A_dim2, alpha, A_data, A_dim1, x_data, one, beta, y_data, one);
}

template<> 
void kp_gemv_core<float>(char* op_A, const int* A_dim1, const int* A_dim2, float* alpha, float* A_data, float* x_data, int* one, float* beta, float* y_data )
{
   sgemv(op_A, A_dim1, A_dim2, alpha, A_data, A_dim1, x_data, one, beta, y_data, one);
}


template<>
void kp_gemm_core<double>(char* op_A, char* op_B, int* opA_dim1, const int* C_dim2, int* opA_dim2, double* alpha, double* A_data, const int*A_dim1, double* B_data, const int* B_dim1, double* beta, double* C_data, const int* C_dim1)
{
   dgemm(op_A, op_B, opA_dim1, C_dim2,  opA_dim2, alpha, A_data, A_dim1,
	 B_data, B_dim1,  beta, C_data, C_dim1);
}

template<>
void kp_gemm_core<float>(char* op_A, char* op_B, int* opA_dim1, const int* C_dim2, int* opA_dim2, float* alpha, float* A_data, const int*A_dim1, float* B_data, const int* B_dim1, float* beta, float* C_data, const int* C_dim1)
{
   sgemm(op_A, op_B, opA_dim1, C_dim2,  opA_dim2, alpha, A_data, A_dim1,
	 B_data, B_dim1,  beta, C_data, C_dim1);
}


template<>
void kp_syevd_core<double>(char* jobz, char* uplo, int* n, double* Q_data, double* DQ_data, double* work, int* lwork, int* iwork, int* liwork, int* info )
{
   DSYEVD(jobz, uplo, n, Q_data, n, DQ_data, work, lwork, iwork, liwork, info);   
}
template<>
void kp_syevd_core<float>(char* jobz, char* uplo, int* n, float* Q_data, float* DQ_data, float* work, int* lwork, int* iwork, int* liwork, int* info )
{
   SSYEVD(jobz, uplo, n, Q_data, n, DQ_data, work, lwork, iwork, liwork, info);   
}

