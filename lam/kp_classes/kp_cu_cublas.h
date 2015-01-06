//kp_cu_cublas.h

#ifndef __SEGER__KP_CU_CUBLAS_H__
#define __SEGER__KP_CU_CUBLAS_H__

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include "kp_KFPP.h"
#include "kernels.h"
#include "kp_cu_cuda.h"
#include "kp_cu_cublas.h"
#include "cublas_v2.h"
#include "cusparse_v2.h"
#include "kp_cu_timer.h"

#include "kp_cu_vector.h"
#include "kp_cu_matrix.h"

template<typename real> class kp_cu_vector;
template<typename real> class kp_cu_matrix;
template<typename real> class kp_cu_smatrix;


// FULL MATRIX (CUBLAS)

// C = alpha*op_A(A) + beta*op_B(B); (possibilite de'avoir &C=&B ou &C=&A)
template<typename real>
void kp_cu_geam_core(cublasHandle_t& handle,
       cublasOperation_t& transa, cublasOperation_t& transb,
       int& m, int& n, const real *alpha,
       const real *A_data, int& lda, const real *beta,
       const real *B_data, int& ldb,
       real *C_data, int& ldc);
template <typename real>
void kp_cu_geam(cublasHandle_t handle,char opA, char opB, real alpha, real beta, const kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B, kp_cu_matrix<real>& C );
template <typename real>
inline void kp_cu_geam(cublasHandle_t handle,char opA, char opB, int alpha, int beta, const kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B, kp_cu_matrix<real>& C )
{kp_cu_geam(handle, opA, opB, (real) alpha, (real) beta, A, B, C);}

// C = alpha*op_A(A)*op_B(B) + beta*C; 
template <typename real>
void kp_cu_gemm_core(cublasHandle_t& handle,
                           cublasOperation_t& transa, cublasOperation_t& transb,
                           int m, int n, int k, const real  *alpha,
                           const real *A, int lda,
                           const real *B, int ldb, const real *beta,
                           real       *C, int ldc);
template <typename real>
void kp_cu_gemm(cublasHandle_t handle, char op_A, char op_B, real alpha, const kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B, real beta, kp_cu_matrix<real>& C);
template <typename real>
inline void kp_cu_gemm(cublasHandle_t handle, char op_A, char op_B, int alpha, const kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B, int beta, kp_cu_matrix<real>& C)
{kp_cu_gemm(handle, op_A, op_B, (real) alpha, A, B, (real) beta, C);}

// y = alpha*op_A(A)*x + beta*y
template <typename real>
void kp_cu_gemv_core(cublasHandle_t& handle, cublasOperation_t& trans,
                           int m, int n, const real *alpha,
                           const real *A_data, int lda,
                           const real *x_data, int incx, const real *beta,
                           real       *y_data, int incy);
template <typename real>
void kp_cu_gemv(cublasHandle_t handle, char op_A, real alpha, const kp_cu_matrix<real>& A, const kp_cu_vector<real>& x, real beta, kp_cu_vector<real>& y);
template <typename real>
void kp_cu_gemv(cublasHandle_t handle, char op_A, int alpha, const kp_cu_matrix<real>& A, const kp_cu_vector<real>& x, int beta, kp_cu_vector<real>& y)
{kp_cu_gemv(handle, op_A, (real) alpha, A, x, (real) beta, y);}

void kp_cu_op2cublastrans(int op, cublasOperation_t* trans);

template <typename real>
void kp_cu_check_op_set_dim(int op, const kp_cu_matrix<real>&M, int& dim1, int& dim2);



// SPARSE MATRIX (CUSPARSE)

// C = alpha * op_A(A) * B + beta * C ; A:sparse B,C:full
template<typename real>
void kp_cu_gemm_core(cusparseHandle_t& handle, cusparseOperation_t transA, 
    int m, int n, int k, int A_nnz, const real *alpha, const cusparseMatDescr_t& descrA, 
    const real *A_values_cu, const int *A_csrRowPtr, const int *A_colind_cu,
    const real *B_data, int ldb, const real *beta, real *C_data,  int ldc);
template <typename real>
void kp_cu_gemm(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix<real>& A, const kp_cu_matrix<real>& B, real beta, kp_cu_matrix<real>& C);
template <typename real>
inline void kp_cu_gemm(cusparseHandle_t handle, char op_A, int alpha, kp_cu_smatrix<real>& A, const kp_cu_matrix<real>& B, int beta, kp_cu_matrix<real>& C)
{kp_cu_gemm(handle, op_A, (real) alpha, A, B, (real) beta, C);}

// y = alpha * op_A(A) * x + beta * y ; A:sparse
template<typename real> 
void kp_cu_gemv_core(cusparseHandle_t& handle, cusparseOperation_t transA, 
               int m, int n, int nnz, const real *alpha, 
               const cusparseMatDescr_t& descrA, const real *A_values,  
	       const int *A_csrRowPtr, const int *A_colind_cu,
               const real *x_data, const real *beta, 
               real *y_data);

template <typename real>
void kp_cu_gemv(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix<real>& A, const kp_cu_vector<real>& x, real beta, kp_cu_vector<real>& y);
template <typename real>
inline void kp_cu_gemv(cusparseHandle_t handle, char op_A, int alpha, kp_cu_smatrix<real>& A, const kp_cu_vector<real>& x, int beta, kp_cu_vector<real>& y)
{kp_cu_gemv(handle, op_A, (real) alpha, A, x, (real) beta, y);}

// C = op_A(A) * op_B(B) ; A,B,C:sparse
template<typename real>
void kp_cu_sgemm_core(cusparseHandle_t& handle,
    cusparseOperation_t transA, cusparseOperation_t transB, int m, int n, int k, 
    const cusparseMatDescr_t& descrA, const int nnzA, const real *csrValA,
    const int *csrRowPtrA, const int *csrColIndA,
    const cusparseMatDescr_t& descrB, const int nnzB, const real *csrValB, 
    const int *csrRowPtrB, const int *csrColIndB,
    const cusparseMatDescr_t& descrC, real *csrValC,
    const int *csrRowPtrC, int *csrColIndC);
template <typename real>
void kp_cu_sgemm(cusparseHandle_t cusparseHandle, char op_A, char op_B, kp_cu_smatrix<real>& A, kp_cu_smatrix<real>& B, kp_cu_smatrix<real>& C);

template <typename real>
void kp_cu_check_op_set_dim(int op, const kp_cu_smatrix<real>&M, int& dim1, int& dim2, cusparseOperation_t* trans);

template <typename real>
void kp_cu_check_op_set_dim(int op, const kp_cu_matrix<real>&M, int& dim1, int& dim2, cusparseOperation_t* trans);

template <typename real>
void kp_cu_check_op_set_dim(int op, const kp_cu_smatrix<real>&M, int& dim1, int& dim2, cusparseOperation_t* trans, cusparseOperation_t* transinv);

#include "kp_cu_cublas.hpp"

#endif
