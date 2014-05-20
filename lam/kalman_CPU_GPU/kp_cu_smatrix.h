//kp_cu_smatrix.h
//kalman project sparse matrix
//we save matrix in coordinate format from MKL (see Sparse Matrix Storage Formats in MKL)
//we use one-based indexing (because one-based indexing would need row first packing for dens arrays in kp_gemm)

#ifndef __SEGER__KP_CU_SMATRIX_H__
#define __SEGER__KP_CU_SMATRIX_H__

#include "kp_cu_vector.h"
#include "kp_cu_matrix.h"
#include "kp_cu_smatrix.h"
#include "kp_smatrix.h"
#include "cusparse_v2.h"
#include "kp_cu_timer.h"

//class kp_smatrix;

class kp_cu_smatrix
{
 public:
	 kp_cu_smatrix();
	 kp_cu_smatrix(const kp_cu_smatrix& M);
	 kp_cu_smatrix(const kp_smatrix& M);
	 ~kp_cu_smatrix();

	 void operator=(const kp_cu_smatrix& M);
	void operator=(const kp_smatrix& M);
	 

	void resize(int nnz_, int dim1_, int dim2_);
	void init_from_transpose(const kp_cu_smatrix& M);
	bool isColumnMajor();
char get_majorDim()const{return majorDim;}
void set_majorDim(char c){majorDim=c;}
 private:
	void _create(int nnz_, int dim1_, int dim2_);
 	void _clear();

	char majorDim;
	bool isCSRconverted;
	bool isCSRconvertedT;
 public:
	int* csrRowPtr;
	int* csrRowPtrT;
 public:
	// ONE-BASED INDEXING
	real* values_cu;
	int*  rowind_cu;
	int*  colind_cu;
	cusparseMatDescr_t descr;
	friend void kp_cu_gemm(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix& A, const kp_cu_matrix& B, real beta, kp_cu_matrix& C);
	//friend void kp_cu_gemm(cusparseHandle_t handle,cublasHandle_t cublasHandle, char op_A, char op_B, real alpha, kp_cu_smatrix& A, kp_cu_matrix& B, real beta, kp_cu_matrix& C);

	friend void kp_cu_sgemm(cusparseHandle_t cusparseHandle, char op_A, char op_B, kp_cu_smatrix& A, kp_cu_smatrix& B, kp_cu_smatrix& C);
	friend void kp_cu_gemv(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix& A, const kp_cu_vector& x, real beta, kp_cu_vector& y);
	friend void kp_cu_cudaFree(int*);
	friend void kp_cu_cudaMalloc(int*);

	int   nnz;
	int   dim1; //nombre de lignes
	int   dim2; //nombre de colones

	
};

void kp_cu_check_op_set_dim(int op, const kp_cu_smatrix&M, int& dim1, int& dim2, cusparseOperation_t* trans);
void kp_cu_check_op_set_dim(int op, const kp_cu_matrix&M, int& dim1, int& dim2, cusparseOperation_t* trans);
void kp_cu_check_op_set_dim(int op, const kp_cu_smatrix&M, int& dim1, int& dim2, cusparseOperation_t* trans, cusparseOperation_t* transinv);

// C = alpha * op_A(A) * B + beta * C 
void kp_cu_gemm(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix& A, const kp_cu_matrix& B, real beta, kp_cu_matrix& C);

// C = alpha * A * B + beta * C 
inline void kp_cu_gemm(cusparseHandle_t handle,cublasHandle_t cublashandle, real alpha, kp_cu_smatrix& A, const kp_cu_matrix& B, real beta, kp_cu_matrix& C)
{kp_cu_gemm(handle, 'N', alpha, A, B, beta, C);}

// y = alpha * op_A(A) * x + beta * y
void kp_cu_gemv(cusparseHandle_t handle, char op_A, real alpha, kp_cu_smatrix& A, const kp_cu_vector& x, real beta, kp_cu_vector& y);

// C = op_A(A) * op_B(B) 
void kp_cu_sgemm(cusparseHandle_t cusparseHandle, char op_A, char op_B, kp_cu_smatrix& A, kp_cu_smatrix& B, kp_cu_smatrix& C);
#endif
