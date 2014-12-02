//kp_cu_smatrix.h
//kalman project sparse matrix
//we save matrix in coordinate format from MKL (see Sparse Matrix Storage Formats in MKL)
//we use one-based indexing (because one-based indexing would need row first packing for dens arrays in kp_gemm)

#ifndef __SEGER__KP_CU_SMATRIX_H__
#define __SEGER__KP_CU_SMATRIX_H__

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
#include <carma_obj.h>

#include "kp_smatrix.h"
#include "kp_cu_vector.h"
#include "kp_cu_matrix.h"
#include <typeinfo>

template<class real> class kp_cu_matrix;
template<class real> class kp_cu_vector;


using namespace std;

template <typename real>
class kp_cu_smatrix
{
	template<typename T> friend class kp_cu_smatrix;

 public:
	 kp_cu_smatrix();
         template <typename T> kp_cu_smatrix(const kp_cu_smatrix<T>& M);
         kp_cu_smatrix(const kp_cu_smatrix<real>& M);
	 template <typename T> kp_cu_smatrix(const kp_smatrix<T>& M);
	 ~kp_cu_smatrix();

	 template <typename T> void operator=(const kp_cu_smatrix<T>& M);
	 void operator=(const kp_cu_smatrix<real>& M);
         template <typename T> void operator=(const kp_smatrix<T>& M);
	 

	 void resize(int nnz_, int dim1_, int dim2_);
	
         void init_from_transpose(const kp_cu_smatrix<real>& M);

         bool isColumnMajor();
         char get_majorDim()const{return majorDim;}
         void set_majorDim(char c){majorDim=c;}
	 void convert2csr(cusparseHandle_t handle);
	 void convert2csrT(cusparseHandle_t handle);
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
	int getDim1()const{return dim1;};
	int getDim2()const{return dim2;};


	// ONE-BASED INDEXING
	real* values_cu;
	int*  rowind_cu;
	int*  colind_cu;
	cusparseMatDescr_t descr;

	template <typename T> friend void kp_cu_gemm(cusparseHandle_t handle, char op_A, T alpha, kp_cu_smatrix<T>& A, const kp_cu_matrix<T>& B, T beta, kp_cu_matrix<T>& C);
	//template <typename T> friend void kp_cu_gemm(cusparseHandle_t handle,cublasHandle_t cublasHandle, char op_A, char op_B, T alpha, kp_cu_smatrix& A, kp_cu_matrix& B, T beta, kp_cu_matrix& C);
	template <typename T> friend void kp_cu_sgemm(cusparseHandle_t cusparseHandle, char op_A, char op_B, kp_cu_smatrix<T>& A, kp_cu_smatrix<T>& B, kp_cu_smatrix<T>& C);
	template <typename T> friend void kp_cu_gemv(cusparseHandle_t handle, char op_A, T alpha, kp_cu_smatrix<T>& A, const kp_cu_vector<T>& x, T beta, kp_cu_vector<T>& y);

	friend void kp_cu_cudaFree(int*);
	friend void kp_cu_cudaMalloc(int*);

	int   nnz;
	int   dim1; //nombre de lignes
	int   dim2; //nombre de colones

	
};
#include "kp_cu_smatrix.hpp"

#endif
