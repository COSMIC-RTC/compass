//kp_smatrix.h
//kalman project sparse matrix
//we save matrix in coordinate format from MKL (see Sparse Matrix Storage Formats in MKL)
//we use one-based indexing (because one-based indexing would need row first packing for dens arrays in kp_gemm)

#ifndef __SEGER__KP_SMATRIX_H__
#define __SEGER__KP_SMATRIX_H__

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include "kp_KFPP.h"

#include <algorithm>

#include "kp_vector.h"
#include "kp_matrix.h"
#include "mkl.h"
#include <typeinfo>

template<class real> class kp_matrix;
template<class real> class kp_vector;

using namespace std;


template <typename real>
class kp_smatrix
{
	template<typename T> friend class kp_smatrix;
 public:
   kp_smatrix();

   template <typename T> kp_smatrix(const kp_smatrix<T>& sm);
   kp_smatrix(const kp_smatrix<real>& sm);
   ~kp_smatrix();
   
   //delete all arrays and create new for nnz=new_nnz
   void resize(int new_nnz, int dim1_, int dim2_); 
   
   template <typename T> void operator=(const kp_smatrix<T>& M);
   void operator=(const kp_smatrix<real>& M);
   
   //we take from M only rows which in rowidx
   //analogue of MATLAB's this = M(idx,:)
   //idx is in normal zero based indexing
   //but we keep sparse matrix in one-based indexing
   void init_from_rowidx(const kp_smatrix<real>& M, const vector<int>& idx);
   
   //init from transpose sparce matrix
   void init_from_transpose(const kp_smatrix<real>& M);
   void check();

   void init_from_matrix(const kp_matrix<real>& B, double epsilon = 0.0);

   void resize2rowMajor();
   void resize2colMajor();
   char get_majorDim()const  {return majorDim;}
   void set_majorDim(char c){majorDim=c;}
   
 private:
   void _create(int nnz_, int dim1_, int dim2_); //create new arrays
   void _clear(); //clear arrays
   
   char majorDim; //U - undefined
                  //R - row major
                  //C - col major
   
 public:

   int getDim1()const{return dim1;};
   int getDim2()const{return dim2;};

   //dim1 -- number of rows    (usually m)
   //dim2 -- number of columns (usually n)
   int dim1, dim2;
   
   //see Sparse Matrix Storage Formats (coordinate format)
   //ONE-BASED!!!!
   real* values;
   int*  rowind;
   int*  colind;
   int   nnz; //size of values, rowind and colind
};

//Multiply sparce matrix by vector
//y := alpha*A*x + beta * y
template <typename real>
void kp_gemv(real alpha, kp_smatrix<real>& A, kp_vector<real>& x, real beta, kp_vector<real>& y);
template <typename real>
inline void kp_gemv(int alpha, kp_smatrix<real>& A, kp_vector<real>& x, int beta, kp_vector<real>& y)
{kp_gemv((real) alpha, A, x, (real) beta, y);}

//Multiply sparce matrix by dense matrix
//C := alpha*op(A)*B + beta * y
//op_A could be N (do nothing) or T (transpose)
template <typename real>
void kp_gemm(char op_A, real alpha, kp_smatrix<real>& A, kp_matrix<real>& B, real beta, kp_matrix<real>& C);
template <typename real>
inline void kp_gemm(char op_A, int alpha, kp_smatrix<real>& A, kp_matrix<real>& B, int beta, kp_matrix<real>& C)
{kp_gemm(op_A, (real) alpha, A, B, (real) beta, C);}

//Multiply sparce matrix by dense matrix
//C := alpha*A*B + beta * y
template <typename real>
inline void kp_gemm(real alpha, kp_smatrix<real>& A, kp_matrix<real>& B, real beta, kp_matrix<real>& C)
{kp_gemm('N', alpha, A, B, beta, C);}

template <typename real>
void kp_check_op_set_dim(int op, const kp_smatrix<real>&M, int& dim1, int& dim2);

template<typename real>
void kp_gemv_core(char *transa, MKL_INT *m, MKL_INT *k, real *alpha, char *matdescra, real *A_values, MKL_INT *A_rowind, MKL_INT *A_colind, MKL_INT *A_nnz, real *x_data, real *beta, real *y_data);

template<typename real>
void kp_gemm_core(char *transa, MKL_INT *m, MKL_INT *n, MKL_INT *k, real *alpha, char *matdescra, real *A_values, MKL_INT *A_rowind, MKL_INT *A_colind, MKL_INT *A_nnz, real *B_data, MKL_INT *ldb, real *beta, real *C_data, MKL_INT *ldc);


#include "kp_smatrix.hpp"


#endif
