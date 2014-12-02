//kp_matrix.h
//Kalman project matrix

#ifndef __SEGER__KP_MATRIX_H__
#define __SEGER__KP_MATRIX_H__

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include "kp_KFPP.h"
#include "carma_host_obj.h"
#include "mkl.h"

using namespace std;

#include "kp_vector.h"
#include "kp_smatrix.h"
#include <typeinfo>

template<typename real> class kp_vector;
template<typename real> class kp_matrix;
template<typename real> class kp_smatrix;

template<typename real>
class kp_matrix
{
	template<typename T> friend class kp_matrix;
#ifndef KP_WITH_CARMA
 public:
   int getDim1()const{return dim1;};
   int getDim2()const{return dim2;};
   real* getData()const{return d;};
   real* getData(int i)const{return &d[i];};
   kp_matrix(const vector< vector<double> >& m);
   
   real&        operator()(int i , int j)       { return d[j * dim1 + i] ;}
   const real&  operator()(int i , int j) const { return d[j * dim1 + i] ;}
   
   //real& at(int i, int j);

   //run dsytrf/dsytri (inverse symmetric matrix);
   void inverse_symmetric();
     
	private:
   real& el(int i , int j) { return d[j * dim1 + i] ;}


#else
 public:
   int getDim1()const{return carma_host_object->getDims(1);};
   int getDim2()const{return carma_host_object->getDims(2);};
   real* getData()const{return carma_host_object->getData();};
   real* getData(int i)const{return carma_host_object->getData();};
   real& operator()(int i , int j) { return *(carma_host_object.getData(j * dim1 + i));}
   const real&  operator()(int i , int j) const { return *(carma_host_object.getData(j * dim1 + i));}

 private:
   carma_host_obj<real>* carma_host_object;
   real& el(int i , int j) { return *(carma_host_object.getData(j * dim1 + i));}
#endif

 public:
// same methods (do not depend on whether KP_WITH_CARMA is defined or not)	
   kp_matrix();
   kp_matrix(int dim1_, int dim2_);
   ~kp_matrix();   
   //analogue of mean(M,1) or mean(M) from MATLAB 
   //return vector of size dim2 contain means of each column
   void cols_mean(kp_vector<real>& rez);
   //analogue of mean(M,2) from MATLAB 
   //return vector of size dim1 contain means of each row
   void rows_mean(kp_vector<real>& rez);
   //calculate trace
   double trace();
   void init_from_matrix(const kp_matrix<real>& M, int r1, int r2, int c1, int c2);
   //we take from M only rows which in rowidx
   //analogue of MATLAB's this = M(rowidx,:)
   void init_from_rowidx(const kp_matrix<real>& M, const vector<int>& rowidx);
   //subM is matrix which was created by init_from_rowidx 
   //analogue of MATLAB's this(rowidx, : ) = subM
   void set_from_rowsubmatrix(const kp_matrix<real>& subM, const vector<int>& rowidx);
   //equivalent of MATLAB
   //find(this(:,col));
   void make_rowidx(int col, vector<int>& rowidx);
   //init from transpose matrix
   void init_from_transpose(const kp_matrix<real>& M); 
   //MATLAB code:
   //this = diag(v) * this
   //we multiplay ith row by v[i]
   void mult_each_row(const kp_vector<real>&v);
   //put submatrix at  r,c -> r + subM.dim1, c + subM.dim2 
   void set_from_submatrix(const kp_matrix<real>& subM, int r, int c);
   void init_from_smatrix(const kp_smatrix<real>& m);
   // MATLAB CODE: M = (M - repmat(v,1,M.dim2)) * val
   // So we substract v from each column of M and after multiply result by val
   void each_column_substract_mult(const kp_vector<real> &v, double val);
   //M(i,j) += v[i]
   void each_column_add(const kp_vector<real> &v);   
   void add_unit_to_diagonal();



   
// each method exists twice : one if KP_WITH_CARMA is not defined ; another one if KP_WITH_CARMA is defined
   template<typename T> kp_matrix(const kp_matrix<T>& m);
   kp_matrix(const kp_matrix<real>& m);
   void resize(int dim1_, int dim2_);
   template<typename T> void operator=( const kp_matrix<T>& m);
   void operator=( const kp_matrix<real>& m);
   void zeros();  //put all elements to 0
   //apply operation to each element
   void operator/= (real val);
   void operator*= (real val);
   void operator+= (real val);
   void operator-= (real val);
   void operator+= (const kp_matrix<real>& M);
   //run dgetrf/dgetri
   void inverse();

   
   
   
 private:
   //lda == dim1 !!!
   int dim1; //dim1 -- number of rows    (usually m)
   int dim2; //dim2 -- number of columns (usually n)
   real* d;
   void _create(int dim1_, int dim2_);
   void _clear(); 



};



//y := alpha*op(A)*x + beta * y
//op_A could be N (do nothing) or T (transpose)
template<typename real>                                                                                   
void kp_gemv(char op_A, real alpha, const kp_matrix<real>& A, const kp_vector<real>& x, real beta, kp_vector<real>& y);

template<typename real>                                                                                   
inline void kp_gemv(char op_A, int alpha, const kp_matrix<real>& A, const kp_vector<real>& x, int beta, kp_vector<real>& y)
{kp_gemv(op_A, (real) alpha, A, x, (real) beta, y);}

//y := alpha*A*x + beta * y
template<typename real>                                                                                   
inline void kp_gemv(real alpha, const kp_matrix<real>& A, const kp_vector<real>& x, real beta, kp_vector<real>& y)
{ kp_gemv('N', alpha, A, x, beta, y); }



//C := alpha*op(A)*op(B) + beta * y
//op_? could be N (do nothing) or T (transpose)
template<typename real>                                                                                   
void kp_gemm(char op_A, char op_B, real alpha, const kp_matrix<real>& A, 
	     const kp_matrix<real>& B, real beta, kp_matrix<real>& C);

template<typename real>                                                                                   
inline void kp_gemm(char op_A, char op_B, int alpha, const kp_matrix<real>& A, 
	     const kp_matrix<real>& B, int beta, kp_matrix<real>& C)
{kp_gemm(op_A, op_B, (real) alpha, A, B, (real) beta, C);}

//C := alpha*A*B + beta * y
template<typename real>                                                                                   
inline void kp_gemm(real alpha, const kp_matrix<real>& A, const kp_matrix<real>& B, real beta, kp_matrix<real>& C)
{ kp_gemm('N', 'N', alpha, A, B, beta, C);}


//analogue of mean in MATLAB 
//if dim == 1 -> calculate means of each columns 
//if dim == 2 -> calculate means of each rows
template<typename real>                                                                                   
void kp_mean(const kp_matrix<real>& M, int dim, kp_vector<real>& rez);

//eigen-value decomposition of symmetric matrix
//Q must be symmetrical matrix (we will not check it)
//Q will be replaced my matrix of eigen vectors
//DQ - vector of eigen values value
//Q_in = OQ_out * DQ * OQ_out'
template<typename real>                                                                                   
void kp_syevd(kp_matrix<real>&Q, kp_vector<real>& DQ);


template<typename real>                                                                                   
void kp_check_op_set_dim(int op, const kp_matrix<real>&M, int& dim1, int& dim2);

template<typename real>                                                                                   
double kp_calc_diff(const kp_matrix<real>& m1, const kp_matrix<real>& m2);

//make vertcat of matrises
template<typename real>                                                                                   
void kp_vertcat(kp_matrix<real>& rez, const vector<kp_matrix<real>*>& ms);
//make horizcat of matrises
template<typename real>                                                                                   
void kp_horizcat(kp_matrix<real>& rez, const vector<kp_matrix<real>*>& ms);

template<typename real>                                                                                   
ostream& operator<<(ostream& out, kp_matrix<real>& m);

template<typename real> void kp_getrf_core(int* dim1, int* dim2, real* data, int* ipiv, int* info);
template<typename real> void kp_getri_core(int* dim1, real* data, int* ipiv, real* wkopt, int* lwork, int* info);
template<typename real> void kp_gemv_core(char* op_A, const int* A_dim1, const int* A_dim2, real* alpha, real* A_data, real* x_data, int* one, real* beta, real* y_data );

template<typename real>
void kp_gemm_core(char* op_A, char* op_B, int* opA_dim1, const int* C_dim2, int* opA_dim2, real* alpha, real* A_data, const int*A_dim1, real* B_data, const int* B_dim1, real* beta, real* C_data, const int* C_dim1);

template<typename real>
void kp_syevd_core(char* jobz, char* uplo, int* n, real* Q_data, real* DQ_data, real* work, int* lwork, int* iwork, int* liwork, int* info );



#include "kp_matrix.hpp"

#endif
