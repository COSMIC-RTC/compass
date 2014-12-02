//kp_cu_matrix.h
//Kalman project matrix

#ifndef __SEGER__KP_CU_MATRIX_H__
#define __SEGER__KP_CU_MATRIX_H__


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
#include <carma_obj.h>
#include "magma.h"
 
#include "kp_matrix.h"
#include "kp_smatrix.h"
#include "kp_cu_smatrix.h"
#include <typeinfo>


using namespace std;

template <typename real>
class kp_cu_matrix
{
	template<typename T> friend class kp_cu_matrix;
#ifndef KP_WITH_CARMA
 public :
	kp_cu_matrix(int dim1_, int dim2_);
	template <typename T> kp_cu_matrix(const kp_matrix<T>& m);
	kp_cu_matrix(const vector< vector<real> >& m);

	int   getDim1()const{return dim1;}
	int   getDim2()const{return dim2;}	
	real* getData()const{return d_cu;}
	real* getData(int idx)const{return &(d_cu[idx]);}

	
	void init_from_transpose(cublasHandle_t handle, const kp_cu_matrix<real>& M);
	void inverse();		
	void init_from_smatrix(const kp_cu_smatrix<real>& M);
	

 private:
	void _create(int dim1_, int dim2_);

#else
 public:
	kp_cu_matrix(carma_context* context_);
	kp_cu_matrix(carma_context* context_, int dim1_, int dim2_);
	kp_cu_matrix(carma_context* context_, const kp_matrix<real>& m);


	int   getDim1()const{return carma_object.getDims(1);}
	int   getDim2()const{return carma_object.getDims(2);}
	real* getData()const{return carma_object.getData();}
	real* getData(int idx)const{return carma_object.getData(idx);}

	void init_from_transpose(const kp_cu_matrix<real>& M);


 private:
	void _create(carma_context* context_, int dim1_, int dim2_);
	carma_obj<real>* carma_object;

#endif

 public:
// each method exists twice : one if KP_WITH_CARMA is not defined ; another one if KP_WITH_CARMA is defined
	kp_cu_matrix() ;
	template <typename T> kp_cu_matrix(const kp_cu_matrix<T>& m);
	kp_cu_matrix(const kp_cu_matrix<real>& m);
	void resize(int dim1_, int dim2_);	
	template <typename T> void operator=(const kp_cu_matrix<T>& m);
	void operator=(const kp_cu_matrix<real>& m);
	template <typename T> void operator=(const kp_matrix<T>& m);
	//if KP_WITH_CARMA is defined, handle (1st argument) is not used, this is the cublasHandle from carma_object which is used
	void gemm(cublasHandle_t handle, char op_A, char op_B, real alpha, const kp_cu_matrix<real>& A, const kp_cu_matrix<real>& B, real beta);

// same methods (do not depend on whether KP_WITH_CARMA is defined or not)	
	~kp_cu_matrix();
	void zeros(); //put all elements to 0
	void init_from_transpose(const kp_matrix<real>& M);
	void init_from_matrix(const kp_matrix<real>& M, int r1, int r2, int c1, int c2);
	void init_from_matrix(const kp_cu_matrix<real>& M, int r1, int r2, int c1, int c2);
	//apply operation to each element
	void operator/=(real val);
	void operator*=(real val);
	void operator+=(real val);
	void operator-=(real val);
	void operator+= (const kp_cu_matrix<real>& M);
	void operator-= (const kp_cu_matrix<real>& M);
	void set_from_submatrix(const kp_matrix<real>& subM, int r, int c);
	void set_from_submatrix(const kp_cu_matrix<real>& subM, int r, int c);
      
 private:
	
	void _clear() ;
	real* d_cu;
	int dim1;
	int dim2;


};

template <typename real>
void kp_cu_getrf_core(int& m, int& n, real* dA, int& ldda, int *ipiv, int *info);

template <typename real>
void kp_cu_getri_core(int& n, real* dA, int& ldda, int* ipiv, real* dwork,
       int& lwork, int* info);

template <typename real> int kp_cu_get_getri_nb_core(int& dim);

template <typename real>
void kp_cu_vertcat(kp_cu_matrix<real>& rez, const vector<kp_matrix<real>*>& ms);

template <typename real>
void kp_cu_vertcat(kp_cu_matrix<real>& rez, const vector<kp_cu_matrix<real>*>& ms);

template <typename real>
void kp_cu_horizcat(kp_cu_matrix<real>& rez, const vector<kp_matrix<real>*>& ms);

template <typename real>
void kp_cu_horizcat(kp_cu_matrix<real>& rez, const vector<kp_cu_matrix<real>*>& ms);



#include "kp_cu_matrix.hpp"

#endif
