//kp_cu_matrix.h
//Kalman project matrix

#ifndef __SEGER__KP_CU_MATRIX_H__
#define __SEGER__KP_CU_MATRIX_H__

#include "kp_matrix.h"
#include "kp_cu_smatrix.h"
#include "cublas_v2.h"
#include "cusparse_v2.h"
#include "magma.h"
#include <iostream>
 
//class kp_cu_smatrix;  //pre declaration of kp_cu_smatrix

class kp_cu_smatrix;
class kp_smatrix;

class kp_cu_matrix
{
 public:
	kp_cu_matrix() ;
	kp_cu_matrix(int dim1_, int dim2_);
	kp_cu_matrix(const kp_cu_matrix& m);
	kp_cu_matrix(const kp_matrix& m);
	kp_cu_matrix(const vector< vector<real> >& m);

	void resize(int dim1_, int dim2_);
	
	void operator=(const kp_cu_matrix& m);
	void operator=(const kp_matrix& m);
	void zeros(); //put all elements to 0

	//apply operation to each element
	void operator/=(real val);
	void operator*=(real val);
	void operator+=(real val);
	void operator-=(real val);
	void operator+= (const kp_cu_matrix& M);
	void operator-= (const kp_cu_matrix& M);

	void inverse();	
	
	//void kp_cu_alloc();
	void init_from_matrix(const kp_matrix& M, int r1, int r2, int c1, int c2);
	void init_from_matrix(const kp_cu_matrix& M, int r1, int r2, int c1, int c2);
	void init_from_transpose(const kp_matrix& M);
	void init_from_transpose(cublasHandle_t handle, const kp_cu_matrix& M);
	void init_from_smatrix(const kp_cu_smatrix& M);

	void set_from_submatrix(const kp_matrix& subM, int r, int c);
	void set_from_submatrix(const kp_cu_matrix& subM, int r, int c);

	char precision;
	
	~kp_cu_matrix();


	
   
      
 private:
	
	void _clear() ;
	void _create(int dim1_, int dim2_);

	 
 public:
	real* d_cu;
	int dim1;
	int dim2;


};

// Calcule C = alpha*op(A) + beta*op(B) (possibilite de'avoir &C=&B ou &C=&A)
void kp_cu_geam(cublasHandle_t hamdle,char opA, char opB, real alpha, real beta, const kp_cu_matrix& A, const kp_cu_matrix& B, kp_cu_matrix& C );

void kp_cu_vertcat(kp_cu_matrix& rez, const vector<kp_matrix*>& ms);
void kp_cu_vertcat(kp_cu_matrix& rez, const vector<kp_cu_matrix*>& ms);
void kp_cu_horizcat(kp_cu_matrix& rez, const vector<kp_matrix*>& ms);
void kp_cu_horizcat(kp_cu_matrix& rez, const vector<kp_cu_matrix*>& ms);


void kp_cu_check_op_set_dim(int op, const kp_cu_matrix&M, int& dim1, int& dim2, cublasOperation_t* trans);

void kp_cu_gemm(cublasHandle_t handle, char op_A, char op_B, real alpha, const kp_cu_matrix& A, const kp_cu_matrix& B, real beta, kp_cu_matrix& C);

void kp_cu_gemv(cublasHandle_t handle, char op_A, real alpha, const kp_cu_matrix& A, const kp_cu_vector& x, real beta, kp_cu_vector& y);

#endif
