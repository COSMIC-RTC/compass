//kp_matrix.cpp

#include "kp_matrix.h"
#include "kp_smatrix.h"
#include "mkl_blas.h"
#include "mkl.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <cmath>
using namespace std;

#ifdef __WITH_FLOPS_CALCULATOR__
#include "kp_flopscalc.h"
#endif

kp_matrix::kp_matrix()
{
   _create(0,0);
}
//                                                                                            
kp_matrix::kp_matrix(int dim1_, int dim2_)
{
   _create(dim1_, dim2_);
}
//                                                                                           
kp_matrix::kp_matrix(const kp_matrix& m)
{
   _create(m.dim1, m.dim2);
   memcpy(d, m.d, sizeof(real) * dim1 * dim2);     
}
//                                                                                           
kp_matrix::kp_matrix(const vector< vector<double> >& m)
{
   size_t d1, d2;   
   if (m.size() == 0)
     {
	d1 = 0;
	d2 = 0;
     }
   else
     {
	d1 = m.size();
	d2 = m[0].size();
	for (size_t i = 0 ; i < m.size(); i++)
	  if (d2 != m[i].size())
	    {
	       cout<<"error | kp_matrix::kp_matrix | vector< vector<double> > is not matrix"<<endl;
	       exit(EXIT_FAILURE);
	    }	
     }
   _create(d1,d2);
   for (size_t i = 0 ; i < m.size() ; i++)
     for (size_t j = 0 ; j < m[i].size() ; j++)
       el(i,j) = m[i][j];
}
//                                                                                           
kp_matrix::~kp_matrix()
{
   _clear();
}
//                                                                                           
void kp_matrix::resize(int dim1_, int dim2_)
{
   if (dim1 * dim2 != dim1_ * dim2_)
     {
	_clear();
	_create(dim1_, dim2_);
     }
   else
     {
	dim1 = dim1_;
	dim2 = dim2_;
     }
}
//                                                                                           
real& kp_matrix::at(int i, int j)
{
   if (i < 0 || i >= dim1 || j < 0 || j >= dim2)
     {
	cerr<<"error | kp_matrix::at | index problem"<<endl;
	exit(EXIT_FAILURE);
     }
   return this->operator()(i,j);
}
//                                                                                           
void kp_matrix::operator=(const kp_matrix&m)
{
   resize(m.dim1, m.dim2);
   memcpy(d, m.d, sizeof(real) * dim1 * dim2);     
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif
}
//                                                                                           
void kp_matrix::zeros()
{
   memset(d, 0, sizeof(real) * dim1 * dim2);
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif
}
//                                                                                           
void kp_matrix::operator/=(real val)
{
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     d[i] /= val;
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif
}
//                                                                                           
void kp_matrix::operator*=(real val)
{
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     d[i] *= val;
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif
}
//                                                                                           
void kp_matrix::operator+=(real val)
{
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     d[i] += val;
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif   
}
//                                                                                           
void kp_matrix::operator-=(real val)
{
   int s = dim1 * dim2;
   for (int i = 0 ; i < s ; i++)
     d[i] -= val;

   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif   
}
//                                                                                           
void kp_matrix::operator+= (const kp_matrix& M)
{
   if (dim1 != M.dim1 || dim2 != M.dim2)
     {
	cerr<<"error | kp_matrix::operator+= | dimension problems"<<endl;
     }
   for (int i = 0 ; i < dim1 * dim2 ; i++)
     d[i] += M.d[i];

   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif   
}
//                                                                                           
void kp_matrix::cols_mean(kp_vector& rez)
{
   rez.resize(dim2);
   rez.zeros();
   if (dim1 == 0 || dim2 == 0)
     return;
   for (size_t j = 0 ; j < (unsigned int)dim2 ; j++) 
     for (size_t i = 0 ; i < (unsigned int)dim1 ; i++)
       {
	  rez[j] += el(i, j);
       }
   rez /= (real) dim1;

   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif   
}
//                                                                                           
void kp_matrix::rows_mean(kp_vector& rez)
{
   rez.resize(dim1);
   rez.zeros();
   if (dim1 == 0 || dim2 == 0)
     return;
   for (size_t j = 0 ; j < (unsigned int)dim2 ; j++) 
     for (size_t i = 0 ; i < (unsigned int)dim1 ; i++)
       {
	  rez[i] += el(i, j);
       }
   rez /= (real) dim2;
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif   
}
//                                                                                           
double kp_matrix::trace()
{
   if (dim1 != dim2)
     {
	cout<<"error | kp_matrix::trace | dim1 != dim2 (matrix must be square)"<<endl;
	exit(EXIT_FAILURE);
     }
   double sum = 0;
   for (size_t i = 0 ; i < (unsigned int)dim1 ; i++)
     sum += el(i,i);
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other(dim1);
   #endif   
   return sum;
}
//                                                                                           
void kp_matrix::init_from_matrix(const kp_matrix& M, int r1, int r2, int c1, int c2)
{
   if (r1 < 0 || r1 > r2 || r2 > M.dim1 || 
       c1 < 0 || c1 > c2 || c2 > M.dim2)
     {
	cout<<r1<<" "<<r2<<" "<<M.dim1<<endl;
	cout<<c1<<" "<<c2<<" "<<M.dim2<<endl;
	cerr<<"Error | kp_matrix::init_from_matrix | index problems"<<endl;
	exit(EXIT_FAILURE);
     }
   resize(r2 - r1, c2 - c1);
   
   for (int i = 0 ; i < r2 - r1 ; i++)
     for (int j = 0 ; j < c2 - c1 ; j++)
       el(i, j) = M(i + r1, j + c1); 
   
   #ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other((long)(r2 - r1) * (long)(c2 - c1));
   #endif
}
//                                                                                           
void kp_matrix::init_from_rowidx(const kp_matrix& M, const vector<int>& rowidx)
{
   if (this == &M)
     {
	cerr<<"Error | kp_matrix::init_from_rowidx | the same matrix"<<endl;
	exit(EXIT_FAILURE);
     }
   resize(rowidx.size(), M.dim2);
   for (unsigned int i = 0 ; i < rowidx.size(); i++)
     {
	if (rowidx[i] < 0 || rowidx[i] >= M.dim1)
	  {
	     cerr<<"error | kp_matrix::init_from_rowidx | index error "<<endl;
	  }
     }
   for (int j = 0 ; j < M.dim2 ; j++)
     for (unsigned int i = 0 ; i < rowidx.size(); i++)
       {
	  el(i,j) = M(rowidx[i], j);
       }
#ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other((long)M.dim2 * rowidx.size());
#endif
}
//                                                                                           
void kp_matrix::set_from_rowsubmatrix(const kp_matrix& subM, const vector<int>& rowidx)
{
   if (this == &subM)
     {
	cerr<<"Error | kp_matrix::kp_matrix::set_from_rowsubmatrix | the same matrix"<<endl;
	exit(EXIT_FAILURE);
     }
   if ((unsigned int)subM.dim1 != rowidx.size() || subM.dim2 != dim2)
     {
	cerr<<"Error | kp_matrix::set_from_rowsubmatrix | dimension probles"<<endl;
	exit(EXIT_FAILURE);
     }
   for (unsigned int i = 0 ; i < rowidx.size(); i++)
     {
	if (rowidx[i] < 0 || rowidx[i] >= dim1)
	  {
	     cerr<<"error | kp_matrix::set_from_rowsubmatrix | index error "<<endl;
	  }
     }
   for (int j = 0 ; j < dim2 ; j++)
     for (unsigned int i = 0 ; i < rowidx.size() ; i++)
       {
	  el(rowidx[i], j) = subM(i,j);	  
       }
#ifdef __WITH_FLOPS_CALCULATOR__
   kp_fpc_global.add_other((long)dim2 * rowidx.size());
#endif

}
//                                                                                           
void kp_matrix::make_rowidx(int col, vector<int>& rowidx)
{
   if (col < 0 || col >= dim2)
     {
	cerr<<"Error | kp_matrix::kp_make_rowidx | index problem"<<endl;
	exit(EXIT_FAILURE);
     }
   rowidx.clear();
   for (int i = 0; i < dim1 ; i++)
     {
	if (el(i, col) != 0)
	  {
	     if (fabs(el(i, col)) < 0.1)
	       {
		  cerr<<"Error | kp_matrix::kp_make_rowidx | looks like strange"<<endl;
		  exit(EXIT_FAILURE);
	       }
	     rowidx.push_back(i);
	  }
     }
   
}
//                                                                                           
void kp_matrix::init_from_transpose(const kp_matrix& M)
{
   resize(M.dim2, M.dim1);
   for (int i = 0 ; i < dim1; i++)
     for (int j = 0 ; j < dim2 ; j++)
       el(i,j) = M(j,i);

   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif   
}
//                                                                                           
void kp_matrix::mult_each_row(const kp_vector& v)
{
   if (v.size() != dim1)
     {
	cerr<<"Error | kp_matrix::mult_each_row | dimension problem"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int j = 0 ; j < dim2 ; j++)
     for (int i = 0 ; i < dim1 ; i++)
       el(i,j) *= v[i];
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif   
}
//                                                                                           
void kp_matrix::set_from_submatrix(const kp_matrix& subM, int r, int c)
{
   if (subM.dim1 + r > dim1 || subM.dim2 + c > dim2 || r < 0 || c < 0)
     {
	cerr<<"Error | kp_matrix::set_from_submatrix | dimension problem"<<endl;
	cerr<<dim1<<"x"<<dim2<<" "<<subM.dim1<<"x"<<subM.dim2<<" "<<r<<" "<<c<<endl;
	exit(EXIT_FAILURE);
     }
   for (int j = 0 ; j < subM.dim2 ; j++)
     for (int i = 0 ; i < subM.dim1 ; i++)
       {
	  el(i + r, j + c) = subM(i,j);
       }
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)subM.dim1 * (long)subM.dim2);
   #endif   
}
//                                                                                           
void kp_matrix::init_from_smatrix(const kp_smatrix&m)
{
   resize(m.dim1, m.dim2);
   zeros();
   for (int i = 0 ; i < m.nnz ; i++)
     {
	if (m.rowind[i] < 1 || m.rowind[i] > dim1 || m.colind[i] < 1 || m.colind[i] > dim2)
	  {
	     cerr<<"error | kp_matrix::init_from_smatrix | smatrix problem"<<endl;
	     exit(EXIT_FAILURE);
	  }
	el(m.rowind[i] - 1 , m.colind[i] - 1) = m.values[i];
     }   
}
//                                                                                           
void kp_matrix::each_column_substract_mult(const kp_vector &v, double val)
{
   kp_matrix& M = *this;
   
   if (v.size() != M.dim1)
     {
	cerr<<"Error | kp_matrix::each_column_substract_mult | dimension problems"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int j = 0 ; j < M.dim2 ; j++)
     for (int i = 0 ; i < M.dim1 ; i++)
       {
	  M(i,j) -= v[i];
	  M(i,j) *= val;
       }
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2 * 2);
   #endif   

}
//                                                                                           
void kp_matrix::each_column_add(const kp_vector &v)
{
   kp_matrix& M = *this;
   if (v.size() != M.dim1)
     {
	cerr<<"Error | kp_matrix::each_column_add | dimension problems"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int j = 0 ; j < M.dim2 ; j++)
     for (int i = 0 ; i < M.dim1 ; i++)
       {
	  M(i,j) += v[i];
       }
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other((long)dim1 * (long)dim2);
   #endif   
}
//                                                                                           
void kp_matrix::add_unit_to_diagonal()
{
   kp_matrix& M = *this;
   if (M.dim1 != M.dim2)
     {
	cerr<<"Error | kp_matrix::add_unit_to_diagonal | matrix isn't square"<<endl;
	exit(EXIT_FAILURE);
     }
   for (int i = 0 ; i < M.dim1 ; i++)
     M(i,i) += 1;
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_other(dim1);
   #endif   
}
//                                                                                           
void kp_matrix::inverse()
{   
   if (dim1 != dim2)
     {
	cerr<<"Error | kp_matrix::inverse | matrix is not square"<<endl;
	exit(EXIT_FAILURE);
     }
   int info;
   int* ipiv = new int[max(1,min(dim1,dim2))];
#ifdef KP_SINGLE   
   sgetrf(&dim1, &dim2, d, &dim1, ipiv, &info );
#else
   dgetrf(&dim1, &dim2, d, &dim1, ipiv, &info );
#endif
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse | some error in dgetrf "<<info<<endl;
	exit(EXIT_FAILURE);
     }

      
   real wkopt;
   int lwork = -1;
   
#ifdef KP_SINGLE   
   sgetri(&dim1, d, &dim1, ipiv, &wkopt, &lwork, &info );
#else
   dgetri(&dim1, d, &dim1, ipiv, &wkopt, &lwork, &info );
#endif
   lwork = (int)wkopt;
   real *work = new real[lwork];
   
#ifdef KP_SINGLE   
   sgetri(&dim1, d, &dim1, ipiv, work, &lwork, &info );
#else
   dgetri(&dim1, d, &dim1, ipiv, work, &lwork, &info );
#endif
   
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse | some error in dgetri "<<info<<endl;
	exit(EXIT_FAILURE);
     }
   delete [] ipiv;
   delete [] work;
   
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   //for dgetrf 2/3 * n^3 (MKL documentation)
   //for dgetri 4/3 * n^3 (MKL documentation)
   kp_fpc_global.add_inverse((long)dim1 * (long)dim1 * (long)dim1 * 2);
   #endif   
   
}
//                                                                                           
void kp_matrix::inverse_symmetric()
{   
   inverse();
/*      if (dim1 != dim2)
     {
	cerr<<"Error | kp_matrix::inverse_symmetric | matrix is not square"<<endl;
	exit(EXIT_FAILURE);
     }
   int* ipiv = new int[max(1,dim1)];
   double wkopt;
   int lwork = -1;
   int info;
   dsytrf( "L", &dim1, d, &dim1, ipiv, &wkopt, &lwork, &info );
   
   lwork = max((int)wkopt, 2 * dim1);
   double *work = new double[lwork];
   
   dsytrf( "L", &dim1, d, &dim1, ipiv, work, &lwork, &info );
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse_symmetric | some error in dsytrf "<<info<<endl;
	exit(EXIT_FAILURE);
     }
   dsytri("L", &dim1, d, &dim1, ipiv, work, &info);
   if(info != 0)
     {
	cerr<<"Error | kp_matrix::inverse_symmetric | some error in dsytri "<<info<<endl;
	exit(EXIT_FAILURE);
     }   
   delete [] ipiv;
   delete [] work;
*/
}
//                                                                                           
void kp_matrix::_create(int dim1_, int dim2_)
{
   dim1 = dim1_;
   dim2 = dim2_;
   d = new real[dim1 * dim2];
}
//                                                                                           
void kp_matrix::_clear()  
{
   if (d == NULL)
     {
	cerr<<"Error | kp_matrix::_clear| d == NULL"<<endl;
	exit(EXIT_FAILURE);
     }
   if (d != NULL)
     {
	delete [] d;	
     }
   d = NULL;
}
//                                                                                           
void kp_gemv(char op_A, real alpha, const kp_matrix& A, const kp_vector& x, real betta, kp_vector& y)
{
   int opA_dim1, opA_dim2;
   kp_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);
//   cout<<"gemv "<<op_A<<" "<<opA_dim1<<"x"<<opA_dim2<<endl;
   if (opA_dim2 != x.size() || opA_dim1 != y.size())
     {
	cerr<<"Error | kp_gemv | dimension problem"<<endl;
	exit(EXIT_FAILURE);
     }
   if (A.dim1 == 0 || x.d == 0 || y.d == 0)
     {
	cerr<<"A.dim1="<<A.dim1<<"x.d="<<x.d<<" y.d="<<y.d<<endl;
	cerr<<"Error | kp_gemv| leading dimension should be > 0"<<endl;
	exit(EXIT_FAILURE);
     }

   int one = 1;
#ifdef KP_SINGLE
   sgemv(&op_A, &A.dim1, &A.dim2, &alpha, A.d, &A.dim1, x.d, &one, &betta, y.d, &one);
#else
   dgemv(&op_A, &A.dim1, &A.dim2, &alpha, A.d, &A.dim1, x.d, &one, &betta, y.d, &one);
#endif
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_gemv((long)2 * (long)A.dim1 * (long)A.dim2);
   #endif   
}
//                                                                                           
void kp_gemm(char op_A, char op_B, real alpha, const kp_matrix& A, const kp_matrix& B, real betta, kp_matrix& C)
{
   int opA_dim1, opA_dim2;
   kp_check_op_set_dim(op_A, A, opA_dim1, opA_dim2);
   
   int opB_dim1, opB_dim2;
   kp_check_op_set_dim(op_B, B, opB_dim1, opB_dim2);
   
   if (C.dim1 != opA_dim1 || opA_dim2 != opB_dim1 || C.dim2 != opB_dim2)
     {
	cerr<<"Error | kp_gemm | diminsion problem"<<endl;
	exit(EXIT_FAILURE);
     }
   if (A.dim1 == 0 || B.dim1 == 0 || C.dim1 == 0)
     {
	cerr<<"A.dim1="<<A.dim1<<" B.dim1="<<B.dim1<<" C.dim1="<<C.dim1<<endl;
	cerr<<"Error | kp_gemm| leading dimension should be > 0"<<endl;
	exit(EXIT_FAILURE);
     }
#ifdef KP_SINGLE
   sgemm(&op_A, &op_B, &opA_dim1, &C.dim2,  &opA_dim2, &alpha, A.d, &A.dim1,
	 B.d, &B.dim1,  &betta, C.d, &C.dim1);
#else
   dgemm(&op_A, &op_B, &opA_dim1, &C.dim2,  &opA_dim2, &alpha, A.d, &A.dim1,
	 B.d, &B.dim1,  &betta, C.d, &C.dim1);
#endif
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_gemm((long)2 * (long)opA_dim1 * (long)C.dim2 * (long)opA_dim2);
   #endif   
}
//                                                                                             
void kp_syevd(kp_matrix& Q, kp_vector& DQ)
{
   if (Q.dim1 != Q.dim2 || Q.dim1 != DQ.size())
     {
	cerr<<"Error | kp_syevd | dimension problems"<<endl;
	exit(EXIT_FAILURE);
     }
   int n = Q.dim1;
   
   int lwork   = 2 * n * n +  6 * n + 1;
   real * work = new real[lwork];
   int liwork  = 5 * n + 3;
   int* iwork  = new int[liwork];
   int info;
   
#ifdef KP_SINGLE
   SSYEVD("V", "U", &n, Q.d, &n, DQ.d, work, &lwork, iwork, &liwork, &info);   
#else
   DSYEVD("V", "U", &n, Q.d, &n, DQ.d, work, &lwork, iwork, &liwork, &info);   
#endif
   delete [] work;
   delete [] iwork;
   
   #ifdef __WITH_FLOPS_CALCULATOR__   
   kp_fpc_global.add_eig((long)2 * (long)n * (long)n * (long)n);
   #endif   

}
//                                                                                            
void kp_check_op_set_dim(int op, const kp_matrix&M, int& dim1, int& dim2)
{
   if (op == 'N')
     {
	dim1 = M.dim1;
	dim2 = M.dim2;
     }
   else if (op == 'T')
     {
	dim1 = M.dim2;
	dim2 = M.dim1;
     }
   else
     {
	cerr<<"Error | kp_check_op_set_dim | op should ge either N either T"<<endl;
	exit(EXIT_FAILURE);
     }
}

//                                                                                            
double kp_calc_diff(const kp_matrix& m1, const kp_matrix& m2)
{
   if (m1.dim1 != m2.dim1 || m1.dim2 != m2.dim2)
     {
	cerr<<"error | kp_calc_diff | dimension problem"<<endl;
	exit(EXIT_FAILURE);
     }
   double rez = 0;
   for (int j = 0 ; j < m1.dim2; j++)
     for (int i = 0 ; i < m1.dim1 ; i++)
       rez  += fabs(m1(i,j) - m2(i,j));
   return rez;
}
//                                                                                            
void kp_vertcat(kp_matrix& rez, const vector<kp_matrix*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim1 += ms[i]->dim1;
	if (ms[i]->dim2 != ms[0]->dim2)
	  {
	     cerr<<"Error | kp_vertcat | inconsistant second dimension: "<<ms[i]->dim2<<" "<<ms[0]->dim2<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(sum_dim1, ms[0]->dim2);
   int dim1 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), dim1, 0);
	dim1 += ms[i]->dim1;
     }
}
//                                                                                  
void kp_horizcat(kp_matrix& rez, const vector<kp_matrix*>& ms)
{
   if (ms.size() == 0)
     {
	rez.resize(0,0);
	return;
     }
   int sum_dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	sum_dim2 += ms[i]->dim2;
	if (ms[i]->dim1 != ms[0]->dim1)
	  {
	     cerr<<"Error | kp_horizcat | inconsistant first dimension: "<<ms[i]->dim1<<" "<<ms[0]->dim1<<endl;
	     exit(EXIT_FAILURE);
	  }
     }
   rez.resize(ms[0]->dim1, sum_dim2);
   int dim2 = 0;
   for (size_t i = 0 ; i < ms.size() ; i++)
     {
	rez.set_from_submatrix(*(ms[i]), 0, dim2);
	dim2 += ms[i]->dim2;
     }
}
//                                                                                           
ostream& operator<<(ostream& out, kp_matrix& m)
{
   for (int i = 0 ; i < m.dim1 ; i++)
     {
	for (int j = 0 ; j < m.dim2 ; j++)
	  {
	     if (j)
	       out<<" ";
	     out<<m(i,j);
	  }
	out<<endl;
     }
   return out;
}
