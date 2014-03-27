//kp_matrix.h
//Kalman project matrix

#ifndef __SEGER__KP_MATRIX_H__
#define __SEGER__KP_MATRIX_H__

#include "kp_vector.h"
#include <iostream>
class kp_smatrix;  //pre declaration of kp_smatrix


class kp_matrix
{
 public:
   kp_matrix();
   kp_matrix(int dim1_, int dim2_);
   kp_matrix(const kp_matrix& m);
   kp_matrix(const vector< vector<double> >& m);
   ~kp_matrix();   
   void resize(int dim1_, int dim2_);
   
   real&        operator()(int i , int j)       { return d[j * dim1 + i] ;}
   const real&  operator()(int i , int j) const { return d[j * dim1 + i] ;}
   real&          el(int i , int j) { return d[j * dim1 + i] ;}
   
   real& at(int i, int j);
   
   void operator=( const kp_matrix& m);
   void zeros();  //put all elements to 0
   
   //apply operation to each element
   void operator /= (real val);
   void operator *= (real val);
   void operator += (real val);
   void operator -= (real val);
   
   
   void operator += (const kp_matrix& M);

   
   //analogue of mean(M,1) or mean(M) from MATLAB 
   //return vector of size dim2 contain means of each column
   void cols_mean(kp_vector& rez);

   //analogue of mean(M,2) from MATLAB 
   //return vector of size dim1 contain means of each row
   void rows_mean(kp_vector& rez);
   
   //calculate trace
   double trace();
   
   void init_from_matrix(const kp_matrix& M, int r1, int r2, int c1, int c2);
   //we take from M only rows which in rowidx
   //analogue of MATLAB's this = M(rowidx,:)
   void init_from_rowidx(const kp_matrix& M, const vector<int>& rowidx);

   //subM is matrix which was created by init_from_rowidx 
   //analogue of MATLAB's this(rowidx, : ) = subM
   void set_from_rowsubmatrix(const kp_matrix& subM, const vector<int>& rowidx);
   
   //equivalent of MATLAB
   //find(this(:,col));
   void make_rowidx(int col, vector<int>& rowidx);
   
   //init from transpose matrix
   void init_from_transpose(const kp_matrix& M); 
   
   //MATLAB code:
   //this = diag(v) * this
   //we multiplay ith row by v[i]
   void mult_each_row(const kp_vector&v);
   
   //put submatrix at  r,c -> r + subM.dim1, c + subM.dim2 
   void set_from_submatrix(const kp_matrix& subM, int r, int c);
   
   void init_from_smatrix(const kp_smatrix& m);
   
   // MATLAB CODE: M = (M - repmat(v,1,M.dim2)) * val
   // So we substract v from each column of M and after multiply result by val
   void each_column_substract_mult(const kp_vector &v, double val);
   
   //M(i,j) += v[i]
   void each_column_add(const kp_vector &v);   
   void add_unit_to_diagonal();

   //run dgetrf/dgetri
   void inverse();
   //run dsytrf/dsytri (inverse symmetric matrix);
   void inverse_symmetric();
      
 private:
   void _create(int dim1_, int dim2_);
   void _clear(); 
 public:
   //dim1 -- number of rows    (usually m)
   //dim2 -- number of columns (usually n)
   int dim1, dim2; 
   //tda == dim1 !!!
   real* d;
};



//y := alpha*op(A)*x + betta * y
//op_A could be N (do nothing) or T (transpose)
void kp_gemv(char op_A, real alpha, const kp_matrix& A, const kp_vector& x, real betta, kp_vector& y);

//y := alpha*A*x + betta * y
inline void kp_gemv(real alpha, const kp_matrix& A, const kp_vector& x, real betta, kp_vector& y)
{ kp_gemv('N', alpha, A, x, betta, y); }



//C := alpha*op(A)*op(B) + betta * y
//op_? could be N (do nothing) or T (transpose)
void kp_gemm(char op_A, char op_B, real alpha, const kp_matrix& A, 
	     const kp_matrix& B, real betta, kp_matrix& C);

//C := alpha*A*B + betta * y
inline void kp_gemm(real alpha, const kp_matrix& A, const kp_matrix& B, real betta, kp_matrix& C)
{ kp_gemm('N', 'N', alpha, A, B, betta, C);}


//analogue of mean in MATLAB 
//if dim == 1 -> calculate means of each columns 
//if dim == 2 -> calculate means of each rows
void kp_mean(const kp_matrix& M, int dim, kp_vector& rez);

//eigen-value decomposition of symmetric matrix
//Q must be symmetrical matrix (we will not check it)
//Q will be replaced my matrix of eigen vectors
//DQ - vector of eigen values value
//Q_in = OQ_out * DQ * OQ_out'
void kp_syevd(kp_matrix&Q, kp_vector& DQ);


void kp_check_op_set_dim(int op, const kp_matrix&M, int& dim1, int& dim2);

double kp_calc_diff(const kp_matrix& m1, const kp_matrix& m2);

//make vertcat of matrises
void kp_vertcat(kp_matrix& rez, const vector<kp_matrix*>& ms);
//make horizcat of matrises
void kp_horizcat(kp_matrix& rez, const vector<kp_matrix*>& ms);

ostream& operator<<(ostream& out, kp_matrix& m);


#endif
