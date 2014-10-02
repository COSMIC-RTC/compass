//kp_smatrix.h
//kalman project sparse matrix
//we save matrix in coordinate format from MKL (see Sparse Matrix Storage Formats in MKL)
//we use one-based indexing (because one-based indexing would need row first packing for dens arrays in kp_gemm)

#ifndef __SEGER__KP_SMATRIX_H__
#define __SEGER__KP_SMATRIX_H__

#include "kp_vector.h"
#include "kp_matrix.h"


class kp_smatrix
{
 public:
   kp_smatrix();
   kp_smatrix(const kp_smatrix& sm);
   ~kp_smatrix();
   
   //delete all arrays and create new for nnz=new_nnz
   void resize(int new_nnz, int dim1_, int dim2_); 
   
   void operator=(const kp_smatrix& M);
   
   //we take from M only rows which in rowidx
   //analogue of MATLAB's this = M(idx,:)
   //idx is in normal zero based indexing
   //but we keep sparse matrix in one-based indexing
   void init_from_rowidx(const kp_smatrix& M, const vector<int>& idx);
   
   //init from transpose sparce matrix
   void init_from_transpose(const kp_smatrix& M);
   void check();
   void init_from_matrix(const kp_matrix& B, double epsilon = 0.0);
   void resize2rowMajor();
   void resize2colMajor();
   char get_majorDim()const  {return majorDim;}
   
 private:
   void _create(int nnz_, int dim1_, int dim2_); //create new arrays
   void _clear(); //clear arrays
   
   char majorDim; //U - undefined
                  //R - row major
                  //C - col major
   
 public:
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
//y := alpha*A*x + betta * y
void kp_gemv(real alpha, kp_smatrix& A, kp_vector& x, real betta, kp_vector& y);

//Multiply sparce matrix by dense matrix
//C := alpha*op(A)*B + betta * y
//op_A could be N (do nothing) or T (transpose)
void kp_gemm(char op_A, real alpha, kp_smatrix& A, kp_matrix& B, real betta, kp_matrix& C);

//Multiply sparce matrix by dense matrix
//C := alpha*A*B + betta * y
inline void kp_gemm(real alpha, kp_smatrix& A, kp_matrix& B, real betta, kp_matrix& C)
{kp_gemm('N', alpha, A, B, betta, C);}

void kp_check_op_set_dim(int op, const kp_smatrix&M, int& dim1, int& dim2);

#endif
