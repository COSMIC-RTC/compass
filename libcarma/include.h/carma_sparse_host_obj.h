/*
 * carmasparsehostobj.h
 *
 *  Created on: Apr 10, 2014
 *      Author: ???
 */

#ifndef CARMASPARSEHOSTOBJ_H_
#define CARMASPARSEHOSTOBJ_H_
#include "carma_host_obj.h"

template<class T_data>
class carma_sparse_host_obj {
public:
  carma_sparse_host_obj();
  carma_sparse_host_obj(carma_sparse_host_obj<T_data>& sm);
  carma_sparse_host_obj(const long *dims, T_data * M, char order);
  virtual ~carma_sparse_host_obj();

  //delete all arrays and create new for nnz=new_nnz
  void resize(int new_nnz, int dim1_, int dim2_);

  void operator=(carma_sparse_host_obj<T_data>& M);

  //we take from M only rows which in rowidx
  //analogue of MATLAB's this = M(idx,:)
  //idx is in normal zero based indexing
  //but we keep sparse matrix in one-based indexing
  void init_from_rowidx(carma_sparse_host_obj<T_data>* M,
      vector<int>* idx);

  //init from transpose sparce matrix
  void init_from_transpose(carma_sparse_host_obj<T_data>* M);
  void check();
  void init_from_matrix(const long *dims, T_data * M, char majorDim);
  void copy_into_matrix(T_data * M, char majorDim);
  void resize2rowMajor();
  void resize2colMajor();
  char get_majorDim() const {
    return majorDim;
  }

  /**< General Utilities */
  operator T_data*() {
    return h_data;
  }
  T_data* operator[](int index) {
    return &h_data[index];
  }
  T_data* getData() {
    return h_data;
  }
  T_data* getData(int index) {
    return &h_data[index];
  }
  const long *getDims() {
    return dims_data;
  }
  long getDims(int i) {
    return dims_data[i];
  }
  int getNzElem() {
    return nz_elem;
  }

private:
  void _create(int nnz_, int dim1_, int dim2_); //create new arrays
  void _clear(); //clear arrays

  char majorDim; //U - undefined
                 //R - row major
                 //C - col major

  MemAlloc mallocType; ///< type of host alloc

public:
  long dims_data[3]; ///< dimensions of the array
  int nz_elem; ///< number of elements in the array

  //see Sparse Matrix Storage Formats (coordinate format)
  //ONE-BASED!!!!
  T_data* h_data;
  int* rowind;
  int* colind;
};

//Multiply sparce matrix by vector
//y := alpha*A*x + betta * y
template<class T_data>
void carma_gemv(T_data alpha, carma_sparse_host_obj<T_data>* A,
    carma_host_obj<T_data>* x, T_data betta, carma_host_obj<T_data>* y,
    void (*ptr_coomv)(char *transa, long *m, long *k, T_data *alpha,
        char *matdescra, T_data *val, int *rowind, int *colind,
        int *nnz, T_data *x, T_data *beta, T_data *y));

//Multiply sparce matrix by dense matrix
//C := alpha*op(A)*B + betta * y
//op_A could be N (do nothing) or T (transpose)
template<class T_data>
void carma_gemm(char op_A, T_data alpha, carma_sparse_host_obj<T_data>* A,
    carma_host_obj<T_data>* B, T_data betta, carma_host_obj<T_data>* C);

//Multiply sparce matrix by dense matrix
//C := alpha*A*B + betta * y
template<class T_data>
inline void carma_gemm(T_data alpha, carma_sparse_host_obj<T_data>* A,
    carma_host_obj<T_data>* B, T_data betta, carma_host_obj<T_data>* C) {
  kp_gemm('N', alpha, A, B, betta, C);
}

#endif /* CARMASPARSEHOSTOBJ_H_ */
