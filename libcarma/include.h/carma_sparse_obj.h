/*
 * carma_sparse_obj.h
 *
 *  Created on: Apr 8, 2014
 *      Author: ???
 */

#ifndef CARMA_SPARSE_OBJ_H_
#define CARMA_SPARSE_OBJ_H_
#include <cusparse_v2.h>
#include "carma_obj.h"

template<class T_data>
class carma_sparse_host_obj;

template<class T_data>
class carma_sparse_obj {
public:
  long dims_data[3]; ///< dimensions of the array
  int nz_elem; ///< number of elements in the array
  int device; ///< device where the carma_obj is allocate
  carma_context *current_context;

  // ZERO-BASED INDEXING CSR-FORMAT
  T_data* d_data; //nz_elem elements
  int* d_rowind;  //dim1+1  elements
  int* d_colind;  //nz_elem elements
  cusparseMatDescr_t descr;

  char majorDim;

private:

public:
  carma_sparse_obj();
  carma_sparse_obj(carma_obj<T_data>* M);
  carma_sparse_obj(carma_sparse_obj<T_data>* M);
  carma_sparse_obj(carma_context *current_context, carma_sparse_host_obj<T_data>* M);
  carma_sparse_obj(carma_context *current_context, const long *dims, T_data * M, bool loadFromHost);
  virtual ~carma_sparse_obj();

  void operator=( carma_sparse_obj<T_data>& M);
  void operator=( carma_sparse_host_obj<T_data>& M);

  void resize(int nnz_, int dim1_, int dim2_);
  void init_from_transpose(carma_sparse_obj<T_data>* M);
  bool isColumnMajor();
  char get_majorDim() const {
    return majorDim;
  }
  void set_majorDim(char c) {
    majorDim = c;
  }

  /**< General Utilities */
  operator T_data*() {
    return d_data;
  }
  T_data* operator[](int index) {
    return &d_data[index];
  }
  T_data* getData() {
    return d_data;
  }
  T_data* getData(int index) {
    return &d_data[index];
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
  carma_context* getContext() {
    return current_context;
  }

  int getDevice() {
    return device;
  }

private:
  void _create(int nnz_, int dim1_, int dim2_);
  void _clear();

  template<
      cusparseStatus_t CUSPARSEAPI (*ptr_nnz)(cusparseHandle_t handle,
          cusparseDirection_t dirA, int m, int n, const cusparseMatDescr_t descrA,
          const T_data *A, int lda, int *nnzPerRowCol, int *nnzTotalDevHostPtr),
      cusparseStatus_t CUSPARSEAPI (*ptr_dense2csr)(cusparseHandle_t handle,
          int m, int n, const cusparseMatDescr_t descrA, const T_data *A,
          int lda, const int *nnzPerRow, T_data *csrValA, int *csrRowPtrA,
          int *csrColIndA)>
  void init_carma_sparse_obj(carma_context *current_context, const long *dims,
      T_data * M, bool loadFromHost);
};

template<class T_data>
cusparseStatus_t carma_gemv(cusparseHandle_t handle, char op_A, T_data alpha,
    carma_sparse_obj<T_data>* A, carma_obj<T_data>* x, T_data beta,
    carma_obj<T_data>* y);

template<class T_data>
cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, T_data alpha,
    carma_sparse_obj<T_data>* A, carma_obj<T_data>* B, T_data beta,
    carma_obj<T_data>* C);


#endif /* CARMA_SPARSE_OBJ_H_ */
