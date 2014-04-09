/*
 * carma_sparse_obj.h
 *
 *  Created on: Apr 8, 2014
 *      Author: ???
 */

#ifndef CARMA_SPARSE_OBJ_H_
#define CARMA_SPARSE_OBJ_H_

template<class T_data>
class carma_sparse_obj {
public:
  int* csrRowPtr;
  int* csrRowPtrT;

  // ONE-BASED INDEXING
  T_data* values_cu;
  int* rowind;
  int* colind;
  cusparseMatDescr_t descr;

  int nnz;
  int dim1; //nombre de lignes
  int dim2; //nombre de colones

private:

  char majorDim;
  bool isCSRconverted;
  bool isCSRconvertedT;

public:
  carma_sparse_obj();
  carma_sparse_obj(const carma_sparse_obj<T_data>* M);
  carma_sparse_obj(const carma_obj<T_data>* M);
  virtual ~carma_sparse_obj();

  void operator=(const carma_sparse_obj<T_data> M);
  void operator=(const carma_obj<T_data> M);

  void resize(int nnz_, int dim1_, int dim2_);
  void init_from_transpose(const carma_sparse_obj<T_data>* M);
  bool isColumnMajor();
  char get_majorDim() const {
    return majorDim;
  }
  void set_majorDim(char c) {
    majorDim = c;
  }

  friend void carma_gemm(cusparseHandle_t handle, char op_A, real alpha,
      carma_sparse_obj<T_data>* A, carma_obj<T_data>* B, real beta, carma_obj<T_data>* C);
  //friend void kp_cu_gemm(cusparseHandle_t handle,cublasHandle_t cublasHandle, char op_A, char op_B, real alpha, kp_cu_smatrix& A, kp_cu_matrix& B, real beta, kp_cu_matrix& C);

  friend void carma_gemv(cusparseHandle_t handle, char op_A, real alpha,
      carma_sparse_obj<T_data>* A, const carma_obj<T_data>* x, real beta,
      carma_obj<T_data>* y);

private:
  void _create(int nnz_, int dim1_, int dim2_);
  void _clear();

};

// C = alpha * op_A(A) * B + beta * C
void carma_gemm(cusparseHandle_t handle, char op_A, real alpha,
    carma_sparse_obj<T_data>* A, const carma_obj<T_data>* B, real beta, carma_obj<T_data>* C);

// C = alpha * A * B + beta * C
inline void carma_gemm(cusparseHandle_t handle, cublasHandle_t cublashandle,
    real alpha, carma_sparse_obj<T_data>* A, const carma_obj<T_data>* B, real beta,
    carma_obj<T_data>* C) {
  carma_gemm(handle, 'N', alpha, A, B, beta, C);
}

// y = alpha * op_A(A) * x + beta * y
void carma_gemv(cusparseHandle_t handle, char op_A, real alpha,
    carma_sparse_obj<T_data>* A, const carma_obj<T_data>* x, real beta, carma_obj<T_data>* y);

#endif /* CARMA_SPARSE_OBJ_H_ */
