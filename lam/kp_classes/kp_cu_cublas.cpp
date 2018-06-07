#include "kp_cu_cublas.h"

template <>
void kp_cu_geam_core<double>(cublasHandle_t &handle, cublasOperation_t &transa,
                             cublasOperation_t &transb, int &m, int &n,
                             const double *alpha, const double *A_data,
                             int &lda, const double *beta, const double *B_data,
                             int &ldb, double *C_data, int &ldc) {
  cublasStatus_t stat;
  stat = cublasDgeam(handle, transa, transb, m, n, alpha, A_data, lda, beta,
                     B_data, ldb, C_data, ldc);

  if (stat != CUBLAS_STATUS_SUCCESS) {
    cerr << "error | kp_cu_matrix::init_from_transpose | erreur lors de la "
            "transposition sur GPU"
         << endl;
    throw KP_CUBLAS_GEAM;
    // exit(EXIT_FAILURE);
  }
}
template <>
void kp_cu_geam_core<float>(cublasHandle_t &handle, cublasOperation_t &transa,
                            cublasOperation_t &transb, int &m, int &n,
                            const float *alpha, const float *A_data, int &lda,
                            const float *beta, const float *B_data, int &ldb,
                            float *C_data, int &ldc) {
  cublasStatus_t stat;
  stat = cublasSgeam(handle, transa, transb, m, n, alpha, A_data, lda, beta,
                     B_data, ldb, C_data, ldc);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    cerr << "error | kp_cu_matrix::init_from_transpose | erreur lors de la "
            "transposition sur GPU"
         << endl;
    throw KP_CUBLAS_GEAM;
    // exit(EXIT_FAILURE);
  }
}

template <>
void kp_cu_gemm_core<double>(cublasHandle_t &handle, cublasOperation_t &transa,
                             cublasOperation_t &transb, int m, int n, int k,
                             const double *alpha, const double *A_data, int lda,
                             const double *B_data, int ldb, const double *beta,
                             double *C_data, int ldc) {
  cublasStatus_t stat1;
  stat1 = cublasDgemm(handle, transa, transb, m, n, k, alpha, A_data, lda,
                      B_data, ldb, beta, C_data, ldc);

  if (stat1 != CUBLAS_STATUS_SUCCESS) {
    cerr << "error | kp_cu_matrix | kp_cu_gemm failed : " << stat1 << endl;
    throw KP_CUBLAS_GEMM;
  }
}

template <>
void kp_cu_gemm_core<float>(cublasHandle_t &handle, cublasOperation_t &transa,
                            cublasOperation_t &transb, int m, int n, int k,
                            const float *alpha, const float *A_data, int lda,
                            const float *B_data, int ldb, const float *beta,
                            float *C_data, int ldc) {
  cublasStatus_t stat1;
  stat1 = cublasSgemm(handle, transa, transb, m, n, k, alpha, A_data, lda,
                      B_data, ldb, beta, C_data, ldc);

  if (stat1 != CUBLAS_STATUS_SUCCESS) {
    cerr << "error | kp_cu_matrix | kp_cu_gemm failed : " << stat1 << endl;
    throw KP_CUBLAS_GEMM;
  }
}

template <>
void kp_cu_gemv_core<double>(cublasHandle_t &handle, cublasOperation_t &trans,
                             int m, int n, const double *alpha,
                             const double *A_data, int lda,
                             const double *x_data, int incx, const double *beta,
                             double *y_data, int incy) {
  cublasStatus_t cublasStat;

  cublasStat = cublasDgemv(handle, trans, m, n, alpha, A_data, lda, x_data,
                           incx, beta, y_data, incy);

  if (cublasStat != CUBLAS_STATUS_SUCCESS) {
    cerr << "error | kp_cu_matrix | kp_cu_gemv failed !" << endl;
    throw KP_CUBLAS_GEMV;
  }
}
template <>
void kp_cu_gemv_core<float>(cublasHandle_t &handle, cublasOperation_t &trans,
                            int m, int n, const float *alpha,
                            const float *A_data, int lda, const float *x_data,
                            int incx, const float *beta, float *y_data,
                            int incy) {
  cublasStatus_t cublasStat;

  cublasStat = cublasSgemv(handle, trans, m, n, alpha, A_data, lda, x_data,
                           incx, beta, y_data, incy);

  if (cublasStat != CUBLAS_STATUS_SUCCESS) {
    cerr << "error | kp_cu_matrix | kp_cu_gemv failed !" << endl;
    throw KP_CUBLAS_GEMV;
  }
}

template <>
void kp_cu_gemv_core<double>(cusparseHandle_t &handle,
                             cusparseOperation_t transA, int m, int n, int nnz,
                             const double *alpha,
                             const cusparseMatDescr_t &descrA,
                             const double *A_values, const int *A_csrRowPtr,
                             const int *A_colind_cu, const double *x_data,
                             const double *beta, double *y_data) {
  cusparseStatus_t status;
  status = cusparseDcsrmv(handle, transA, m, n, nnz, alpha, descrA, A_values,
                          A_csrRowPtr, A_colind_cu, x_data, beta, y_data);

  if (status != CUSPARSE_STATUS_SUCCESS) {
    cerr << "Error | kp_cu_gemv (sparse) | Matrix-vector multiplication failed"
         << endl;
    throw KP_CUSPARSE_GEMV;
    // exit(EXIT_FAILURE);
  }
}
template <>
void kp_cu_gemv_core<float>(cusparseHandle_t &handle,
                            cusparseOperation_t transA, int m, int n, int nnz,
                            const float *alpha,
                            const cusparseMatDescr_t &descrA,
                            const float *A_values, const int *A_csrRowPtr,
                            const int *A_colind_cu, const float *x_data,
                            const float *beta, float *y_data) {
  cusparseStatus_t status;
  status = cusparseScsrmv(handle, transA, m, n, nnz, alpha, descrA, A_values,
                          A_csrRowPtr, A_colind_cu, x_data, beta, y_data);

  if (status != CUSPARSE_STATUS_SUCCESS) {
    cerr << "Error | kp_cu_gemv (sparse) | Matrix-vector multiplication failed"
         << endl;
    throw KP_CUSPARSE_GEMV;
    // exit(EXIT_FAILURE);
  }
}

template <>
void kp_cu_gemm_core<double>(cusparseHandle_t &handle,
                             cusparseOperation_t transA, int m, int n, int k,
                             int A_nnz, const double *alpha,
                             const cusparseMatDescr_t &descrA,
                             const double *A_values_cu, const int *A_csrRowPtr,
                             const int *A_colind_cu, const double *B_data,
                             int ldb, const double *beta, double *C_data,
                             int ldc) {
  cusparseStatus_t status;
  status =
      cusparseDcsrmm(handle, transA, m, n, k, A_nnz, alpha, descrA, A_values_cu,
                     A_csrRowPtr, A_colind_cu, B_data, ldb, beta, C_data, ldc);

  if (status != CUSPARSE_STATUS_SUCCESS) {
    cerr << "Error | kp_cu_gemm (sparse) | Matrix-matrix multiplication failed"
         << endl;
    throw KP_CUSPARSE_GEMM;
    // exit(EXIT_FAILURE);
  }
}
template <>
void kp_cu_gemm_core<float>(
    cusparseHandle_t &handle, cusparseOperation_t transA, int m, int n, int k,
    int A_nnz, const float *alpha, const cusparseMatDescr_t &descrA,
    const float *A_values_cu, const int *A_csrRowPtr, const int *A_colind_cu,
    const float *B_data, int ldb, const float *beta, float *C_data, int ldc) {
  cusparseStatus_t status;
  status =
      cusparseScsrmm(handle, transA, m, n, k, A_nnz, alpha, descrA, A_values_cu,
                     A_csrRowPtr, A_colind_cu, B_data, ldb, beta, C_data, ldc);
  if (status != CUSPARSE_STATUS_SUCCESS) {
    cerr << "Error | kp_cu_gemm (sparse) | Matrix-matrix multiplication failed"
         << endl;
    throw KP_CUSPARSE_GEMM;
    // exit(EXIT_FAILURE);
  }
}

template <>
void kp_cu_sgemm_core<double>(
    cusparseHandle_t &handle, cusparseOperation_t transA,
    cusparseOperation_t transB, int m, int n, int k,
    const cusparseMatDescr_t &descrA, const int A_nnz,
    const double *A_values_cu, const int *A_csrRowPtr, const int *A_colind_cu,
    const cusparseMatDescr_t &descrB, const int B_nnz,
    const double *B_values_cu, const int *B_csrRowPtr, const int *B_colind_cu,
    const cusparseMatDescr_t &descrC, double *C_values_cu,
    const int *C_csrRowPtr, int *C_colind_cu) {
  cusparseStatus_t status;
  status = cusparseDcsrgemm(handle, transA, transB, m, n, k, descrA, A_nnz,
                            A_values_cu, A_csrRowPtr, A_colind_cu, descrB,
                            B_nnz, B_values_cu, B_csrRowPtr, B_colind_cu,
                            descrC, C_values_cu, C_csrRowPtr, C_colind_cu);

  if (status != CUSPARSE_STATUS_SUCCESS) {
    cerr << "Error | kp_cu_sgemm (sparse) | Conversion from COO to CSR format "
            "failed"
         << endl;
    // throw KP_CUSPARSE_COO2CSR;
    // exit(EXIT_FAILURE);
  }
}
template <>
void kp_cu_sgemm_core<float>(
    cusparseHandle_t &handle, cusparseOperation_t transA,
    cusparseOperation_t transB, int m, int n, int k,
    const cusparseMatDescr_t &descrA, const int A_nnz, const float *A_values_cu,
    const int *A_csrRowPtr, const int *A_colind_cu,
    const cusparseMatDescr_t &descrB, const int B_nnz, const float *B_values_cu,
    const int *B_csrRowPtr, const int *B_colind_cu,
    const cusparseMatDescr_t &descrC, float *C_values_cu,
    const int *C_csrRowPtr, int *C_colind_cu) {
  cusparseStatus_t status;
  status = cusparseScsrgemm(handle, transA, transB, m, n, k, descrA, A_nnz,
                            A_values_cu, A_csrRowPtr, A_colind_cu, descrB,
                            B_nnz, B_values_cu, B_csrRowPtr, B_colind_cu,
                            descrC, C_values_cu, C_csrRowPtr, C_colind_cu);

  if (status != CUSPARSE_STATUS_SUCCESS) {
    cerr << "Error | kp_cu_sgemm (sparse) | Conversion from COO to CSR format "
            "failed"
         << endl;
    // throw KP_CUSPARSE_COO2CSR;
    // exit(EXIT_FAILURE);
  }
}

void kp_cu_op2cublastrans(int op, cublasOperation_t *trans) {
  if (op == 'N')
    *trans = CUBLAS_OP_N;
  else if (op == 'T')
    *trans = CUBLAS_OP_T;
  else {
    cerr << "Error | kp_cu_op2cublastrans | op should be either N or T" << endl;
    exit(EXIT_FAILURE);
  }
}
