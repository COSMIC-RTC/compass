#include "kp_smatrix.h"

template <>
void kp_gemm_core<double>(char *transa, MKL_INT *m, MKL_INT *n, MKL_INT *k,
                          double *alpha, char *matdescra, double *A_values,
                          MKL_INT *A_rowind, MKL_INT *A_colind, MKL_INT *A_nnz,
                          double *B_data, MKL_INT *ldb, double *beta,
                          double *C_data, MKL_INT *ldc) {
  mkl_dcoomm(transa, m, n, k, alpha, matdescra, A_values, A_rowind, A_colind,
             A_nnz, B_data, ldb, beta, C_data, ldc);
}
template <>
void kp_gemm_core<float>(char *transa, MKL_INT *m, MKL_INT *n, MKL_INT *k,
                         float *alpha, char *matdescra, float *A_values,
                         MKL_INT *A_rowind, MKL_INT *A_colind, MKL_INT *A_nnz,
                         float *B_data, MKL_INT *ldb, float *beta,
                         float *C_data, MKL_INT *ldc) {
  mkl_scoomm(transa, m, n, k, alpha, matdescra, A_values, A_rowind, A_colind,
             A_nnz, B_data, ldb, beta, C_data, ldc);
}

template <>
void kp_gemv_core<double>(char *transa, MKL_INT *m, MKL_INT *k, double *alpha,
                          char *matdescra, double *A_values, MKL_INT *A_rowind,
                          MKL_INT *A_colind, MKL_INT *A_nnz, double *x_data,
                          double *beta, double *y_data) {
  mkl_dcoomv(transa, m, k, alpha, matdescra, A_values, A_rowind, A_colind,
             A_nnz, x_data, beta, y_data);
}
template <>
void kp_gemv_core<float>(char *transa, MKL_INT *m, MKL_INT *k, float *alpha,
                         char *matdescra, float *A_values, MKL_INT *A_rowind,
                         MKL_INT *A_colind, MKL_INT *A_nnz, float *x_data,
                         float *beta, float *y_data) {
  mkl_scoomv(transa, m, k, alpha, matdescra, A_values, A_rowind, A_colind,
             A_nnz, x_data, beta, y_data);
}
