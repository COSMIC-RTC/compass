// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_cusolver.hpp
//! \ingroup   libcarma
//! \brief     this file provides wrappers to the cuSolver functions
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _CARMA_CUSOLVER_H_
#define _CARMA_CUSOLVER_H_

#include <carma_host_obj.hpp>
#include <carma_obj.hpp>
#include <carma_sparse_obj.hpp>
#include <cusolverDn.h>

#ifndef SOLVER_EIG_MODE_VECTOR
  #define SOLVER_EIG_MODE_VECTOR 'V'
#endif
#ifndef SOLVER_EIG_MODE_NOVECTOR
  #define SOLVER_EIG_MODE_NOVECTOR 'N'
#endif

cusolverStatus_t carma_init_cusolver(cusolverDnHandle_t *cusolver_handle);
cusolverStatus_t carma_shutdown_cusolver(cusolverDnHandle_t cusolver_handle);

template <class T>
int32_t carma_syevd(char jobz, CarmaObj<T> *mat,
                CarmaObj<T> *eigenvals);
// template <class T, int32_t method>
// int32_t carma_syevd(cusolverDnHandle_t cusolver_handle, char jobz, CarmaObj<T> *mat, CarmaHostObj<T> *eigenvals);
// template <class T>
// int32_t carma_syevd_m(cusolverDnHandle_t cusolver_handle, int64_t ngpu, char jobz, int64_t N, T *mat, T *eigenvals);
// template <class T>
// int32_t carma_syevd_m(cusolverDnHandle_t cusolver_handle, int64_t ngpu, char jobz, CarmaHostObj<T> *mat,
//                   CarmaHostObj<T> *eigenvals);
// template <class T>
// int32_t carma_syevd_m(cusolverDnHandle_t cusolver_handle, int64_t ngpu, char jobz, CarmaHostObj<T> *mat,
//                   CarmaHostObj<T> *eigenvals, CarmaHostObj<T> *U);
// template <class T>
// int32_t carma_getri(CarmaObj<T> *d_iA);
template <class T>
int32_t carma_potr_inv(CarmaObj<T> *d_iA);
// template <class T>
// int32_t carma_potr_inv_m(int64_t num_gpus, CarmaHostObj<T> *h_A, CarmaObj<T> *d_iA);

template <class T>
int32_t carma_syevd(cusolverDnHandle_t cusolver_handle, char jobz, int64_t N, T *mat, T *eigenvals);
template <class T>
int32_t carma_potr_inv(cusolverDnHandle_t cusolver_handle, int64_t N, T *h_A);
// template <class T>
// int32_t carma_potr_inv_m(cusolverDnHandle_t cusolver_handle, int64_t num_gpus, int64_t N, T *h_A, T *d_iA);

// template <class T_data>
// int32_t carma_svd_cpu(CarmaHostObj<T_data> *imat,
//                   CarmaHostObj<T_data> *eigenvals,
//                   CarmaHostObj<T_data> *mod2act,
//                   CarmaHostObj<T_data> *mes2mod);
// template <class T>
// int32_t carma_getri_cpu(CarmaHostObj<T> *h_A);
// template <class T>
// int32_t carma_potr_inv_cpu(CarmaHostObj<T> *h_A);
// template <class T>
// int32_t carma_syevd_cpu(char jobz, CarmaHostObj<T> *h_A,
//                     CarmaHostObj<T> *eigenvals);

// template <class T>
// int32_t carma_getri_cpu(int64_t N, T *h_A);
// template <class T>
// int32_t carma_potr_inv_cpu(int64_t N, T *h_A);
// template <class T>
// int32_t carma_syevd_cpu(char jobz, int64_t N, T *h_A, T *eigenvals);
// template <class T>
// int32_t carma_axpy_cpu(int64_t N, T alpha, T *h_X, int64_t incX, T *h_Y, int64_t incY);
// template <class T>
// int32_t carma_gemm_cpu(char transa, char transb, int64_t m, int64_t n, int64_t k, T alpha,
//                    T *A, int64_t lda, T *B, int64_t ldb, T beta, T *C, int64_t ldc);

// template <class T_data>
// int32_t carma_cusolver_csr2ell(CarmaSparseObj<T_data> *dA);

// template <class T_data>
// int32_t carma_cusolver_spmv(T_data alpha, CarmaSparseObj<T_data> *dA,
//                      CarmaObj<T_data> *dx, T_data beta, CarmaObj<T_data>
//                      *dy);

// template <class T_data>
// int32_t carma_sparse_free(CarmaSparseObj<T_data> *dA);

#endif  // _CARMA_CUSOLVER_H_
