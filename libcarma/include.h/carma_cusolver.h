// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_cusolver.h
//! \ingroup   libcarma
//! \brief     this file provides wrappers to the cuSolver functions
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#ifndef _CARMA_CUSOLVER_H_
#define _CARMA_CUSOLVER_H_

#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>
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
int carma_syevd(char jobz, CarmaObj<T> *mat,
                CarmaObj<T> *eigenvals);
// template <class T, int method>
// int carma_syevd(cusolverDnHandle_t cusolver_handle, char jobz, CarmaObj<T> *mat, CarmaHostObj<T> *eigenvals);
// template <class T>
// int carma_syevd_m(cusolverDnHandle_t cusolver_handle, long ngpu, char jobz, long N, T *mat, T *eigenvals);
// template <class T>
// int carma_syevd_m(cusolverDnHandle_t cusolver_handle, long ngpu, char jobz, CarmaHostObj<T> *mat,
//                   CarmaHostObj<T> *eigenvals);
// template <class T>
// int carma_syevd_m(cusolverDnHandle_t cusolver_handle, long ngpu, char jobz, CarmaHostObj<T> *mat,
//                   CarmaHostObj<T> *eigenvals, CarmaHostObj<T> *U);
// template <class T>
// int carma_getri(CarmaObj<T> *d_iA);
template <class T>
int carma_potr_inv(CarmaObj<T> *d_iA);
// template <class T>
// int carma_potr_inv_m(long num_gpus, CarmaHostObj<T> *h_A, CarmaObj<T> *d_iA);

// MAGMA functions (direct access)
template <class T>
int carma_syevd(cusolverDnHandle_t cusolver_handle, char jobz, long N, T *mat, T *eigenvals);
template <class T>
int carma_potr_inv(cusolverDnHandle_t cusolver_handle, long N, T *h_A);
// template <class T>
// int carma_potr_inv_m(cusolverDnHandle_t cusolver_handle, long num_gpus, long N, T *h_A, T *d_iA);

// template <class T_data>
// int carma_svd_cpu(CarmaHostObj<T_data> *imat,
//                   CarmaHostObj<T_data> *eigenvals,
//                   CarmaHostObj<T_data> *mod2act,
//                   CarmaHostObj<T_data> *mes2mod);
// template <class T>
// int carma_getri_cpu(CarmaHostObj<T> *h_A);
// template <class T>
// int carma_potr_inv_cpu(CarmaHostObj<T> *h_A);
// template <class T>
// int carma_syevd_cpu(char jobz, CarmaHostObj<T> *h_A,
//                     CarmaHostObj<T> *eigenvals);

// // MAGMA functions (direct access)
// // template <class T>
// // int carma_svd_cpu(long N, long M, T *imat, T *eigenvals, T *mod2act,
// //                   T *mes2mod);
// template <class T>
// int carma_getri_cpu(long N, T *h_A);
// template <class T>
// int carma_potr_inv_cpu(long N, T *h_A);
// template <class T>
// int carma_syevd_cpu(char jobz, long N, T *h_A, T *eigenvals);
// template <class T>
// int carma_axpy_cpu(long N, T alpha, T *h_X, long incX, T *h_Y, long incY);
// template <class T>
// int carma_gemm_cpu(char transa, char transb, long m, long n, long k, T alpha,
//                    T *A, long lda, T *B, long ldb, T beta, T *C, long ldc);

// template <class T_data>
// int carma_cusolver_csr2ell(CarmaSparseObj<T_data> *dA);

// template <class T_data>
// int carma_cusolver_spmv(T_data alpha, CarmaSparseObj<T_data> *dA,
//                      CarmaObj<T_data> *dx, T_data beta, CarmaObj<T_data>
//                      *dy);

// template <class T_data>
// int carma_sparse_free(CarmaSparseObj<T_data> *dA);

#endif  // _CARMA_CUSOLVER_H_
