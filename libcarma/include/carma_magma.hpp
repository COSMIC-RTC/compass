// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_magma.hpp
//! \ingroup   libcarma
//! \brief     this file provides wrappers to the magma functions
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _CARMA_MAGMA_H_
#define _CARMA_MAGMA_H_

#include <carma_host_obj.hpp>
#include <carma_obj.hpp>
#include <carma_sparse_obj.hpp>

// MAGMA functions
int32_t carma_magma_disabled();
// template <class T>
// int32_t carma_svd(CarmaObj<T> *imat, CarmaObj<T> *eigenvals,
//               CarmaObj<T> *mod2act, CarmaObj<T> *mes2mod);
template <class T>
int32_t carma_magma_syevd(char jobz, CarmaObj<T> *mat,
                      CarmaHostObj<T> *eigenvals);
// template <class T, int32_t method>
// int32_t carma_magma_syevd(char jobz, CarmaObj<T> *mat, CarmaHostObj<T>
// *eigenvals);
template <class T>
int32_t carma_magma_syevd_m(int64_t ngpu, char jobz, int64_t N, T *mat, T *eigenvals);
template <class T>
int32_t carma_magma_syevd_m(int64_t ngpu, char jobz, CarmaHostObj<T> *mat,
                        CarmaHostObj<T> *eigenvals);
template <class T>
int32_t carma_magma_syevd_m(int64_t ngpu, char jobz, CarmaHostObj<T> *mat,
                        CarmaHostObj<T> *eigenvals, CarmaHostObj<T> *U);
template <class T>
int32_t carma_magma_getri(CarmaObj<T> *d_iA);
template <class T>
int32_t carma_magma_potr_inv(CarmaObj<T> *d_iA);
template <class T>
int32_t carma_magma_potr_inv_m(int64_t num_gpus, CarmaHostObj<T> *h_A,
                        CarmaObj<T> *d_iA);

// MAGMA functions (direct access)
template <class T>
int32_t carma_magma_syevd(char jobz, int64_t N, T *mat, T *eigenvals);
// template <class T>
// int32_t carma_magma_potr_inv_m(int64_t num_gpus, int64_t N, T *h_A, T *d_iA);

template <class T_data>
int32_t carma_magma_svd_cpu(CarmaHostObj<T_data> *imat,
                        CarmaHostObj<T_data> *eigenvals,
                        CarmaHostObj<T_data> *mod2act,
                        CarmaHostObj<T_data> *mes2mod);
template <class T>
int32_t carma_magma_getri_cpu(CarmaHostObj<T> *h_A);
template <class T>
int32_t carma_magma_potr_inv_cpu(CarmaHostObj<T> *h_A);
template <class T>
int32_t carma_magma_syevd_cpu(char jobz, CarmaHostObj<T> *h_A,
                          CarmaHostObj<T> *eigenvals);

// MAGMA functions (direct access)
// template <class T>
// int32_t carma_magma_svd_cpu(int64_t N, int64_t M, T *imat, T *eigenvals, T *mod2act,
//                   T *mes2mod);
template <class T>
int32_t carma_magma_getri_cpu(int64_t N, T *h_A);
template <class T>
int32_t carma_magma_potr_inv_cpu(int64_t N, T *h_A);
template <class T>
int32_t carma_magma_syevd_cpu(char jobz, int64_t N, T *h_A, T *eigenvals);
template <class T>
int32_t carma_magma_axpy_cpu(int64_t N, T alpha, T *h_X, int64_t incX, T *h_Y, int64_t incY);
template <class T>
int32_t carma_gemm_cpu(char transa, char transb, int64_t m, int64_t n, int64_t k, T alpha,
                   T *A, int64_t lda, T *B, int64_t ldb, T beta, T *C, int64_t ldc);

template <class T_data>
int32_t carma_magma_gemv(char trans, int32_t m, int32_t n, T_data alpha, T_data *matA,
                     int32_t lda, T_data *vectx, int32_t incx, T_data beta,
                     T_data *vecty, int32_t incy);

template <class T_data>
int32_t carma_magma_csr2ell(CarmaSparseObj<T_data> *dA);

template <class T_data>
int32_t carma_magma_spmv(T_data alpha, CarmaSparseObj<T_data> *dA,
                     CarmaObj<T_data> *dx, T_data beta, CarmaObj<T_data> *dy);

template <class T_data>
int32_t carma_magma_sparse_free(CarmaSparseObj<T_data> *dA);

#endif  // _CARMA_MAGMA_H_
