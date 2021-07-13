// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the
//  terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for
//  the simulation of AO systems.
//
//  The final product includes a software package for simulating all the
//  critical subcomponents of AO, particularly in the context of the ELT and a
//  real-time core based on several control approaches, with performances
//  consistent with its integration into an instrument. Taking advantage of the
//  specific hardware architecture of the GPU, the COMPASS tool allows to
//  achieve adequate execution speeds to conduct large simulation campaigns
//  called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to
//  both testspecific components of AO of the E-ELT (such as wavefront analysis
//  device with a pyramid or elongated Laser star), and various systems
//  configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//  details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with COMPASS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      carma_magma.h
//! \ingroup   libcarma
//! \brief     this file provides wrappers to the magma functions
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _CARMA_MAGMA_H_
#define _CARMA_MAGMA_H_

#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>

// MAGMA functions
int carma_magma_disabled();
// template <class T>
// int carma_svd(CarmaObj<T> *imat, CarmaObj<T> *eigenvals,
//               CarmaObj<T> *mod2act, CarmaObj<T> *mes2mod);
template <class T>
int carma_magma_syevd(char jobz, CarmaObj<T> *mat,
                      CarmaHostObj<T> *eigenvals);
// template <class T, int method>
// int carma_magma_syevd(char jobz, CarmaObj<T> *mat, CarmaHostObj<T>
// *eigenvals);
template <class T>
int carma_magma_syevd_m(long ngpu, char jobz, long N, T *mat, T *eigenvals);
template <class T>
int carma_magma_syevd_m(long ngpu, char jobz, CarmaHostObj<T> *mat,
                        CarmaHostObj<T> *eigenvals);
template <class T>
int carma_magma_syevd_m(long ngpu, char jobz, CarmaHostObj<T> *mat,
                        CarmaHostObj<T> *eigenvals, CarmaHostObj<T> *U);
template <class T>
int carma_magma_getri(CarmaObj<T> *d_iA);
template <class T>
int carma_magma_potr_inv(CarmaObj<T> *d_iA);
template <class T>
int carma_magma_potr_inv_m(long num_gpus, CarmaHostObj<T> *h_A,
                        CarmaObj<T> *d_iA);

// MAGMA functions (direct access)
template <class T>
int carma_magma_syevd(char jobz, long N, T *mat, T *eigenvals);
// template <class T>
// int carma_magma_potr_inv_m(long num_gpus, long N, T *h_A, T *d_iA);

template <class T_data>
int carma_magma_svd_cpu(CarmaHostObj<T_data> *imat,
                        CarmaHostObj<T_data> *eigenvals,
                        CarmaHostObj<T_data> *mod2act,
                        CarmaHostObj<T_data> *mes2mod);
template <class T>
int carma_magma_getri_cpu(CarmaHostObj<T> *h_A);
template <class T>
int carma_magma_potr_inv_cpu(CarmaHostObj<T> *h_A);
template <class T>
int carma_magma_syevd_cpu(char jobz, CarmaHostObj<T> *h_A,
                          CarmaHostObj<T> *eigenvals);

// MAGMA functions (direct access)
// template <class T>
// int carma_magma_svd_cpu(long N, long M, T *imat, T *eigenvals, T *mod2act,
//                   T *mes2mod);
template <class T>
int carma_magma_getri_cpu(long N, T *h_A);
template <class T>
int carma_magma_potr_inv_cpu(long N, T *h_A);
template <class T>
int carma_magma_syevd_cpu(char jobz, long N, T *h_A, T *eigenvals);
template <class T>
int carma_magma_axpy_cpu(long N, T alpha, T *h_X, long incX, T *h_Y, long incY);
template <class T>
int carma_gemm_cpu(char transa, char transb, long m, long n, long k, T alpha,
                   T *A, long lda, T *B, long ldb, T beta, T *C, long ldc);

template <class T_data>
int carma_magma_gemv(char trans, int m, int n, T_data alpha, T_data *matA,
                     int lda, T_data *vectx, int incx, T_data beta,
                     T_data *vecty, int incy);

template <class T_data>
int carma_magma_csr2ell(CarmaSparseObj<T_data> *dA);

template <class T_data>
int carma_magma_spmv(T_data alpha, CarmaSparseObj<T_data> *dA,
                     CarmaObj<T_data> *dx, T_data beta, CarmaObj<T_data> *dy);

template <class T_data>
int carma_magma_sparse_free(CarmaSparseObj<T_data> *dA);

#endif  // _CARMA_MAGMA_H_
