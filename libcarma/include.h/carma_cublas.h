// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser 
//  General Public License as published by the Free Software Foundation, either version 3 of the License, 
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration 
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems. 
//  
//  The final product includes a software package for simulating all the critical subcomponents of AO, 
//  particularly in the context of the ELT and a real-time core based on several control approaches, 
//  with performances consistent with its integration into an instrument. Taking advantage of the specific 
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT. 
//  
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components 
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and 
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      carma_cublas.h
//! \ingroup   libcarma
//! \brief     this file provides the cublas features to carma_obj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#ifndef CARMA_CUBLAS_H_
#define CARMA_CUBLAS_H_

#include <cuda_runtime_api.h>
/* Using updated (v2) interfaces to cublas */
#include <cublas_v2.h>
#include <string>

#define carma_checkCublasStatus(status) \
  __carma_checkCublasStatus(status, __LINE__, __FILE__)

cublasStatus_t __carma_checkCublasStatus(cublasStatus_t status, int line,
                                         std::string file);

cublasStatus_t carma_initCublas(cublasHandle_t *cublas_handle);
cublasStatus_t carma_shutdownCublas(cublasHandle_t cublas_handle);

cublasOperation_t carma_char2cublasOperation(char operation);

/*
 * _____ _____ __  __ ____  _        _  _____ _____ ____
 *|_   _| ____|  \/  |  _ \| |      / \|_   _| ____/ ___|
 *  | | |  _| | |\/| | |_) | |     / _ \ | | |  _| \___ \
 *  | | | |___| |  | |  __/| |___ / ___ \| | | |___ ___) |
 *  |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|____/
 *
 */

template <class T_data>
int carma_where_amax(cublasHandle_t cublas_handle, int n, const T_data *vect,
                     int incx);

template <class T_data>
int carma_where_amin(cublasHandle_t cublas_handle, int n, const T_data *vect,
                     int incx);

template <class T_data>
T_data carma_getasum(cublasHandle_t cublas_handle, int n, const T_data *vect,
                     int incx);

template <class T_data>
cublasStatus_t carma_axpy(cublasHandle_t cublas_handle, int n,
                          const T_data alpha, const T_data *vectx, int incx,
                          T_data *vecty, int incy);

template <class T_data>
T_data carma_dot(cublasHandle_t cublas_handle, int n, T_data *vectx, int incx,
                 T_data *vecty, int incy);

template <class T_data>
T_data carma_nrm2(cublasHandle_t cublas_handle, int n, T_data *vect, int incx);

template <class T_data>
cublasStatus_t carma_rot(cublasHandle_t cublas_handle, int n, T_data *vectx,
                         int incx, T_data *vecty, int incy, T_data sc,
                         T_data ss);

template <class T_data>
cublasStatus_t carma_scal(cublasHandle_t cublas_handle, int n, T_data alpha,
                          T_data *vectx, int incx);

template <class T_data>
cublasStatus_t carma_swap(cublasHandle_t cublas_handle, int n, T_data *vectx,
                          int incx, T_data *vecty, int incy);

template <class T_data>
cublasStatus_t carma_copy(cublasHandle_t cublas_handle, int n,
                          const T_data *vectx, int incx, T_data *vecty,
                          int incy);

template <class T_data>
cublasStatus_t carma_gemv(cublasHandle_t cublas_handle, char trans, int m,
                          int n, T_data alpha, T_data *matA, int lda,
                          T_data *vectx, int incx, T_data beta, T_data *vecty,
                          int incy);

template <class T_data>
cublasStatus_t carma_symv(cublasHandle_t cublas_handle, char uplo, int n,
                          T_data alpha, T_data *matA, int lda, T_data *vectx,
                          int incx, T_data beta, T_data *vecty, int incy);

template <class T_data>
cublasStatus_t carma_ger(cublasHandle_t cublas_handle, int m, int n,
                         T_data alpha, T_data *vectx, int incx, T_data *vecty,
                         int incy, T_data *matA, int lda);

template <class T_data>
cublasStatus_t carma_gemm(cublasHandle_t cublas_handle, char transa,
                          char transb, int m, int n, int k, T_data alpha,
                          T_data *matA, int lda, T_data *matB, int ldb,
                          T_data beta, T_data *matC, int ldc);

template <class T_data>
cublasStatus_t carma_symm(cublasHandle_t cublas_handle, char side, char uplo,
                          int m, int n, T_data alpha, T_data *matA, int lda,
                          T_data *matB, int ldb, T_data beta, T_data *matC,
                          int ldc);

template <class T_data>
cublasStatus_t carma_syrk(cublasHandle_t cublas_handle, char uplo, char transa,
                          int n, int k, T_data alpha, T_data *matA, int lda,
                          T_data beta, T_data *matC, int ldc);

template <class T_data>
cublasStatus_t carma_syrkx(cublasHandle_t cublas_handle, char uplo, char transa,
                           int n, int k, T_data alpha, T_data *matA, int lda,
                           T_data *matB, int ldb, T_data beta, T_data *matC,
                           int ldc);

template <class T_data>
cublasStatus_t carma_geam(cublasHandle_t cublas_handle, char transa,
                          char transb, int m, int n, T_data alpha, T_data *matA,
                          int lda, T_data beta, T_data *matB, int ldb,
                          T_data *matC, int ldc);

template <class T_data>
cublasStatus_t carma_dgmm(cublasHandle_t cublas_handle, char side, int m, int n,
                          const T_data *matA, int lda, const T_data *vectx,
                          int incx, T_data *matC, int ldc);

#endif /* CARMA_CUBLAS_H_ */
