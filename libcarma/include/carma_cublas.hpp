// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma_cublas.hpp
//! \ingroup   libcarma
//! \brief     this file provides the cublas features to CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef CARMA_CUBLAS_H_
#define CARMA_CUBLAS_H_

#include <cuda_runtime_api.h>
/* Using updated (v2) interfaces to cublas */
#include <cublas_v2.h>
#include <string>

#define carma_checkCublasStatus(status) \
  __carma_checkCublasStatus(status, __LINE__, __FILE__)

cublasStatus_t __carma_checkCublasStatus(cublasStatus_t status, int32_t line,
                                         std::string file);

cublasStatus_t carma_init_cublas(cublasHandle_t *cublas_handle);
cublasStatus_t carma_shutdown_cublas(cublasHandle_t cublas_handle);

cublasOperation_t carma_char2cublas_operation(char operation);

/*
 * _____ _____ __  __ ____  _        _  _____ _____ ____
 *|_   _| ____|  \/  |  _ \| |      / \|_   _| ____/ ___|
 *  | | |  _| | |\/| | |_) | |     / _ \ | | |  _| \___ \
 *  | | | |___| |  | |  __/| |___ / ___ \| | | |___ ___) |
 *  |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|____/
 *
 */

template <class T_data>
int32_t carma_where_amax(cublasHandle_t cublas_handle, int32_t n, const T_data *vect,
                     int32_t incx);

template <class T_data>
int32_t carma_where_amin(cublasHandle_t cublas_handle, int32_t n, const T_data *vect,
                     int32_t incx);

template <class T_data>
T_data carma_getasum(cublasHandle_t cublas_handle, int32_t n, const T_data *vect,
                     int32_t incx);

template <class T_data>
cublasStatus_t carma_axpy(cublasHandle_t cublas_handle, int32_t n,
                          const T_data alpha, const T_data *vectx, int32_t incx,
                          T_data *vecty, int32_t incy);

template <class T_data>
T_data carma_dot(cublasHandle_t cublas_handle, int32_t n, T_data *vectx, int32_t incx,
                 T_data *vecty, int32_t incy);

template <class T_data>
T_data carma_nrm2(cublasHandle_t cublas_handle, int32_t n, T_data *vect, int32_t incx);

template <class T_data>
cublasStatus_t carma_rot(cublasHandle_t cublas_handle, int32_t n, T_data *vectx,
                         int32_t incx, T_data *vecty, int32_t incy, T_data sc,
                         T_data ss);

template <class T_data>
cublasStatus_t carma_scal(cublasHandle_t cublas_handle, int32_t n, T_data alpha,
                          T_data *vectx, int32_t incx);

template <class T_data>
cublasStatus_t carma_swap(cublasHandle_t cublas_handle, int32_t n, T_data *vectx,
                          int32_t incx, T_data *vecty, int32_t incy);

template <class T_data>
cublasStatus_t carma_copy(cublasHandle_t cublas_handle, int32_t n,
                          const T_data *vectx, int32_t incx, T_data *vecty,
                          int32_t incy);

template <class T_data>
cublasStatus_t carma_gemv(cublasHandle_t cublas_handle, char trans, int32_t m,
                          int32_t n, T_data alpha, T_data *matA, int32_t lda,
                          T_data *vectx, int32_t incx, T_data beta, T_data *vecty,
                          int32_t incy);

template <class T_data>
cublasStatus_t carma_symv(cublasHandle_t cublas_handle, char uplo, int32_t n,
                          T_data alpha, T_data *matA, int32_t lda, T_data *vectx,
                          int32_t incx, T_data beta, T_data *vecty, int32_t incy);

template <class T_data>
cublasStatus_t carma_ger(cublasHandle_t cublas_handle, int32_t m, int32_t n,
                         T_data alpha, T_data *vectx, int32_t incx, T_data *vecty,
                         int32_t incy, T_data *matA, int32_t lda);

template <class T_data>
cublasStatus_t carma_gemm(cublasHandle_t cublas_handle, char transa,
                          char transb, int32_t m, int32_t n, int32_t k, T_data alpha,
                          T_data *matA, int32_t lda, T_data *matB, int32_t ldb,
                          T_data beta, T_data *matC, int32_t ldc);

template <class T_data>
cublasStatus_t carma_gemm_strided_batched(cublasHandle_t cublas_handle,
                    char transa, char transb, int32_t m, int32_t n, int32_t k, T_data alpha,
                    T_data *matsA, int32_t lda, int64_t strideA,
                    T_data *matsB, int32_t ldb, int64_t strideB,
                    T_data beta, T_data *matsC, int32_t ldc, int64_t strideC,
                    int32_t batch_count);

template <class T_data>
cublasStatus_t carma_symm(cublasHandle_t cublas_handle, char side, char uplo,
                          int32_t m, int32_t n, T_data alpha, T_data *matA, int32_t lda,
                          T_data *matB, int32_t ldb, T_data beta, T_data *matC,
                          int32_t ldc);

template <class T_data>
cublasStatus_t carma_syrk(cublasHandle_t cublas_handle, char uplo, char transa,
                          int32_t n, int32_t k, T_data alpha, T_data *matA, int32_t lda,
                          T_data beta, T_data *matC, int32_t ldc);

template <class T_data>
cublasStatus_t carma_syrkx(cublasHandle_t cublas_handle, char uplo, char transa,
                           int32_t n, int32_t k, T_data alpha, T_data *matA, int32_t lda,
                           T_data *matB, int32_t ldb, T_data beta, T_data *matC,
                           int32_t ldc);

template <class T_data>
cublasStatus_t carma_geam(cublasHandle_t cublas_handle, char transa,
                          char transb, int32_t m, int32_t n, T_data alpha, T_data *matA,
                          int32_t lda, T_data beta, T_data *matB, int32_t ldb,
                          T_data *matC, int32_t ldc);

template <class T_data>
cublasStatus_t carma_dgmm(cublasHandle_t cublas_handle, char side, int32_t m, int32_t n,
                          const T_data *matA, int32_t lda, const T_data *vectx,
                          int32_t incx, T_data *matC, int32_t ldc);

template <class T_data>
cublasStatus_t carma_sbmv(cublasHandle_t cublas_handle, char uplo, int32_t N, int32_t K, T_data alpha, T_data *A,
                          int32_t lda, T_data *x, int32_t incx, T_data beta, T_data *y, int32_t incy);
#endif /* CARMA_CUBLAS_H_ */