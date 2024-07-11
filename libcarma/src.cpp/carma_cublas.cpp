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

//! \file      carma_cublas.cpp
//! \ingroup   libcarma
//! \brief     this file provides the cublas features to CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include <carma_cublas.hpp>
#include <carma_utils.hpp>
#include <iostream>
#include <stdexcept>
#include <string>

#include <type_list.hpp>

#ifdef CAN_DO_HALF
using TypeListObj = GenericTypeList<int32_t, uint32_t, uint16_t, float, double,
                                    half, cuFloatComplex,
                                    cuDoubleComplex>;  // , tuple_t<float>>;
#else
using TypeListObj =
    GenericTypeList<int32_t, uint32_t, uint16_t, float, double, cuFloatComplex,
                    cuDoubleComplex>;  // , tuple_t<float>>;
#endif

#define CARMA_NIY                                            \
  {                                                          \
    DEBUG_TRACE("Method not implemented yet!");              \
    throw std::runtime_error("Method not implemented yet!"); \
  }

cublasStatus_t __carma_checkCublasStatus(cublasStatus_t status, int32_t line,
                                         std::string file) {
  /**< Generic CUBLAS check status routine */
  switch (status) {
    case CUBLAS_STATUS_SUCCESS:
      return status;
    case CUBLAS_STATUS_NOT_INITIALIZED:
      std::cerr << "CUBLAS error : Unsupported numerical value was passed to "
                   "function."
                << std::endl;
      break;
    case CUBLAS_STATUS_ALLOC_FAILED:
      std::cerr << "CUBLAS error : Resource allocation failed." << std::endl;
      break;
    case CUBLAS_STATUS_INVALID_VALUE:
      std::cerr << "CUBLAS error : CUBLAS_STATUS_ALLOC_FAILED !!!!!"
                << std::endl;
      break;
    case CUBLAS_STATUS_ARCH_MISMATCH:
      std::cerr << "CUBLAS error : Resource allocation failed." << std::endl;
      break;
    case CUBLAS_STATUS_MAPPING_ERROR:
      std::cerr << "CUBLAS error : Unsupported numerical value was passed to "
                   "function."
                << std::endl;
      break;
    case CUBLAS_STATUS_EXECUTION_FAILED:
      std::cerr << "CUBLAS error : Resource allocation failed." << std::endl;
      break;
    case CUBLAS_STATUS_INTERNAL_ERROR:
      std::cerr << "!CUBLAS error : An internal CUBLAS operation failed."
                << std::endl;
      break;
// Define not supported status for pre-6.0 compatibility.
#if CUDA_VERSION > 6000
    case CUBLAS_STATUS_NOT_SUPPORTED:
      std::cerr << "CUBLAS error : Unsupported numerical value was passed to "
                   "function."
                << std::endl;
      break;
#endif
#if CUDA_VERSION > 6000
    case CUBLAS_STATUS_LICENSE_ERROR:
      std::cerr << "CUBLAS error : The functionnality requested requires some "
                   "license and an error was detected when trying to check the "
                   "current licensing. This error can happen if the license is "
                   "not present or is expired or if the environment variable "
                   "NVIDIA_LICENSE_FILE is not set properly."
                << std::endl;
#endif
  }
  std::cerr << "CUBLAS error in " << file << "@" << line << std::endl;
  return status;
}

cublasStatus_t carma_init_cublas(cublasHandle_t *cublas_handle) {
  /**< Generic CUBLAS init routine */
  return carma_checkCublasStatus(cublasCreate(cublas_handle));
}

cublasStatus_t carma_shutdown_cublas(cublasHandle_t cublas_handle) {
  /**< Generic CUBLAS shutdown routine */
  return carma_checkCublasStatus(cublasDestroy(cublas_handle));
}

cublasOperation_t carma_char2cublas_operation(char operation) {
  switch (operation) {
    case 't':
    case 'T':
      return CUBLAS_OP_T;
    case 'c':
    case 'C':
      return CUBLAS_OP_C;
    default:
      return CUBLAS_OP_N;
  }
}

cublasSideMode_t carma_char2cublasSide(char operation) {
  switch (operation) {
    case 'l':
    case 'L':
      return CUBLAS_SIDE_LEFT;
    default:
      return CUBLAS_SIDE_RIGHT;
  }
}

cublasFillMode_t carma_char2cublasFillMode(char operation) {
  switch (operation) {
    case 'l':
    case 'L':
      return CUBLAS_FILL_MODE_LOWER;
    default:
      return CUBLAS_FILL_MODE_UPPER;
  }
}
/*
 *  _____ _____ __  __ ____  _        _  _____ _____ ____
 * |_   _| ____|  \/  |  _ \| |      / \|_   _| ____/ ___|
 *   | | |  _| | |\/| | |_) | |     / _ \ | | |  _| \___ \
 *   | | | |___| |  | |  __/| |___ / ___ \| | | |___ ___) |
 *   |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|____/
 *
 */

/** These templates are used to select the proper Iamax and Iamin executable
 * from T_data*/
template <class T_data,
          cublasStatus_t (*afunc)(cublasHandle_t handle, int32_t n, const T_data *x,
                                  int32_t incx, int32_t *result)>
int32_t carma_where_gen(cublasHandle_t cublas_handle, int32_t n, const T_data *vect,
                    int32_t incx) {
  int32_t result = 0;
  carma_checkCublasStatus(afunc(cublas_handle, n, vect, incx, &result));
  return result - 1;
}

template <class T>
int32_t carma_where_amax(cublasHandle_t cublas_handle, int32_t n, const T *vect,
                     int32_t incx) CARMA_NIY;
/**< Specialized template for carma_amax executable selection */
template <>
int32_t carma_where_amax(cublasHandle_t cublas_handle, int32_t n, const float *vect,
                     int32_t incx) {
  return carma_where_gen<float, cublasIsamax>(cublas_handle, n, vect, incx);
}
template <>
int32_t carma_where_amax(cublasHandle_t cublas_handle, int32_t n, const double *vect,
                     int32_t incx) {
  return carma_where_gen<double, cublasIdamax>(cublas_handle, n, vect, incx);
}
template <>
int32_t carma_where_amax(cublasHandle_t cublas_handle, int32_t n,
                     const cuFloatComplex *vect, int32_t incx) {
  return carma_where_gen<cuFloatComplex, cublasIcamax>(cublas_handle, n, vect,
                                                       incx);
}
template <>
int32_t carma_where_amax(cublasHandle_t cublas_handle, int32_t n,
                     const cuDoubleComplex *vect, int32_t incx) {
  return carma_where_gen<cuDoubleComplex, cublasIzamax>(cublas_handle, n, vect,
                                                        incx);
}

template <class T>
int32_t carma_where_amin(cublasHandle_t cublas_handle, int32_t n, const T *vect,
                     int32_t incx) CARMA_NIY;
/**< Specialized template for carma_amin executable selection */
template <>
int32_t carma_where_amin(cublasHandle_t cublas_handle, int32_t n, const float *vect,
                     int32_t incx) {
  return carma_where_gen<float, cublasIsamin>(cublas_handle, n, vect, incx);
}
template <>
int32_t carma_where_amin(cublasHandle_t cublas_handle, int32_t n, const double *vect,
                     int32_t incx) {
  return carma_where_gen<double, cublasIdamin>(cublas_handle, n, vect, incx);
}
template <>
int32_t carma_where_amin(cublasHandle_t cublas_handle, int32_t n,
                     const cuFloatComplex *vect, int32_t incx) {
  return carma_where_gen<cuFloatComplex, cublasIcamin>(cublas_handle, n, vect,
                                                       incx);
}
template <>
int32_t carma_where_amin(cublasHandle_t cublas_handle, int32_t n,
                     const cuDoubleComplex *vect, int32_t incx) {
  return carma_where_gen<cuDoubleComplex, cublasIzamin>(cublas_handle, n, vect,
                                                        incx);
}

/** These templates are used to select the proper asum executable from T_data*/
template <class T_data,
          cublasStatus_t (*afunc)(cublasHandle_t handle, int32_t n, const T_data *x,
                                  int32_t incx, T_data *result)>
T_data carma_asum_gen(cublasHandle_t cublas_handle, int32_t n, const T_data *vect,
                      int32_t incx) {
  T_data result = 0;
  carma_checkCublasStatus(afunc(cublas_handle, n, vect, incx, &result));
  return result;
}

template <class T>
T carma_getasum(cublasHandle_t cublas_handle, int32_t n, const T *vect,
                int32_t incx) CARMA_NIY;
/**< Specialized template for carma_getasum executable selection */
template <>
float carma_getasum(cublasHandle_t cublas_handle, int32_t n, const float *vect,
                    int32_t incx) {
  return carma_asum_gen<float, cublasSasum>(cublas_handle, n, vect, incx);
}
template <>
double carma_getasum(cublasHandle_t cublas_handle, int32_t n, const double *vect,
                     int32_t incx) {
  return carma_asum_gen<double, cublasDasum>(cublas_handle, n, vect, incx);
}

/** These templates are used to select the proper axpy executable from T_data*/
template <class T_data>
cublasStatus_t carma_axpy(cublasHandle_t cublas_handle, int32_t n,
                          const T_data alpha, const T_data *vectx, int32_t incx,
                          T_data *vecty, int32_t incy) CARMA_NIY;
/**< Specialized template for carma_axpy executable selection */
template <>
cublasStatus_t carma_axpy<float>(cublasHandle_t cublas_handle, int32_t n,
                                 const float alpha, const float *vectx,
                                 int32_t incx, float *vecty, int32_t incy) {
  return carma_checkCublasStatus(
      cublasSaxpy(cublas_handle, n, &alpha, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_axpy<double>(cublasHandle_t cublas_handle, int32_t n,
                                  const double alpha, const double *vectx,
                                  int32_t incx, double *vecty, int32_t incy) {
  return carma_checkCublasStatus(
      cublasDaxpy(cublas_handle, n, &alpha, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_axpy<cuFloatComplex>(cublasHandle_t cublas_handle, int32_t n,
                                          const cuFloatComplex alpha,
                                          const cuFloatComplex *vectx, int32_t incx,
                                          cuFloatComplex *vecty, int32_t incy) {
  return carma_checkCublasStatus(
      cublasCaxpy(cublas_handle, n, &alpha, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_axpy<cuDoubleComplex>(cublasHandle_t cublas_handle, int32_t n,
                                           const cuDoubleComplex alpha,
                                           const cuDoubleComplex *vectx,
                                           int32_t incx, cuDoubleComplex *vecty,
                                           int32_t incy) {
  return carma_checkCublasStatus(
      cublasZaxpy(cublas_handle, n, &alpha, vectx, incx, vecty, incy));
}

/** These templates are used to select the proper copy
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_copy(cublasHandle_t cublas_handle, int32_t n,
                          const T_data *vectx, int32_t incx, T_data *vecty,
                          int32_t incy) CARMA_NIY;
/**< Specialized template for carma_copy executable selection */
template <>
cublasStatus_t carma_copy<float>(cublasHandle_t cublas_handle, int32_t n,
                                 const float *vectx, int32_t incx, float *vecty,
                                 int32_t incy) {
  return carma_checkCublasStatus(
      cublasScopy(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_copy<double>(cublasHandle_t cublas_handle, int32_t n,
                                  const double *vectx, int32_t incx, double *vecty,
                                  int32_t incy) {
  return carma_checkCublasStatus(
      cublasDcopy(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_copy<cuFloatComplex>(cublasHandle_t cublas_handle, int32_t n,
                                          const cuFloatComplex *vectx, int32_t incx,
                                          cuFloatComplex *vecty, int32_t incy) {
  return carma_checkCublasStatus(
      cublasCcopy(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_copy<cuDoubleComplex>(cublasHandle_t cublas_handle, int32_t n,
                                           const cuDoubleComplex *vectx,
                                           int32_t incx, cuDoubleComplex *vecty,
                                           int32_t incy) {
  return carma_checkCublasStatus(
      cublasZcopy(cublas_handle, n, vectx, incx, vecty, incy));
}

/** These templates are used to select the proper dot
 * executable from T_data*/
template <class T_data>
T_data carma_dot(cublasHandle_t cublas_handle, int32_t n, T_data *vectx, int32_t incx,
                 T_data *vecty, int32_t incy) CARMA_NIY;
/**< Specialized template for carma_dot executable selection */
template <>
float carma_dot<float>(cublasHandle_t cublas_handle, int32_t n, float *vectx,
                       int32_t incx, float *vecty, int32_t incy) {
  float result = 0;
  carma_checkCublasStatus(
      cublasSdot(cublas_handle, n, vectx, incx, vecty, incy, &result));
  return result;
}
template <>
double carma_dot<double>(cublasHandle_t cublas_handle, int32_t n, double *vectx,
                         int32_t incx, double *vecty, int32_t incy) {
  double result = 0;
  carma_checkCublasStatus(
      cublasDdot(cublas_handle, n, vectx, incx, vecty, incy, &result));
  return result;
}
template <>
cuFloatComplex carma_dot<cuFloatComplex>(cublasHandle_t cublas_handle, int32_t n,
                                         cuFloatComplex *vectx, int32_t incx,
                                         cuFloatComplex *vecty, int32_t incy) {
  cuFloatComplex result;
  result.x = 0;
  result.y = 0;
  carma_checkCublasStatus(
      cublasCdotu(cublas_handle, n, vectx, incx, vecty, incy, &result));
  return result;
}
template <>
cuDoubleComplex carma_dot<cuDoubleComplex>(cublasHandle_t cublas_handle, int32_t n,
                                           cuDoubleComplex *vectx, int32_t incx,
                                           cuDoubleComplex *vecty, int32_t incy) {
  cuDoubleComplex result;
  result.x = 0;
  result.y = 0;
  carma_checkCublasStatus(
      cublasZdotu(cublas_handle, n, vectx, incx, vecty, incy, &result));
  return result;
}

/** These templates are used to select the proper nrm2
 * executable from T_data*/
template <class T_data>
T_data carma_nrm2(cublasHandle_t cublas_handle, int32_t n, T_data *vect,
                  int32_t incx) CARMA_NIY;
/**< Specialized template for carma_nrm2 executable selection */
template <>
float carma_nrm2<float>(cublasHandle_t cublas_handle, int32_t n, float *vect,
                        int32_t incx) {
  float result = 0;
  carma_checkCublasStatus(cublasSnrm2(cublas_handle, n, vect, incx, &result));
  return result;
}
template <>
double carma_nrm2<double>(cublasHandle_t cublas_handle, int32_t n, double *vect,
                          int32_t incx) {
  double result = 0;
  carma_checkCublasStatus(cublasDnrm2(cublas_handle, n, vect, incx, &result));
  return result;
}
/*add
 template<>
 float carma_nrm2<cuFloatComplex>(cublasHandle_t
 cublas_handle, int32_t n, cuFloatComplex *vect, int32_t incx)
 { float result = 0;
 carma_checkCublasStatus(cublasScnrm2(cublas_handle,n,
 vect,incx,&result)); return result;
 }*/

/** These templates are used to select the proper rot
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_rot(cublasHandle_t cublas_handle, int32_t n, T_data *vectx,
                         int32_t incx, T_data *vecty, int32_t incy, T_data sc,
                         T_data ss) CARMA_NIY;
/**< Specialized template for carma_rot executable selection */
template <>
cublasStatus_t carma_rot<float>(cublasHandle_t cublas_handle, int32_t n,
                                float *vectx, int32_t incx, float *vecty, int32_t incy,
                                float sc, float ss) {
  return carma_checkCublasStatus(
      cublasSrot(cublas_handle, n, vectx, incx, vecty, incy, &sc, &ss));
}
template <>
cublasStatus_t carma_rot<double>(cublasHandle_t cublas_handle, int32_t n,
                                 double *vectx, int32_t incx, double *vecty,
                                 int32_t incy, double sc, double ss) {
  return carma_checkCublasStatus(
      cublasDrot(cublas_handle, n, vectx, incx, vecty, incy, &sc, &ss));
}

/** These templates are used to select the proper scal
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_scal(cublasHandle_t cublas_handle, int32_t n, T_data alpha,
                          T_data *vectx, int32_t incx) CARMA_NIY;
/**< Specialized template for carma_scal executable selection */
template <>
cublasStatus_t carma_scal<float>(cublasHandle_t cublas_handle, int32_t n,
                                 float alpha, float *vectx, int32_t incx) {
  return carma_checkCublasStatus(
      cublasSscal(cublas_handle, n, &alpha, vectx, incx));
}
template <>
cublasStatus_t carma_scal<double>(cublasHandle_t cublas_handle, int32_t n,
                                  double alpha, double *vectx, int32_t incx) {
  return carma_checkCublasStatus(
      cublasDscal(cublas_handle, n, &alpha, vectx, incx));
}
template <>
cublasStatus_t carma_scal<cuFloatComplex>(cublasHandle_t cublas_handle, int32_t n,
                                          cuFloatComplex alpha,
                                          cuFloatComplex *vectx, int32_t incx) {
  return carma_checkCublasStatus(
      cublasCscal(cublas_handle, n, &alpha, vectx, incx));
}
template <>
cublasStatus_t carma_scal<cuDoubleComplex>(cublasHandle_t cublas_handle, int32_t n,
                                           cuDoubleComplex alpha,
                                           cuDoubleComplex *vectx, int32_t incx) {
  return carma_checkCublasStatus(
      cublasZscal(cublas_handle, n, &alpha, vectx, incx));
}

/** These templates are used to select the proper swap
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_swap(cublasHandle_t cublas_handle, int32_t n, T_data *vectx,
                          int32_t incx, T_data *vecty, int32_t incy) CARMA_NIY;
/**< Specialized template for carma_swap executable selection */
template <>
cublasStatus_t carma_swap<float>(cublasHandle_t cublas_handle, int32_t n,
                                 float *vectx, int32_t incx, float *vecty,
                                 int32_t incy) {
  return carma_checkCublasStatus(
      cublasSswap(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_swap<double>(cublasHandle_t cublas_handle, int32_t n,
                                  double *vectx, int32_t incx, double *vecty,
                                  int32_t incy) {
  return carma_checkCublasStatus(
      cublasDswap(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_swap<cuFloatComplex>(cublasHandle_t cublas_handle, int32_t n,
                                          cuFloatComplex *vectx, int32_t incx,
                                          cuFloatComplex *vecty, int32_t incy) {
  return carma_checkCublasStatus(
      cublasCswap(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_swap<cuDoubleComplex>(cublasHandle_t cublas_handle, int32_t n,
                                           cuDoubleComplex *vectx, int32_t incx,
                                           cuDoubleComplex *vecty, int32_t incy) {
  return carma_checkCublasStatus(
      cublasZswap(cublas_handle, n, vectx, incx, vecty, incy));
}

/** These templates are used to select the proper gemv
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_gemv(cublasHandle_t cublas_handle, char trans, int32_t m,
                          int32_t n, T_data alpha, T_data *matA, int32_t lda,
                          T_data *vectx, int32_t incx, T_data beta, T_data *vecty,
                          int32_t incy) CARMA_NIY;
/**< Specialized template for carma_gemv executable selection */
template <>
cublasStatus_t carma_gemv<float>(cublasHandle_t cublas_handle, char trans,
                                 int32_t m, int32_t n, float alpha, float *matA,
                                 int32_t lda, float *vectx, int32_t incx, float beta,
                                 float *vecty, int32_t incy) {
  cublasOperation_t trans2 = carma_char2cublas_operation(trans);
  return carma_checkCublasStatus(cublasSgemv(cublas_handle, trans2, m, n,
                                             &alpha, matA, lda, vectx, incx,
                                             &beta, vecty, incy));
}
template <>
cublasStatus_t carma_gemv<double>(cublasHandle_t cublas_handle, char trans,
                                  int32_t m, int32_t n, double alpha, double *matA,
                                  int32_t lda, double *vectx, int32_t incx, double beta,
                                  double *vecty, int32_t incy) {
  cublasOperation_t trans2 = carma_char2cublas_operation(trans);
  return carma_checkCublasStatus(cublasDgemv(cublas_handle, trans2, m, n,
                                             &alpha, matA, lda, vectx, incx,
                                             &beta, vecty, incy));
}
template <>
cublasStatus_t carma_gemv<cuFloatComplex>(
    cublasHandle_t cublas_handle, char trans, int32_t m, int32_t n,
    cuFloatComplex alpha, cuFloatComplex *matA, int32_t lda, cuFloatComplex *vectx,
    int32_t incx, cuFloatComplex beta, cuFloatComplex *vecty, int32_t incy) {
  cublasOperation_t trans2 = carma_char2cublas_operation(trans);
  return carma_checkCublasStatus(cublasCgemv(cublas_handle, trans2, m, n,
                                             &alpha, matA, lda, vectx, incx,
                                             &beta, vecty, incy));
}
template <>
cublasStatus_t carma_gemv<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char trans, int32_t m, int32_t n,
                                           cuDoubleComplex alpha,
                                           cuDoubleComplex *matA, int32_t lda,
                                           cuDoubleComplex *vectx, int32_t incx,
                                           cuDoubleComplex beta,
                                           cuDoubleComplex *vecty, int32_t incy) {
  cublasOperation_t trans2 = carma_char2cublas_operation(trans);
  return carma_checkCublasStatus(cublasZgemv(cublas_handle, trans2, m, n,
                                             &alpha, matA, lda, vectx, incx,
                                             &beta, vecty, incy));
}
#ifdef CAN_DO_HALF
template <>
cublasStatus_t carma_gemv<half>(cublasHandle_t cublas_handle, char trans, int32_t m,
                                int32_t n, half alpha, half *matA, int32_t lda,
                                half *vectx, int32_t incx, half beta, half *vecty,
                                int32_t incy) {
  int32_t k = (((trans == 'N') || (trans == 'n')) ? n : m);
  cublasOperation_t transa2 = carma_char2cublas_operation(trans);
  cublasOperation_t transb2 = carma_char2cublas_operation('N');
  return carma_checkCublasStatus(cublasHgemm(cublas_handle, transa2, transb2, m,
                                             1, k, &alpha, matA, lda, vectx, n,
                                             &beta, vecty, n));
}
#endif
/** These templates are used to select the proper ger executable from T_data*/
template <class T_data>
cublasStatus_t carma_ger(cublasHandle_t cublas_handle, int32_t m, int32_t n,
                         T_data alpha, T_data *vectx, int32_t incx, T_data *vecty,
                         int32_t incy, T_data *matA, int32_t lda) CARMA_NIY;
/**< Specialized template for carma_ger executable selection */
template <>
cublasStatus_t carma_ger<float>(cublasHandle_t cublas_handle, int32_t m, int32_t n,
                                float alpha, float *vectx, int32_t incx,
                                float *vecty, int32_t incy, float *matA, int32_t lda) {
  return carma_checkCublasStatus(cublasSger(cublas_handle, m, n, &alpha, vectx,
                                            incx, vecty, incy, matA, lda));
}
template <>
cublasStatus_t carma_ger<double>(cublasHandle_t cublas_handle, int32_t m, int32_t n,
                                 double alpha, double *vectx, int32_t incx,
                                 double *vecty, int32_t incy, double *matA,
                                 int32_t lda) {
  return carma_checkCublasStatus(cublasDger(cublas_handle, m, n, &alpha, vectx,
                                            incx, vecty, incy, matA, lda));
}
template <>
cublasStatus_t carma_ger<cuFloatComplex>(cublasHandle_t cublas_handle, int32_t m,
                                         int32_t n, cuFloatComplex alpha,
                                         cuFloatComplex *vectx, int32_t incx,
                                         cuFloatComplex *vecty, int32_t incy,
                                         cuFloatComplex *matA, int32_t lda) {
  return carma_checkCublasStatus(cublasCgeru(cublas_handle, m, n, &alpha, vectx,
                                             incx, vecty, incy, matA, lda));
}
template <>
cublasStatus_t carma_ger<cuDoubleComplex>(cublasHandle_t cublas_handle, int32_t m,
                                          int32_t n, cuDoubleComplex alpha,
                                          cuDoubleComplex *vectx, int32_t incx,
                                          cuDoubleComplex *vecty, int32_t incy,
                                          cuDoubleComplex *matA, int32_t lda) {
  return carma_checkCublasStatus(cublasZgeru(cublas_handle, m, n, &alpha, vectx,
                                             incx, vecty, incy, matA, lda));
}

/** These templates are used to select the proper symv
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_symv(cublasHandle_t cublas_handle, char uplo, int32_t n,
                          T_data alpha, T_data *matA, int32_t lda, T_data *vectx,
                          int32_t incx, T_data beta, T_data *vecty,
                          int32_t incy) CARMA_NIY;
/**< Specialized template for carma_symv executable selection */
template <>
cublasStatus_t carma_symv<float>(cublasHandle_t cublas_handle, char uplo, int32_t n,
                                 float alpha, float *matA, int32_t lda,
                                 float *vectx, int32_t incx, float beta,
                                 float *vecty, int32_t incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);

  return carma_checkCublasStatus(cublasSsymv(cublas_handle, filla, n, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}
template <>
cublasStatus_t carma_symv<double>(cublasHandle_t cublas_handle, char uplo,
                                  int32_t n, double alpha, double *matA, int32_t lda,
                                  double *vectx, int32_t incx, double beta,
                                  double *vecty, int32_t incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasDsymv(cublas_handle, filla, n, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}
template <>
cublasStatus_t carma_symv<cuFloatComplex>(
    cublasHandle_t cublas_handle, char uplo, int32_t n, cuFloatComplex alpha,
    cuFloatComplex *matA, int32_t lda, cuFloatComplex *vectx, int32_t incx,
    cuFloatComplex beta, cuFloatComplex *vecty, int32_t incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasCsymv(cublas_handle, filla, n, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}
template <>
cublasStatus_t carma_symv<cuDoubleComplex>(
    cublasHandle_t cublas_handle, char uplo, int32_t n, cuDoubleComplex alpha,
    cuDoubleComplex *matA, int32_t lda, cuDoubleComplex *vectx, int32_t incx,
    cuDoubleComplex beta, cuDoubleComplex *vecty, int32_t incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasZsymv(cublas_handle, filla, n, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}

/** These templates are used to select the proper gemm
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_gemm(cublasHandle_t cublas_handle, char transa,
                          char transb, int32_t m, int32_t n, int32_t k, T_data alpha,
                          T_data *matA, int32_t lda, T_data *matB, int32_t ldb,
                          T_data beta, T_data *matC, int32_t ldc) CARMA_NIY;
/**< Specialized template for carma_gemm executable selection */
template <>
cublasStatus_t carma_gemm<float>(cublasHandle_t cublas_handle, char transa,
                                 char transb, int32_t m, int32_t n, int32_t k, float alpha,
                                 float *matA, int32_t lda, float *matB, int32_t ldb,
                                 float beta, float *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasSgemm(cublas_handle, transa2, transb2, m,
                                             n, k, &alpha, matA, lda, matB, ldb,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_gemm<double>(cublasHandle_t cublas_handle, char transa,
                                  char transb, int32_t m, int32_t n, int32_t k,
                                  double alpha, double *matA, int32_t lda,
                                  double *matB, int32_t ldb, double beta,
                                  double *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasDgemm(cublas_handle, transa2, transb2, m,
                                             n, k, &alpha, matA, lda, matB, ldb,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_gemm<cuFloatComplex>(
    cublasHandle_t cublas_handle, char transa, char transb, int32_t m, int32_t n, int32_t k,
    cuFloatComplex alpha, cuFloatComplex *matA, int32_t lda, cuFloatComplex *matB,
    int32_t ldb, cuFloatComplex beta, cuFloatComplex *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasCgemm(cublas_handle, transa2, transb2, m,
                                             n, k, &alpha, matA, lda, matB, ldb,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_gemm<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char transa, char transb, int32_t m,
                                           int32_t n, int32_t k, cuDoubleComplex alpha,
                                           cuDoubleComplex *matA, int32_t lda,
                                           cuDoubleComplex *matB, int32_t ldb,
                                           cuDoubleComplex beta,
                                           cuDoubleComplex *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasZgemm(cublas_handle, transa2, transb2, m,
                                             n, k, &alpha, matA, lda, matB, ldb,
                                             &beta, matC, ldc));
}

#ifdef CAN_DO_HALF
template <>
cublasStatus_t carma_gemm<half>(cublasHandle_t cublas_handle, char transa,
                                char transb, int32_t m, int32_t n, int32_t k, half alpha,
                                half *matA, int32_t lda, half *matB, int32_t ldb,
                                half beta, half *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasHgemm(cublas_handle, transa2, transb2, m,
                                             n, k, &alpha, matA, lda, matB, ldb,
                                             &beta, matC, ldc));
}
#endif

/** These templates are used to select the proper batched gemm
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_gemm_strided_batched(cublasHandle_t cublas_handle,
                    char transa, char transb, int32_t m, int32_t n, int32_t k, T_data alpha,
                    T_data *matsA, int32_t lda, int64_t strideA,
                    T_data *matsB, int32_t ldb, int64_t strideB, T_data beta,
                    T_data *matsC, int32_t ldc, int64_t strideC,
                    int32_t batch_count) CARMA_NIY;
/**< Specialized template for carma_gemm_batched executable selection */
template <>
cublasStatus_t carma_gemm_strided_batched<float>(cublasHandle_t cublas_handle,
                    char transa, char transb, int32_t m, int32_t n, int32_t k, float alpha,
                    float *matsA, int32_t lda, int64_t strideA,
                    float *matsB, int32_t ldb, int64_t strideB, float beta,
                    float *matsC, int32_t ldc, int64_t strideC,
                    int32_t batch_count) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasSgemmStridedBatched(cublas_handle,
                      transa2, transb2, m, n, k, &alpha, matsA, lda, strideA,
                      matsB, ldb, strideB, &beta, matsC, ldc, strideC,
                      batch_count));
}
template <>
cublasStatus_t carma_gemm_strided_batched<double>(cublasHandle_t cublas_handle,
                    char transa, char transb, int32_t m, int32_t n, int32_t k, double alpha,
                    double *matsA, int32_t lda, int64_t strideA,
                    double *matsB, int32_t ldb, int64_t strideB, double beta,
                    double *matsC, int32_t ldc, int64_t strideC,
                    int32_t batch_count) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasDgemmStridedBatched(cublas_handle,
                      transa2, transb2, m, n, k, &alpha, matsA, lda, strideA,
                      matsB, ldb, strideB, &beta, matsC, ldc, strideC,
                      batch_count));
}
template <>
cublasStatus_t carma_gemm_strided_batched<cuFloatComplex>(
    cublasHandle_t cublas_handle, char transa, char transb, int32_t m, int32_t n, int32_t k,
    cuFloatComplex alpha, cuFloatComplex *matsA, int32_t lda, int64_t strideA,
    cuFloatComplex *matsB, int32_t ldb, int64_t strideB, cuFloatComplex beta,
    cuFloatComplex *matsC, int32_t ldc, int64_t strideC, int32_t batch_count) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasCgemmStridedBatched(cublas_handle,
                      transa2, transb2, m, n, k, &alpha, matsA, lda, strideA,
                      matsB, ldb, strideB, &beta, matsC, ldc, strideC,
                      batch_count));
}
template <>
cublasStatus_t carma_gemm_strided_batched<cuDoubleComplex>(
  cublasHandle_t cublas_handle, char transa, char transb, int32_t m, int32_t n, int32_t k,
  cuDoubleComplex alpha, cuDoubleComplex *matsA, int32_t lda, int64_t strideA,
  cuDoubleComplex *matsB, int32_t ldb, int64_t strideB, cuDoubleComplex beta,
  cuDoubleComplex *matsC, int32_t ldc, int64_t strideC, int32_t batch_count) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasZgemmStridedBatched(cublas_handle,
                      transa2, transb2, m, n, k, &alpha, matsA, lda, strideA,
                      matsB, ldb, strideB, &beta, matsC, ldc, strideC,
                      batch_count));
}

#ifdef CAN_DO_HALF
template <>
cublasStatus_t carma_gemm_strided_batched<half>(cublasHandle_t cublas_handle,
                      char transa, char transb, int32_t m, int32_t n, int32_t k, half alpha,
                      half *matsA, int32_t lda, int64_t strideA,
                      half *matsB, int32_t ldb, int64_t strideB, half beta,
                      half *matsC, int32_t ldc, int64_t strideC,
                      int32_t batch_count) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasHgemmStridedBatched(cublas_handle,
                      transa2, transb2, m, n, k, &alpha, matsA, lda, strideA,
                      matsB, ldb, strideB, &beta, matsC, ldc, strideC,
                      batch_count));
}
#endif

/** These templates are used to select the proper symm
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_symm(cublasHandle_t cublas_handle, char side, char uplo,
                          int32_t m, int32_t n, T_data alpha, T_data *matA, int32_t lda,
                          T_data *matB, int32_t ldb, T_data beta, T_data *matC,
                          int32_t ldc) CARMA_NIY;
/**< Specialized template for carma_symm executable selection */
template <>
cublasStatus_t carma_symm<float>(cublasHandle_t cublas_handle, char side,
                                 char uplo, int32_t m, int32_t n, float alpha,
                                 float *matA, int32_t lda, float *matB, int32_t ldb,
                                 float beta, float *matC, int32_t ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);

  return carma_checkCublasStatus(cublasSsymm(cublas_handle, side_cublas,
                                             uplo_cublas, m, n, &alpha, matA,
                                             lda, matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_symm<double>(cublasHandle_t cublas_handle, char side,
                                  char uplo, int32_t m, int32_t n, double alpha,
                                  double *matA, int32_t lda, double *matB, int32_t ldb,
                                  double beta, double *matC, int32_t ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasDsymm(cublas_handle, side_cublas,
                                             uplo_cublas, m, n, &alpha, matA,
                                             lda, matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_symm<cuFloatComplex>(
    cublasHandle_t cublas_handle, char side, char uplo, int32_t m, int32_t n,
    cuFloatComplex alpha, cuFloatComplex *matA, int32_t lda, cuFloatComplex *matB,
    int32_t ldb, cuFloatComplex beta, cuFloatComplex *matC, int32_t ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasCsymm(cublas_handle, side_cublas,
                                             uplo_cublas, m, n, &alpha, matA,
                                             lda, matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_symm<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char side, char uplo, int32_t m, int32_t n,
                                           cuDoubleComplex alpha,
                                           cuDoubleComplex *matA, int32_t lda,
                                           cuDoubleComplex *matB, int32_t ldb,
                                           cuDoubleComplex beta,
                                           cuDoubleComplex *matC, int32_t ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasZsymm(cublas_handle, side_cublas,
                                             uplo_cublas, m, n, &alpha, matA,
                                             lda, matB, ldb, &beta, matC, ldc));
}

/** These templates are used to select the proper syrk
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_syrk(cublasHandle_t cublas_handle, char uplo, char transa,
                          int32_t n, int32_t k, T_data alpha, T_data *matA, int32_t lda,
                          T_data beta, T_data *matC, int32_t ldc) CARMA_NIY;
/**< Specialized template for carma_syrk executable selection */
template <>
cublasStatus_t carma_syrk<float>(cublasHandle_t cublas_handle, char uplo,
                                 char transa, int32_t n, int32_t k, float alpha,
                                 float *matA, int32_t lda, float beta, float *matC,
                                 int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasSsyrk(cublas_handle, uplo_cublas,
                                             transa2, n, k, &alpha, matA, lda,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrk<double>(cublasHandle_t cublas_handle, char uplo,
                                  char transa, int32_t n, int32_t k, double alpha,
                                  double *matA, int32_t lda, double beta,
                                  double *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasDsyrk(cublas_handle, uplo_cublas,
                                             transa2, n, k, &alpha, matA, lda,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrk<cuFloatComplex>(cublasHandle_t cublas_handle,
                                          char uplo, char transa, int32_t n, int32_t k,
                                          cuFloatComplex alpha,
                                          cuFloatComplex *matA, int32_t lda,
                                          cuFloatComplex beta,
                                          cuFloatComplex *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasCsyrk(cublas_handle, uplo_cublas,
                                             transa2, n, k, &alpha, matA, lda,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrk<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char uplo, char transa, int32_t n, int32_t k,
                                           cuDoubleComplex alpha,
                                           cuDoubleComplex *matA, int32_t lda,
                                           cuDoubleComplex beta,
                                           cuDoubleComplex *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasZsyrk(cublas_handle, uplo_cublas,
                                             transa2, n, k, &alpha, matA, lda,
                                             &beta, matC, ldc));
}

/** These templates are used to select the proper
 * syrkx executable from T_data*/
template <class T_data>
cublasStatus_t carma_syrkx(cublasHandle_t cublas_handle, char uplo, char transa,
                           int32_t n, int32_t k, T_data alpha, T_data *matA, int32_t lda,
                           T_data *matB, int32_t ldb, T_data beta, T_data *matC,
                           int32_t ldc) CARMA_NIY;
/**< Specialized template for carma_syrkx executable selection */
template <>
cublasStatus_t carma_syrkx<float>(cublasHandle_t cublas_handle, char uplo,
                                  char transa, int32_t n, int32_t k, float alpha,
                                  float *matA, int32_t lda, float *matB, int32_t ldb,
                                  float beta, float *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasSsyrkx(cublas_handle, uplo_cublas,
                                              transa2, n, k, &alpha, matA, lda,
                                              matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrkx<double>(cublasHandle_t cublas_handle, char uplo,
                                   char transa, int32_t n, int32_t k, double alpha,
                                   double *matA, int32_t lda, double *matB, int32_t ldb,
                                   double beta, double *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasDsyrkx(cublas_handle, uplo_cublas,
                                              transa2, n, k, &alpha, matA, lda,
                                              matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrkx<cuFloatComplex>(
    cublasHandle_t cublas_handle, char uplo, char transa, int32_t n, int32_t k,
    cuFloatComplex alpha, cuFloatComplex *matA, int32_t lda, cuFloatComplex *matB,
    int32_t ldb, cuFloatComplex beta, cuFloatComplex *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasCsyrkx(cublas_handle, uplo_cublas,
                                              transa2, n, k, &alpha, matA, lda,
                                              matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrkx<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                            char uplo, char transa, int32_t n,
                                            int32_t k, cuDoubleComplex alpha,
                                            cuDoubleComplex *matA, int32_t lda,
                                            cuDoubleComplex *matB, int32_t ldb,
                                            cuDoubleComplex beta,
                                            cuDoubleComplex *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasZsyrkx(cublas_handle, uplo_cublas,
                                              transa2, n, k, &alpha, matA, lda,
                                              matB, ldb, &beta, matC, ldc));
}

/** These templates are used to select the proper geam
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_geam(cublasHandle_t cublas_handle, char transa,
                          char transb, int32_t m, int32_t n, T_data alpha, T_data *matA,
                          int32_t lda, T_data beta, T_data *matB, int32_t ldb,
                          T_data *matC, int32_t ldc) CARMA_NIY;
/**< Specialized template for carma_geam executable selection */
template <>
cublasStatus_t carma_geam<float>(cublasHandle_t cublas_handle, char transa,
                                 char transb, int32_t m, int32_t n, float alpha,
                                 float *matA, int32_t lda, float beta, float *matB,
                                 int32_t ldb, float *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasSgeam(cublas_handle, transa2, transb2, m,
                                             n, &alpha, matA, lda, &beta, matB,
                                             ldb, matC, ldc));
}
template <>
cublasStatus_t carma_geam<double>(cublasHandle_t cublas_handle, char transa,
                                  char transb, int32_t m, int32_t n, double alpha,
                                  double *matA, int32_t lda, double beta,
                                  double *matB, int32_t ldb, double *matC,
                                  int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasDgeam(cublas_handle, transa2, transb2, m,
                                             n, &alpha, matA, lda, &beta, matB,
                                             ldb, matC, ldc));
}
template <>
cublasStatus_t carma_geam<cuFloatComplex>(
    cublasHandle_t cublas_handle, char transa, char transb, int32_t m, int32_t n,
    cuFloatComplex alpha, cuFloatComplex *matA, int32_t lda, cuFloatComplex beta,
    cuFloatComplex *matB, int32_t ldb, cuFloatComplex *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasCgeam(cublas_handle, transa2, transb2, m,
                                             n, &alpha, matA, lda, &beta, matB,
                                             ldb, matC, ldc));
}
template <>
cublasStatus_t carma_geam<cuDoubleComplex>(
    cublasHandle_t cublas_handle, char transa, char transb, int32_t m, int32_t n,
    cuDoubleComplex alpha, cuDoubleComplex *matA, int32_t lda, cuDoubleComplex beta,
    cuDoubleComplex *matB, int32_t ldb, cuDoubleComplex *matC, int32_t ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasZgeam(cublas_handle, transa2, transb2, m,
                                             n, &alpha, matA, lda, &beta, matB,
                                             ldb, matC, ldc));
}

/** These templates are used to select the proper dgmm
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_dgmm(cublasHandle_t cublas_handle, char side, int32_t m, int32_t n,
                          const T_data *matA, int32_t lda, const T_data *vectx,
                          int32_t incx, T_data *matC, int32_t ldc) CARMA_NIY;
/**< Specialized template for carma_dgmm executable selection */
template <>
cublasStatus_t carma_dgmm<int32_t>(cublasHandle_t cublas_handle, char side, int32_t m,
                               int32_t n, const int32_t *matA, int32_t lda,
                               const int32_t *vectx, int32_t incx, int32_t *matC, int32_t ldc) {
  DEBUG_TRACE("Not implemented");
  return carma_checkCublasStatus(CUBLAS_STATUS_NOT_SUPPORTED);
}

#ifdef CAN_DO_HALF
template <>
cublasStatus_t carma_dgmm<half>(cublasHandle_t cublas_handle, char side, int32_t m,
                                int32_t n, const half *matA, int32_t lda,
                                const half *vectx, int32_t incx, half *matC,
                                int32_t ldc) {
  DEBUG_TRACE("Not implemented");
  return carma_checkCublasStatus(CUBLAS_STATUS_NOT_SUPPORTED);
}
#endif

template <>
cublasStatus_t carma_dgmm<float>(cublasHandle_t cublas_handle, char side, int32_t m,
                                 int32_t n, const float *matA, int32_t lda,
                                 const float *vectx, int32_t incx, float *matC,
                                 int32_t ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  return carma_checkCublasStatus(cublasSdgmm(
      cublas_handle, side_cublas, m, n, matA, lda, vectx, incx, matC, ldc));
}
template <>
cublasStatus_t carma_dgmm<double>(cublasHandle_t cublas_handle, char side,
                                  int32_t m, int32_t n, const double *matA, int32_t lda,
                                  const double *vectx, int32_t incx, double *matC,
                                  int32_t ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  return carma_checkCublasStatus(cublasDdgmm(
      cublas_handle, side_cublas, m, n, matA, lda, vectx, incx, matC, ldc));
}
template <>
cublasStatus_t carma_dgmm<cuFloatComplex>(cublasHandle_t cublas_handle,
                                          char side, int32_t m, int32_t n,
                                          const cuFloatComplex *matA, int32_t lda,
                                          const cuFloatComplex *vectx, int32_t incx,
                                          cuFloatComplex *matC, int32_t ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  return carma_checkCublasStatus(cublasCdgmm(
      cublas_handle, side_cublas, m, n, matA, lda, vectx, incx, matC, ldc));
}
template <>
cublasStatus_t carma_dgmm<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char side, int32_t m, int32_t n,
                                           const cuDoubleComplex *matA, int32_t lda,
                                           const cuDoubleComplex *vectx,
                                           int32_t incx, cuDoubleComplex *matC,
                                           int32_t ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  return carma_checkCublasStatus(cublasZdgmm(
      cublas_handle, side_cublas, m, n, matA, lda, vectx, incx, matC, ldc));
}

/** These templates are used to select the proper sbmv
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_sbmv(cublasHandle_t cublas_handle, char uplo, int32_t n, int32_t k,
                          T_data alpha, T_data *matA, int32_t lda, T_data *vectx,
                          int32_t incx, T_data beta, T_data *vecty,
                          int32_t incy) CARMA_NIY;
/**< Specialized template for carma_sbmv executable selection */
template<>
cublasStatus_t carma_sbmv<float>(cublasHandle_t cublas_handle, char uplo, int32_t n, int32_t k,
                                 float alpha, float *matA, int32_t lda,
                                 float *vectx, int32_t incx, float beta,
                                 float *vecty, int32_t incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);

  return carma_checkCublasStatus(cublasSsbmv(cublas_handle, filla, n, k, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}
template<>
cublasStatus_t carma_sbmv<double>(cublasHandle_t cublas_handle, char uplo, int32_t n, int32_t k,
                                  double alpha, double *matA, int32_t lda,
                                  double *vectx, int32_t incx, double beta,
                                  double *vecty, int32_t incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasDsbmv(cublas_handle, filla, n, k, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}


struct CarmaCublasInterfacer {
  template <typename T_data>
  static void call() {
    force_keep(&carma_where_amax<T_data>);
    force_keep(&carma_where_amin<T_data>);
    force_keep(&carma_getasum<T_data>);
    force_keep(&carma_axpy<T_data>);
    force_keep(&carma_copy<T_data>);
    force_keep(&carma_dot<T_data>);
    force_keep(&carma_nrm2<T_data>);
    force_keep(&carma_rot<T_data>);
    force_keep(&carma_scal<T_data>);
    force_keep(&carma_swap<T_data>);
    force_keep(&carma_gemv<T_data>);
    force_keep(&carma_ger<T_data>);
    force_keep(&carma_symv<T_data>);
    force_keep(&carma_gemm<T_data>);
    force_keep(&carma_symm<T_data>);
    force_keep(&carma_syrk<T_data>);
    force_keep(&carma_syrkx<T_data>);
    force_keep(&carma_geam<T_data>);
    force_keep(&carma_dgmm<T_data>);
    force_keep(&carma_sbmv<T_data>);
  }
};

void declare_carma_cublas() { apply<CarmaCublasInterfacer, TypeListObj>(); }
