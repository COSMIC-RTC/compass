// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_cublas.cpp
//! \ingroup   libcarma
//! \brief     this file provides the cublas features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24

#include <carma_cublas.h>
#include <carma_utils.h>
#include <iostream>
#include <stdexcept>
#include <string>

#include <type_list.hpp>

#ifdef CAN_DO_HALF
using TypeListObj = GenericTypeList<int, unsigned int, uint16_t, float, double,
                                    half, cuFloatComplex,
                                    cuDoubleComplex>;  // , tuple_t<float>>;
#else
using TypeListObj =
    GenericTypeList<int, unsigned int, uint16_t, float, double, cuFloatComplex,
                    cuDoubleComplex>;  // , tuple_t<float>>;
#endif

#define CARMA_NIY                                            \
  {                                                          \
    DEBUG_TRACE("Method not implemented yet!");              \
    throw std::runtime_error("Method not implemented yet!"); \
  }

cublasStatus_t __carma_checkCublasStatus(cublasStatus_t status, int line,
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
          cublasStatus_t (*afunc)(cublasHandle_t handle, int n, const T_data *x,
                                  int incx, int *result)>
int carma_where_gen(cublasHandle_t cublas_handle, int n, const T_data *vect,
                    int incx) {
  int result = 0;
  carma_checkCublasStatus(afunc(cublas_handle, n, vect, incx, &result));
  return result - 1;
}

template <class T>
int carma_where_amax(cublasHandle_t cublas_handle, int n, const T *vect,
                     int incx) CARMA_NIY;
/**< Specialized template for carma_amax executable selection */
template <>
int carma_where_amax(cublasHandle_t cublas_handle, int n, const float *vect,
                     int incx) {
  return carma_where_gen<float, cublasIsamax>(cublas_handle, n, vect, incx);
}
template <>
int carma_where_amax(cublasHandle_t cublas_handle, int n, const double *vect,
                     int incx) {
  return carma_where_gen<double, cublasIdamax>(cublas_handle, n, vect, incx);
}
template <>
int carma_where_amax(cublasHandle_t cublas_handle, int n,
                     const cuFloatComplex *vect, int incx) {
  return carma_where_gen<cuFloatComplex, cublasIcamax>(cublas_handle, n, vect,
                                                       incx);
}
template <>
int carma_where_amax(cublasHandle_t cublas_handle, int n,
                     const cuDoubleComplex *vect, int incx) {
  return carma_where_gen<cuDoubleComplex, cublasIzamax>(cublas_handle, n, vect,
                                                        incx);
}

template <class T>
int carma_where_amin(cublasHandle_t cublas_handle, int n, const T *vect,
                     int incx) CARMA_NIY;
/**< Specialized template for carma_amin executable selection */
template <>
int carma_where_amin(cublasHandle_t cublas_handle, int n, const float *vect,
                     int incx) {
  return carma_where_gen<float, cublasIsamin>(cublas_handle, n, vect, incx);
}
template <>
int carma_where_amin(cublasHandle_t cublas_handle, int n, const double *vect,
                     int incx) {
  return carma_where_gen<double, cublasIdamin>(cublas_handle, n, vect, incx);
}
template <>
int carma_where_amin(cublasHandle_t cublas_handle, int n,
                     const cuFloatComplex *vect, int incx) {
  return carma_where_gen<cuFloatComplex, cublasIcamin>(cublas_handle, n, vect,
                                                       incx);
}
template <>
int carma_where_amin(cublasHandle_t cublas_handle, int n,
                     const cuDoubleComplex *vect, int incx) {
  return carma_where_gen<cuDoubleComplex, cublasIzamin>(cublas_handle, n, vect,
                                                        incx);
}

/** These templates are used to select the proper asum executable from T_data*/
template <class T_data,
          cublasStatus_t (*afunc)(cublasHandle_t handle, int n, const T_data *x,
                                  int incx, T_data *result)>
T_data carma_asum_gen(cublasHandle_t cublas_handle, int n, const T_data *vect,
                      int incx) {
  T_data result = 0;
  carma_checkCublasStatus(afunc(cublas_handle, n, vect, incx, &result));
  return result;
}

template <class T>
T carma_getasum(cublasHandle_t cublas_handle, int n, const T *vect,
                int incx) CARMA_NIY;
/**< Specialized template for carma_getasum executable selection */
template <>
float carma_getasum(cublasHandle_t cublas_handle, int n, const float *vect,
                    int incx) {
  return carma_asum_gen<float, cublasSasum>(cublas_handle, n, vect, incx);
}
template <>
double carma_getasum(cublasHandle_t cublas_handle, int n, const double *vect,
                     int incx) {
  return carma_asum_gen<double, cublasDasum>(cublas_handle, n, vect, incx);
}

/** These templates are used to select the proper axpy executable from T_data*/
template <class T_data>
cublasStatus_t carma_axpy(cublasHandle_t cublas_handle, int n,
                          const T_data alpha, const T_data *vectx, int incx,
                          T_data *vecty, int incy) CARMA_NIY;
/**< Specialized template for carma_axpy executable selection */
template <>
cublasStatus_t carma_axpy<float>(cublasHandle_t cublas_handle, int n,
                                 const float alpha, const float *vectx,
                                 int incx, float *vecty, int incy) {
  return carma_checkCublasStatus(
      cublasSaxpy(cublas_handle, n, &alpha, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_axpy<double>(cublasHandle_t cublas_handle, int n,
                                  const double alpha, const double *vectx,
                                  int incx, double *vecty, int incy) {
  return carma_checkCublasStatus(
      cublasDaxpy(cublas_handle, n, &alpha, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_axpy<cuFloatComplex>(cublasHandle_t cublas_handle, int n,
                                          const cuFloatComplex alpha,
                                          const cuFloatComplex *vectx, int incx,
                                          cuFloatComplex *vecty, int incy) {
  return carma_checkCublasStatus(
      cublasCaxpy(cublas_handle, n, &alpha, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_axpy<cuDoubleComplex>(cublasHandle_t cublas_handle, int n,
                                           const cuDoubleComplex alpha,
                                           const cuDoubleComplex *vectx,
                                           int incx, cuDoubleComplex *vecty,
                                           int incy) {
  return carma_checkCublasStatus(
      cublasZaxpy(cublas_handle, n, &alpha, vectx, incx, vecty, incy));
}

/** These templates are used to select the proper copy
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_copy(cublasHandle_t cublas_handle, int n,
                          const T_data *vectx, int incx, T_data *vecty,
                          int incy) CARMA_NIY;
/**< Specialized template for carma_copy executable selection */
template <>
cublasStatus_t carma_copy<float>(cublasHandle_t cublas_handle, int n,
                                 const float *vectx, int incx, float *vecty,
                                 int incy) {
  return carma_checkCublasStatus(
      cublasScopy(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_copy<double>(cublasHandle_t cublas_handle, int n,
                                  const double *vectx, int incx, double *vecty,
                                  int incy) {
  return carma_checkCublasStatus(
      cublasDcopy(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_copy<cuFloatComplex>(cublasHandle_t cublas_handle, int n,
                                          const cuFloatComplex *vectx, int incx,
                                          cuFloatComplex *vecty, int incy) {
  return carma_checkCublasStatus(
      cublasCcopy(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_copy<cuDoubleComplex>(cublasHandle_t cublas_handle, int n,
                                           const cuDoubleComplex *vectx,
                                           int incx, cuDoubleComplex *vecty,
                                           int incy) {
  return carma_checkCublasStatus(
      cublasZcopy(cublas_handle, n, vectx, incx, vecty, incy));
}

/** These templates are used to select the proper dot
 * executable from T_data*/
template <class T_data>
T_data carma_dot(cublasHandle_t cublas_handle, int n, T_data *vectx, int incx,
                 T_data *vecty, int incy) CARMA_NIY;
/**< Specialized template for carma_dot executable selection */
template <>
float carma_dot<float>(cublasHandle_t cublas_handle, int n, float *vectx,
                       int incx, float *vecty, int incy) {
  float result = 0;
  carma_checkCublasStatus(
      cublasSdot(cublas_handle, n, vectx, incx, vecty, incy, &result));
  return result;
}
template <>
double carma_dot<double>(cublasHandle_t cublas_handle, int n, double *vectx,
                         int incx, double *vecty, int incy) {
  double result = 0;
  carma_checkCublasStatus(
      cublasDdot(cublas_handle, n, vectx, incx, vecty, incy, &result));
  return result;
}
template <>
cuFloatComplex carma_dot<cuFloatComplex>(cublasHandle_t cublas_handle, int n,
                                         cuFloatComplex *vectx, int incx,
                                         cuFloatComplex *vecty, int incy) {
  cuFloatComplex result;
  result.x = 0;
  result.y = 0;
  carma_checkCublasStatus(
      cublasCdotu(cublas_handle, n, vectx, incx, vecty, incy, &result));
  return result;
}
template <>
cuDoubleComplex carma_dot<cuDoubleComplex>(cublasHandle_t cublas_handle, int n,
                                           cuDoubleComplex *vectx, int incx,
                                           cuDoubleComplex *vecty, int incy) {
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
T_data carma_nrm2(cublasHandle_t cublas_handle, int n, T_data *vect,
                  int incx) CARMA_NIY;
/**< Specialized template for carma_nrm2 executable selection */
template <>
float carma_nrm2<float>(cublasHandle_t cublas_handle, int n, float *vect,
                        int incx) {
  float result = 0;
  carma_checkCublasStatus(cublasSnrm2(cublas_handle, n, vect, incx, &result));
  return result;
}
template <>
double carma_nrm2<double>(cublasHandle_t cublas_handle, int n, double *vect,
                          int incx) {
  double result = 0;
  carma_checkCublasStatus(cublasDnrm2(cublas_handle, n, vect, incx, &result));
  return result;
}
/*add
 template<>
 float carma_nrm2<cuFloatComplex>(cublasHandle_t
 cublas_handle, int n, cuFloatComplex *vect, int incx)
 { float result = 0;
 carma_checkCublasStatus(cublasScnrm2(cublas_handle,n,
 vect,incx,&result)); return result;
 }*/

/** These templates are used to select the proper rot
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_rot(cublasHandle_t cublas_handle, int n, T_data *vectx,
                         int incx, T_data *vecty, int incy, T_data sc,
                         T_data ss) CARMA_NIY;
/**< Specialized template for carma_rot executable selection */
template <>
cublasStatus_t carma_rot<float>(cublasHandle_t cublas_handle, int n,
                                float *vectx, int incx, float *vecty, int incy,
                                float sc, float ss) {
  return carma_checkCublasStatus(
      cublasSrot(cublas_handle, n, vectx, incx, vecty, incy, &sc, &ss));
}
template <>
cublasStatus_t carma_rot<double>(cublasHandle_t cublas_handle, int n,
                                 double *vectx, int incx, double *vecty,
                                 int incy, double sc, double ss) {
  return carma_checkCublasStatus(
      cublasDrot(cublas_handle, n, vectx, incx, vecty, incy, &sc, &ss));
}

/** These templates are used to select the proper scal
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_scal(cublasHandle_t cublas_handle, int n, T_data alpha,
                          T_data *vectx, int incx) CARMA_NIY;
/**< Specialized template for carma_scal executable selection */
template <>
cublasStatus_t carma_scal<float>(cublasHandle_t cublas_handle, int n,
                                 float alpha, float *vectx, int incx) {
  return carma_checkCublasStatus(
      cublasSscal(cublas_handle, n, &alpha, vectx, incx));
}
template <>
cublasStatus_t carma_scal<double>(cublasHandle_t cublas_handle, int n,
                                  double alpha, double *vectx, int incx) {
  return carma_checkCublasStatus(
      cublasDscal(cublas_handle, n, &alpha, vectx, incx));
}
template <>
cublasStatus_t carma_scal<cuFloatComplex>(cublasHandle_t cublas_handle, int n,
                                          cuFloatComplex alpha,
                                          cuFloatComplex *vectx, int incx) {
  return carma_checkCublasStatus(
      cublasCscal(cublas_handle, n, &alpha, vectx, incx));
}
template <>
cublasStatus_t carma_scal<cuDoubleComplex>(cublasHandle_t cublas_handle, int n,
                                           cuDoubleComplex alpha,
                                           cuDoubleComplex *vectx, int incx) {
  return carma_checkCublasStatus(
      cublasZscal(cublas_handle, n, &alpha, vectx, incx));
}

/** These templates are used to select the proper swap
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_swap(cublasHandle_t cublas_handle, int n, T_data *vectx,
                          int incx, T_data *vecty, int incy) CARMA_NIY;
/**< Specialized template for carma_swap executable selection */
template <>
cublasStatus_t carma_swap<float>(cublasHandle_t cublas_handle, int n,
                                 float *vectx, int incx, float *vecty,
                                 int incy) {
  return carma_checkCublasStatus(
      cublasSswap(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_swap<double>(cublasHandle_t cublas_handle, int n,
                                  double *vectx, int incx, double *vecty,
                                  int incy) {
  return carma_checkCublasStatus(
      cublasDswap(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_swap<cuFloatComplex>(cublasHandle_t cublas_handle, int n,
                                          cuFloatComplex *vectx, int incx,
                                          cuFloatComplex *vecty, int incy) {
  return carma_checkCublasStatus(
      cublasCswap(cublas_handle, n, vectx, incx, vecty, incy));
}
template <>
cublasStatus_t carma_swap<cuDoubleComplex>(cublasHandle_t cublas_handle, int n,
                                           cuDoubleComplex *vectx, int incx,
                                           cuDoubleComplex *vecty, int incy) {
  return carma_checkCublasStatus(
      cublasZswap(cublas_handle, n, vectx, incx, vecty, incy));
}

/** These templates are used to select the proper gemv
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_gemv(cublasHandle_t cublas_handle, char trans, int m,
                          int n, T_data alpha, T_data *matA, int lda,
                          T_data *vectx, int incx, T_data beta, T_data *vecty,
                          int incy) CARMA_NIY;
/**< Specialized template for carma_gemv executable selection */
template <>
cublasStatus_t carma_gemv<float>(cublasHandle_t cublas_handle, char trans,
                                 int m, int n, float alpha, float *matA,
                                 int lda, float *vectx, int incx, float beta,
                                 float *vecty, int incy) {
  cublasOperation_t trans2 = carma_char2cublas_operation(trans);
  return carma_checkCublasStatus(cublasSgemv(cublas_handle, trans2, m, n,
                                             &alpha, matA, lda, vectx, incx,
                                             &beta, vecty, incy));
}
template <>
cublasStatus_t carma_gemv<double>(cublasHandle_t cublas_handle, char trans,
                                  int m, int n, double alpha, double *matA,
                                  int lda, double *vectx, int incx, double beta,
                                  double *vecty, int incy) {
  cublasOperation_t trans2 = carma_char2cublas_operation(trans);
  return carma_checkCublasStatus(cublasDgemv(cublas_handle, trans2, m, n,
                                             &alpha, matA, lda, vectx, incx,
                                             &beta, vecty, incy));
}
template <>
cublasStatus_t carma_gemv<cuFloatComplex>(
    cublasHandle_t cublas_handle, char trans, int m, int n,
    cuFloatComplex alpha, cuFloatComplex *matA, int lda, cuFloatComplex *vectx,
    int incx, cuFloatComplex beta, cuFloatComplex *vecty, int incy) {
  cublasOperation_t trans2 = carma_char2cublas_operation(trans);
  return carma_checkCublasStatus(cublasCgemv(cublas_handle, trans2, m, n,
                                             &alpha, matA, lda, vectx, incx,
                                             &beta, vecty, incy));
}
template <>
cublasStatus_t carma_gemv<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char trans, int m, int n,
                                           cuDoubleComplex alpha,
                                           cuDoubleComplex *matA, int lda,
                                           cuDoubleComplex *vectx, int incx,
                                           cuDoubleComplex beta,
                                           cuDoubleComplex *vecty, int incy) {
  cublasOperation_t trans2 = carma_char2cublas_operation(trans);
  return carma_checkCublasStatus(cublasZgemv(cublas_handle, trans2, m, n,
                                             &alpha, matA, lda, vectx, incx,
                                             &beta, vecty, incy));
}
#ifdef CAN_DO_HALF
template <>
cublasStatus_t carma_gemv<half>(cublasHandle_t cublas_handle, char trans, int m,
                                int n, half alpha, half *matA, int lda,
                                half *vectx, int incx, half beta, half *vecty,
                                int incy) {
  int k = (((trans == 'N') || (trans == 'n')) ? n : m);
  cublasOperation_t transa2 = carma_char2cublas_operation(trans);
  cublasOperation_t transb2 = carma_char2cublas_operation('N');
  return carma_checkCublasStatus(cublasHgemm(cublas_handle, transa2, transb2, m,
                                             1, k, &alpha, matA, lda, vectx, n,
                                             &beta, vecty, n));
}
#endif
/** These templates are used to select the proper ger executable from T_data*/
template <class T_data>
cublasStatus_t carma_ger(cublasHandle_t cublas_handle, int m, int n,
                         T_data alpha, T_data *vectx, int incx, T_data *vecty,
                         int incy, T_data *matA, int lda) CARMA_NIY;
/**< Specialized template for carma_ger executable selection */
template <>
cublasStatus_t carma_ger<float>(cublasHandle_t cublas_handle, int m, int n,
                                float alpha, float *vectx, int incx,
                                float *vecty, int incy, float *matA, int lda) {
  return carma_checkCublasStatus(cublasSger(cublas_handle, m, n, &alpha, vectx,
                                            incx, vecty, incy, matA, lda));
}
template <>
cublasStatus_t carma_ger<double>(cublasHandle_t cublas_handle, int m, int n,
                                 double alpha, double *vectx, int incx,
                                 double *vecty, int incy, double *matA,
                                 int lda) {
  return carma_checkCublasStatus(cublasDger(cublas_handle, m, n, &alpha, vectx,
                                            incx, vecty, incy, matA, lda));
}
template <>
cublasStatus_t carma_ger<cuFloatComplex>(cublasHandle_t cublas_handle, int m,
                                         int n, cuFloatComplex alpha,
                                         cuFloatComplex *vectx, int incx,
                                         cuFloatComplex *vecty, int incy,
                                         cuFloatComplex *matA, int lda) {
  return carma_checkCublasStatus(cublasCgeru(cublas_handle, m, n, &alpha, vectx,
                                             incx, vecty, incy, matA, lda));
}
template <>
cublasStatus_t carma_ger<cuDoubleComplex>(cublasHandle_t cublas_handle, int m,
                                          int n, cuDoubleComplex alpha,
                                          cuDoubleComplex *vectx, int incx,
                                          cuDoubleComplex *vecty, int incy,
                                          cuDoubleComplex *matA, int lda) {
  return carma_checkCublasStatus(cublasZgeru(cublas_handle, m, n, &alpha, vectx,
                                             incx, vecty, incy, matA, lda));
}

/** These templates are used to select the proper symv
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_symv(cublasHandle_t cublas_handle, char uplo, int n,
                          T_data alpha, T_data *matA, int lda, T_data *vectx,
                          int incx, T_data beta, T_data *vecty,
                          int incy) CARMA_NIY;
/**< Specialized template for carma_symv executable selection */
template <>
cublasStatus_t carma_symv<float>(cublasHandle_t cublas_handle, char uplo, int n,
                                 float alpha, float *matA, int lda,
                                 float *vectx, int incx, float beta,
                                 float *vecty, int incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);

  return carma_checkCublasStatus(cublasSsymv(cublas_handle, filla, n, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}
template <>
cublasStatus_t carma_symv<double>(cublasHandle_t cublas_handle, char uplo,
                                  int n, double alpha, double *matA, int lda,
                                  double *vectx, int incx, double beta,
                                  double *vecty, int incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasDsymv(cublas_handle, filla, n, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}
template <>
cublasStatus_t carma_symv<cuFloatComplex>(
    cublasHandle_t cublas_handle, char uplo, int n, cuFloatComplex alpha,
    cuFloatComplex *matA, int lda, cuFloatComplex *vectx, int incx,
    cuFloatComplex beta, cuFloatComplex *vecty, int incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasCsymv(cublas_handle, filla, n, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}
template <>
cublasStatus_t carma_symv<cuDoubleComplex>(
    cublasHandle_t cublas_handle, char uplo, int n, cuDoubleComplex alpha,
    cuDoubleComplex *matA, int lda, cuDoubleComplex *vectx, int incx,
    cuDoubleComplex beta, cuDoubleComplex *vecty, int incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasZsymv(cublas_handle, filla, n, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}

/** These templates are used to select the proper gemm
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_gemm(cublasHandle_t cublas_handle, char transa,
                          char transb, int m, int n, int k, T_data alpha,
                          T_data *matA, int lda, T_data *matB, int ldb,
                          T_data beta, T_data *matC, int ldc) CARMA_NIY;
/**< Specialized template for carma_gemm executable selection */
template <>
cublasStatus_t carma_gemm<float>(cublasHandle_t cublas_handle, char transa,
                                 char transb, int m, int n, int k, float alpha,
                                 float *matA, int lda, float *matB, int ldb,
                                 float beta, float *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasSgemm(cublas_handle, transa2, transb2, m,
                                             n, k, &alpha, matA, lda, matB, ldb,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_gemm<double>(cublasHandle_t cublas_handle, char transa,
                                  char transb, int m, int n, int k,
                                  double alpha, double *matA, int lda,
                                  double *matB, int ldb, double beta,
                                  double *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasDgemm(cublas_handle, transa2, transb2, m,
                                             n, k, &alpha, matA, lda, matB, ldb,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_gemm<cuFloatComplex>(
    cublasHandle_t cublas_handle, char transa, char transb, int m, int n, int k,
    cuFloatComplex alpha, cuFloatComplex *matA, int lda, cuFloatComplex *matB,
    int ldb, cuFloatComplex beta, cuFloatComplex *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasCgemm(cublas_handle, transa2, transb2, m,
                                             n, k, &alpha, matA, lda, matB, ldb,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_gemm<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char transa, char transb, int m,
                                           int n, int k, cuDoubleComplex alpha,
                                           cuDoubleComplex *matA, int lda,
                                           cuDoubleComplex *matB, int ldb,
                                           cuDoubleComplex beta,
                                           cuDoubleComplex *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasZgemm(cublas_handle, transa2, transb2, m,
                                             n, k, &alpha, matA, lda, matB, ldb,
                                             &beta, matC, ldc));
}

#ifdef CAN_DO_HALF
template <>
cublasStatus_t carma_gemm<half>(cublasHandle_t cublas_handle, char transa,
                                char transb, int m, int n, int k, half alpha,
                                half *matA, int lda, half *matB, int ldb,
                                half beta, half *matC, int ldc) {
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
                    char transa, char transb, int m, int n, int k, T_data alpha,
                    T_data *matsA, int lda, long long int strideA, 
                    T_data *matsB, int ldb, long long int strideB, T_data beta,
                    T_data *matsC, int ldc, long long int strideC, 
                    int batch_count) CARMA_NIY;
/**< Specialized template for carma_gemm_batched executable selection */
template <>
cublasStatus_t carma_gemm_strided_batched<float>(cublasHandle_t cublas_handle, 
                    char transa, char transb, int m, int n, int k, float alpha, 
                    float *matsA, int lda, long long int strideA,
                    float *matsB, int ldb, long long int strideB, float beta,
                    float *matsC, int ldc, long long int strideC,
                    int batch_count) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasSgemmStridedBatched(cublas_handle, 
                      transa2, transb2, m, n, k, &alpha, matsA, lda, strideA,
                      matsB, ldb, strideB, &beta, matsC, ldc, strideC, 
                      batch_count));
}
template <>
cublasStatus_t carma_gemm_strided_batched<double>(cublasHandle_t cublas_handle,
                    char transa, char transb, int m, int n, int k, double alpha, 
                    double *matsA, int lda, long long int strideA,
                    double *matsB, int ldb, long long int strideB, double beta,
                    double *matsC, int ldc, long long int strideC,
                    int batch_count) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasDgemmStridedBatched(cublas_handle, 
                      transa2, transb2, m, n, k, &alpha, matsA, lda, strideA,
                      matsB, ldb, strideB, &beta, matsC, ldc, strideC, 
                      batch_count));
}
template <>
cublasStatus_t carma_gemm_strided_batched<cuFloatComplex>(
    cublasHandle_t cublas_handle, char transa, char transb, int m, int n, int k,
    cuFloatComplex alpha, cuFloatComplex *matsA, int lda, long long int strideA,
    cuFloatComplex *matsB, int ldb, long long int strideB, cuFloatComplex beta,
    cuFloatComplex *matsC, int ldc, long long int strideC, int batch_count) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasCgemmStridedBatched(cublas_handle,
                      transa2, transb2, m, n, k, &alpha, matsA, lda, strideA,
                      matsB, ldb, strideB, &beta, matsC, ldc, strideC,
                      batch_count));
}
template <>
cublasStatus_t carma_gemm_strided_batched<cuDoubleComplex>(
  cublasHandle_t cublas_handle, char transa, char transb, int m, int n, int k,
  cuDoubleComplex alpha, cuDoubleComplex *matsA, int lda, long long int strideA,
  cuDoubleComplex *matsB, int ldb, long long int strideB, cuDoubleComplex beta,
  cuDoubleComplex *matsC, int ldc, long long int strideC, int batch_count) {
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
                      char transa, char transb, int m, int n, int k, half alpha,
                      half *matsA, int lda, long long int strideA,
                      half *matsB, int ldb, long long int strideB, half beta, 
                      half *matsC, int ldc, long long int strideC, 
                      int batch_count) {
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
                          int m, int n, T_data alpha, T_data *matA, int lda,
                          T_data *matB, int ldb, T_data beta, T_data *matC,
                          int ldc) CARMA_NIY;
/**< Specialized template for carma_symm executable selection */
template <>
cublasStatus_t carma_symm<float>(cublasHandle_t cublas_handle, char side,
                                 char uplo, int m, int n, float alpha,
                                 float *matA, int lda, float *matB, int ldb,
                                 float beta, float *matC, int ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);

  return carma_checkCublasStatus(cublasSsymm(cublas_handle, side_cublas,
                                             uplo_cublas, m, n, &alpha, matA,
                                             lda, matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_symm<double>(cublasHandle_t cublas_handle, char side,
                                  char uplo, int m, int n, double alpha,
                                  double *matA, int lda, double *matB, int ldb,
                                  double beta, double *matC, int ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasDsymm(cublas_handle, side_cublas,
                                             uplo_cublas, m, n, &alpha, matA,
                                             lda, matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_symm<cuFloatComplex>(
    cublasHandle_t cublas_handle, char side, char uplo, int m, int n,
    cuFloatComplex alpha, cuFloatComplex *matA, int lda, cuFloatComplex *matB,
    int ldb, cuFloatComplex beta, cuFloatComplex *matC, int ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasCsymm(cublas_handle, side_cublas,
                                             uplo_cublas, m, n, &alpha, matA,
                                             lda, matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_symm<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char side, char uplo, int m, int n,
                                           cuDoubleComplex alpha,
                                           cuDoubleComplex *matA, int lda,
                                           cuDoubleComplex *matB, int ldb,
                                           cuDoubleComplex beta,
                                           cuDoubleComplex *matC, int ldc) {
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
                          int n, int k, T_data alpha, T_data *matA, int lda,
                          T_data beta, T_data *matC, int ldc) CARMA_NIY;
/**< Specialized template for carma_syrk executable selection */
template <>
cublasStatus_t carma_syrk<float>(cublasHandle_t cublas_handle, char uplo,
                                 char transa, int n, int k, float alpha,
                                 float *matA, int lda, float beta, float *matC,
                                 int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasSsyrk(cublas_handle, uplo_cublas,
                                             transa2, n, k, &alpha, matA, lda,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrk<double>(cublasHandle_t cublas_handle, char uplo,
                                  char transa, int n, int k, double alpha,
                                  double *matA, int lda, double beta,
                                  double *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasDsyrk(cublas_handle, uplo_cublas,
                                             transa2, n, k, &alpha, matA, lda,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrk<cuFloatComplex>(cublasHandle_t cublas_handle,
                                          char uplo, char transa, int n, int k,
                                          cuFloatComplex alpha,
                                          cuFloatComplex *matA, int lda,
                                          cuFloatComplex beta,
                                          cuFloatComplex *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasCsyrk(cublas_handle, uplo_cublas,
                                             transa2, n, k, &alpha, matA, lda,
                                             &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrk<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char uplo, char transa, int n, int k,
                                           cuDoubleComplex alpha,
                                           cuDoubleComplex *matA, int lda,
                                           cuDoubleComplex beta,
                                           cuDoubleComplex *matC, int ldc) {
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
                           int n, int k, T_data alpha, T_data *matA, int lda,
                           T_data *matB, int ldb, T_data beta, T_data *matC,
                           int ldc) CARMA_NIY;
/**< Specialized template for carma_syrkx executable selection */
template <>
cublasStatus_t carma_syrkx<float>(cublasHandle_t cublas_handle, char uplo,
                                  char transa, int n, int k, float alpha,
                                  float *matA, int lda, float *matB, int ldb,
                                  float beta, float *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasSsyrkx(cublas_handle, uplo_cublas,
                                              transa2, n, k, &alpha, matA, lda,
                                              matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrkx<double>(cublasHandle_t cublas_handle, char uplo,
                                   char transa, int n, int k, double alpha,
                                   double *matA, int lda, double *matB, int ldb,
                                   double beta, double *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasDsyrkx(cublas_handle, uplo_cublas,
                                              transa2, n, k, &alpha, matA, lda,
                                              matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrkx<cuFloatComplex>(
    cublasHandle_t cublas_handle, char uplo, char transa, int n, int k,
    cuFloatComplex alpha, cuFloatComplex *matA, int lda, cuFloatComplex *matB,
    int ldb, cuFloatComplex beta, cuFloatComplex *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasFillMode_t uplo_cublas = carma_char2cublasFillMode(uplo);
  return carma_checkCublasStatus(cublasCsyrkx(cublas_handle, uplo_cublas,
                                              transa2, n, k, &alpha, matA, lda,
                                              matB, ldb, &beta, matC, ldc));
}
template <>
cublasStatus_t carma_syrkx<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                            char uplo, char transa, int n,
                                            int k, cuDoubleComplex alpha,
                                            cuDoubleComplex *matA, int lda,
                                            cuDoubleComplex *matB, int ldb,
                                            cuDoubleComplex beta,
                                            cuDoubleComplex *matC, int ldc) {
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
                          char transb, int m, int n, T_data alpha, T_data *matA,
                          int lda, T_data beta, T_data *matB, int ldb,
                          T_data *matC, int ldc) CARMA_NIY;
/**< Specialized template for carma_geam executable selection */
template <>
cublasStatus_t carma_geam<float>(cublasHandle_t cublas_handle, char transa,
                                 char transb, int m, int n, float alpha,
                                 float *matA, int lda, float beta, float *matB,
                                 int ldb, float *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasSgeam(cublas_handle, transa2, transb2, m,
                                             n, &alpha, matA, lda, &beta, matB,
                                             ldb, matC, ldc));
}
template <>
cublasStatus_t carma_geam<double>(cublasHandle_t cublas_handle, char transa,
                                  char transb, int m, int n, double alpha,
                                  double *matA, int lda, double beta,
                                  double *matB, int ldb, double *matC,
                                  int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasDgeam(cublas_handle, transa2, transb2, m,
                                             n, &alpha, matA, lda, &beta, matB,
                                             ldb, matC, ldc));
}
template <>
cublasStatus_t carma_geam<cuFloatComplex>(
    cublasHandle_t cublas_handle, char transa, char transb, int m, int n,
    cuFloatComplex alpha, cuFloatComplex *matA, int lda, cuFloatComplex beta,
    cuFloatComplex *matB, int ldb, cuFloatComplex *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasCgeam(cublas_handle, transa2, transb2, m,
                                             n, &alpha, matA, lda, &beta, matB,
                                             ldb, matC, ldc));
}
template <>
cublasStatus_t carma_geam<cuDoubleComplex>(
    cublasHandle_t cublas_handle, char transa, char transb, int m, int n,
    cuDoubleComplex alpha, cuDoubleComplex *matA, int lda, cuDoubleComplex beta,
    cuDoubleComplex *matB, int ldb, cuDoubleComplex *matC, int ldc) {
  cublasOperation_t transa2 = carma_char2cublas_operation(transa);
  cublasOperation_t transb2 = carma_char2cublas_operation(transb);
  return carma_checkCublasStatus(cublasZgeam(cublas_handle, transa2, transb2, m,
                                             n, &alpha, matA, lda, &beta, matB,
                                             ldb, matC, ldc));
}

/** These templates are used to select the proper dgmm
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_dgmm(cublasHandle_t cublas_handle, char side, int m, int n,
                          const T_data *matA, int lda, const T_data *vectx,
                          int incx, T_data *matC, int ldc) CARMA_NIY;
/**< Specialized template for carma_dgmm executable selection */
template <>
cublasStatus_t carma_dgmm<int>(cublasHandle_t cublas_handle, char side, int m,
                               int n, const int *matA, int lda,
                               const int *vectx, int incx, int *matC, int ldc) {
  DEBUG_TRACE("Not implemented");
  return carma_checkCublasStatus(CUBLAS_STATUS_NOT_SUPPORTED);
}

#ifdef CAN_DO_HALF
template <>
cublasStatus_t carma_dgmm<half>(cublasHandle_t cublas_handle, char side, int m,
                                int n, const half *matA, int lda,
                                const half *vectx, int incx, half *matC,
                                int ldc) {
  DEBUG_TRACE("Not implemented");
  return carma_checkCublasStatus(CUBLAS_STATUS_NOT_SUPPORTED);
}
#endif

template <>
cublasStatus_t carma_dgmm<float>(cublasHandle_t cublas_handle, char side, int m,
                                 int n, const float *matA, int lda,
                                 const float *vectx, int incx, float *matC,
                                 int ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  return carma_checkCublasStatus(cublasSdgmm(
      cublas_handle, side_cublas, m, n, matA, lda, vectx, incx, matC, ldc));
}
template <>
cublasStatus_t carma_dgmm<double>(cublasHandle_t cublas_handle, char side,
                                  int m, int n, const double *matA, int lda,
                                  const double *vectx, int incx, double *matC,
                                  int ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  return carma_checkCublasStatus(cublasDdgmm(
      cublas_handle, side_cublas, m, n, matA, lda, vectx, incx, matC, ldc));
}
template <>
cublasStatus_t carma_dgmm<cuFloatComplex>(cublasHandle_t cublas_handle,
                                          char side, int m, int n,
                                          const cuFloatComplex *matA, int lda,
                                          const cuFloatComplex *vectx, int incx,
                                          cuFloatComplex *matC, int ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  return carma_checkCublasStatus(cublasCdgmm(
      cublas_handle, side_cublas, m, n, matA, lda, vectx, incx, matC, ldc));
}
template <>
cublasStatus_t carma_dgmm<cuDoubleComplex>(cublasHandle_t cublas_handle,
                                           char side, int m, int n,
                                           const cuDoubleComplex *matA, int lda,
                                           const cuDoubleComplex *vectx,
                                           int incx, cuDoubleComplex *matC,
                                           int ldc) {
  cublasSideMode_t side_cublas = carma_char2cublasSide(side);
  return carma_checkCublasStatus(cublasZdgmm(
      cublas_handle, side_cublas, m, n, matA, lda, vectx, incx, matC, ldc));
}

/** These templates are used to select the proper sbmv
 * executable from T_data*/
template <class T_data>
cublasStatus_t carma_sbmv(cublasHandle_t cublas_handle, char uplo, int n, int k,
                          T_data alpha, T_data *matA, int lda, T_data *vectx,
                          int incx, T_data beta, T_data *vecty,
                          int incy) CARMA_NIY;
/**< Specialized template for carma_sbmv executable selection */
template<>
cublasStatus_t carma_sbmv<float>(cublasHandle_t cublas_handle, char uplo, int n, int k,
                                 float alpha, float *matA, int lda,
                                 float *vectx, int incx, float beta,
                                 float *vecty, int incy) {
  cublasFillMode_t filla = carma_char2cublasFillMode(uplo);

  return carma_checkCublasStatus(cublasSsbmv(cublas_handle, filla, n, k, &alpha,
                                             matA, lda, vectx, incx, &beta,
                                             vecty, incy));
}
template<>
cublasStatus_t carma_sbmv<double>(cublasHandle_t cublas_handle, char uplo, int n, int k,
                                  double alpha, double *matA, int lda,
                                  double *vectx, int incx, double beta,
                                  double *vecty, int incy) {
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
