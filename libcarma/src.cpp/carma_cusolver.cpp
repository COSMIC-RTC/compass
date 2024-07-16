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

//! \file      carma_cusolver.cpp
//! \ingroup   libcarma
//! \brief     this file provides wrappers to the cuSolver functions
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <carma_cusolver.hpp>

#include <type_list.hpp>

using TypeListObj =
    GenericTypeList<uint16_t, int32_t, uint32_t, float, double, cuFloatComplex,
                    cuDoubleComplex>;  // , tuple_t<float>>;

#define CHECK_CUSOLVER(fct)                                             \
  do {                                                                  \
    cusolverStatus_t info = fct;                                        \
    if (info != CUSOLVER_STATUS_SUCCESS) {                              \
      printf("%s@%d cuSolver returned error %d.\n", __FILE__, __LINE__, \
             (int32_t)info);                                                \
    }                                                                   \
  } while (0)

cusolverStatus_t carma_init_cusolver(cusolverDnHandle_t *cusolver_handle) {
  cusolverStatus_t status;
  CHECK_CUSOLVER(status = cusolverDnCreate(cusolver_handle));
  return status;
}
cusolverStatus_t carma_shutdown_cusolver(cusolverDnHandle_t cusolver_handle) {
  cusolverStatus_t status;
  CHECK_CUSOLVER(status = cusolverDnDestroy(cusolver_handle));
  return status;
}

template <class T, typename Fn_bufferSize, typename Fn>
int32_t carma_syevd_gen(Fn_bufferSize const &ptr_syevd_gpu_bufferSize,
                    Fn const &ptr_syevd_gpu, char job, int32_t N, int32_t lda, T *d_mat,
                    T *d_eigenvals, CarmaDevice *device) {
  int32_t *devInfo = NULL;
  T *d_work = NULL;
  int32_t lwork = 0;

  int32_t info_gpu = 0;

  // step 2: query working space of syevd
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  cusolverEigMode_t jobz =
      (job == SOLVER_EIG_MODE_VECTOR ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR);
  CHECK_CUSOLVER(ptr_syevd_gpu_bufferSize(device->get_cusolver_handle(), jobz, uplo, N, d_mat,
                                          lda, d_eigenvals, &lwork));
  cudaMalloc((void **)&d_work, sizeof(T) * lwork);
  cudaMalloc((void **)&devInfo, sizeof(int32_t));

  // step 3: compute spectrum
  CHECK_CUSOLVER(ptr_syevd_gpu(device->get_cusolver_handle(), jobz, uplo, N, d_mat, lda,
                               d_eigenvals, d_work, lwork, devInfo));

  if (devInfo) cudaFree(devInfo);
  if (d_work) cudaFree(d_work);
  cudaDeviceSynchronize();
  return EXIT_SUCCESS;
}

template <class T>
int32_t carma_syevd(char jobz, int64_t N, T *mat, T *eigenvals, CarmaDevice *device) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int32_t carma_syevd<float>(char jobz, int64_t N, float *mat, float *eigenvals,
                       CarmaDevice *device) {
  return carma_syevd_gen(cusolverDnSsyevd_bufferSize, cusolverDnSsyevd, jobz, N,
                         N, mat, eigenvals, device);
}

template <>
int32_t carma_syevd<double>(char jobz, int64_t N, double *mat, double *eigenvals,
                        CarmaDevice *device) {
  return carma_syevd_gen(cusolverDnDsyevd_bufferSize, cusolverDnDsyevd, jobz, N,
                         N, mat, eigenvals, device);
}

template <class T>
int32_t carma_syevd(char jobz, CarmaObj<T> *mat, CarmaObj<T> *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int32_t carma_syevd<float>(char jobz, CarmaObjS *mat, CarmaObjS *eigenvals) {
  int64_t N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  int32_t device_id = mat->get_device();
  CarmaDevice *device = mat->get_context()->get_device(device_id);

  return carma_syevd_gen(cusolverDnSsyevd_bufferSize, cusolverDnSsyevd, jobz, N,
                         N, mat->get_data(), eigenvals->get_data(), device);
}

template <>
int32_t carma_syevd<double>(char jobz, CarmaObjD *mat, CarmaObjD *eigenvals) {
  int64_t N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  int32_t device_id = mat->get_device();
  CarmaDevice *device = mat->get_context()->get_device(device_id);
  return carma_syevd_gen(cusolverDnDsyevd_bufferSize, cusolverDnDsyevd, jobz, N,
                         N, mat->get_data(), eigenvals->get_data(), device);
}

template <class T, typename Fn_bufferSize, typename Fn>
int32_t carma_potr_inv_gen(Fn_bufferSize const &ptr_potri_gpu_bufferSize,
                       Fn const &ptr_potrf_gpu, Fn const &ptr_potri_gpu, int32_t N,
                       int32_t lda, T *d_mat, CarmaDevice *device) {
  int32_t *devInfo = NULL;
  T *d_work = NULL;
  int32_t lwork = 0;

  int32_t info_gpu = 0;

  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  CHECK_CUSOLVER(ptr_potri_gpu_bufferSize(device->get_cusolver_handle(), uplo,
                                          N, d_mat, lda, &lwork));
  cudaMalloc((void **)&d_work, sizeof(T) * lwork);
  cudaMalloc((void **)&devInfo, sizeof(int32_t));

  CHECK_CUSOLVER(ptr_potrf_gpu(device->get_cusolver_handle(), uplo, N, d_mat,
                               lda, d_work, lwork, devInfo));

  CHECK_CUSOLVER(ptr_potri_gpu(device->get_cusolver_handle(), uplo, N, d_mat,
                               lda, d_work, lwork, devInfo));

  fill_sym_matrix<T>('L', d_mat, N, N * lda, device);

  if (devInfo) cudaFree(devInfo);
  if (d_work) cudaFree(d_work);
  cudaDeviceSynchronize();
  return EXIT_SUCCESS;
}

template <class T>
int32_t carma_potr_inv(int64_t N, T *mat, CarmaDevice *device) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int32_t carma_potr_inv<float>(int64_t N, float *mat, CarmaDevice *device) {
  return carma_potr_inv_gen(cusolverDnSpotri_bufferSize, cusolverDnSpotrf,
                            cusolverDnSpotri, N, N, mat, device);
}

template <>
int32_t carma_potr_inv<double>(int64_t N, double *mat, CarmaDevice *device) {
  return carma_potr_inv_gen(cusolverDnDpotri_bufferSize, cusolverDnDpotrf,
                            cusolverDnDpotri, N, N, mat, device);
}

template <class T>
int32_t carma_potr_inv(CarmaObj<T> *mat) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int32_t carma_potr_inv<float>(CarmaObjS *mat) {
  int64_t N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  int32_t device_id = mat->get_device();
  CarmaDevice *device = mat->get_context()->get_device(device_id);
  return carma_potr_inv_gen(cusolverDnSpotri_bufferSize, cusolverDnSpotrf,
                            cusolverDnSpotri, N, N, mat->get_data(), device);
}

template <>
int32_t carma_potr_inv<double>(CarmaObjD *mat) {
  int64_t N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  int32_t device_id = mat->get_device();
  CarmaDevice *device = mat->get_context()->get_device(device_id);
  return carma_potr_inv_gen(cusolverDnDpotri_bufferSize, cusolverDnDpotrf,
                            cusolverDnDpotri, N, N, mat->get_data(), device);
}

struct CarmaCusolverInterfacer {
  template <typename T_data>
  static void call() {
    force_keep((int32_t (*)(char, int64_t, T_data *, T_data *, CarmaDevice *)) &
               carma_syevd<T_data>);
    force_keep((int32_t (*)(char, CarmaObj<T_data> *, CarmaObj<T_data> *)) &
               carma_syevd<T_data>);
    force_keep((int32_t (*)(int64_t, T_data *, CarmaDevice *)) &
               carma_potr_inv<T_data>);
    force_keep((int32_t (*)(CarmaObj<T_data> *)) & carma_potr_inv<T_data>);
  }
};

void declare_carma() { apply<CarmaCusolverInterfacer, TypeListObj>(); }
