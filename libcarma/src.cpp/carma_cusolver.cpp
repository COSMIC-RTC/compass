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

//! \file      carma_cusolver.cpp
//! \ingroup   libcarma
//! \brief     this file provides wrappers to the cuSolver functions
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2019/03/28
//! \copyright GNU Lesser General Public License

#include <carma_cusolver.h>

#include <type_list.hpp>

#ifdef CAN_DO_HALF
using TypeListObj = GenericTypeList<uint16_t, int, unsigned int, float, double,
                                    half, cuFloatComplex,
                                    cuDoubleComplex>;  // , tuple_t<float>>;
#else
using TypeListObj =
    GenericTypeList<uint16_t, int, unsigned int, float, double, cuFloatComplex,
                    cuDoubleComplex>;  // , tuple_t<float>>;
#endif

#define CHECK_CUSOLVER(fct)                                             \
  do {                                                                  \
    cusolverStatus_t info = fct;                                        \
    if (info != CUSOLVER_STATUS_SUCCESS) {                              \
      printf("%s@%d cuSolver returned error %d.\n", __FILE__, __LINE__, \
             (int)info);                                                \
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
int carma_syevd_gen(Fn_bufferSize const &ptr_syevd_gpu_bufferSize,
                    Fn const &ptr_syevd_gpu, char job, int N, int lda, T *d_mat,
                    T *d_eigenvals, CarmaDevice *device) {
  int *devInfo = NULL;
  T *d_work = NULL;
  int lwork = 0;

  int info_gpu = 0;

  // step 2: query working space of syevd
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  cusolverEigMode_t jobz =
      (job == SOLVER_EIG_MODE_VECTOR ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR);
  CHECK_CUSOLVER(ptr_syevd_gpu_bufferSize(device->get_cusolver_handle(), jobz, uplo, N, d_mat,
                                          lda, d_eigenvals, &lwork));
  cudaMalloc((void **)&d_work, sizeof(T) * lwork);
  cudaMalloc((void **)&devInfo, sizeof(int));

  // step 3: compute spectrum
  CHECK_CUSOLVER(ptr_syevd_gpu(device->get_cusolver_handle(), jobz, uplo, N, d_mat, lda,
                               d_eigenvals, d_work, lwork, devInfo));

  if (devInfo) cudaFree(devInfo);
  if (d_work) cudaFree(d_work);
  cudaDeviceSynchronize();
  return EXIT_SUCCESS;
}

template <class T>
int carma_syevd(char jobz, long N, T *mat, T *eigenvals, CarmaDevice *device) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_syevd<float>(char jobz, long N, float *mat, float *eigenvals,
                       CarmaDevice *device) {
  return carma_syevd_gen(cusolverDnSsyevd_bufferSize, cusolverDnSsyevd, jobz, N,
                         N, mat, eigenvals, device);
}

template <>
int carma_syevd<double>(char jobz, long N, double *mat, double *eigenvals,
                        CarmaDevice *device) {
  return carma_syevd_gen(cusolverDnDsyevd_bufferSize, cusolverDnDsyevd, jobz, N,
                         N, mat, eigenvals, device);
}

template <class T>
int carma_syevd(char jobz, CarmaObj<T> *mat, CarmaObj<T> *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_syevd<float>(char jobz, CarmaObjS *mat, CarmaObjS *eigenvals) {
  long N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  int device_id = mat->get_device();
  CarmaDevice *device = mat->get_context()->get_device(device_id);

  return carma_syevd_gen(cusolverDnSsyevd_bufferSize, cusolverDnSsyevd, jobz, N,
                         N, mat->get_data(), eigenvals->get_data(), device);
}

template <>
int carma_syevd<double>(char jobz, CarmaObjD *mat, CarmaObjD *eigenvals) {
  long N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  int device_id = mat->get_device();
  CarmaDevice *device = mat->get_context()->get_device(device_id);
  return carma_syevd_gen(cusolverDnDsyevd_bufferSize, cusolverDnDsyevd, jobz, N,
                         N, mat->get_data(), eigenvals->get_data(), device);
}

template <class T, typename Fn_bufferSize, typename Fn>
int carma_potr_inv_gen(Fn_bufferSize const &ptr_potri_gpu_bufferSize,
                       Fn const &ptr_potrf_gpu, Fn const &ptr_potri_gpu, int N,
                       int lda, T *d_mat, CarmaDevice *device) {
  int *devInfo = NULL;
  T *d_work = NULL;
  int lwork = 0;

  int info_gpu = 0;

  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  CHECK_CUSOLVER(ptr_potri_gpu_bufferSize(device->get_cusolver_handle(), uplo,
                                          N, d_mat, lda, &lwork));
  cudaMalloc((void **)&d_work, sizeof(T) * lwork);
  cudaMalloc((void **)&devInfo, sizeof(int));

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
int carma_potr_inv(long N, T *mat, CarmaDevice *device) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_potr_inv<float>(long N, float *mat, CarmaDevice *device) {
  return carma_potr_inv_gen(cusolverDnSpotri_bufferSize, cusolverDnSpotrf,
                            cusolverDnSpotri, N, N, mat, device);
}

template <>
int carma_potr_inv<double>(long N, double *mat, CarmaDevice *device) {
  return carma_potr_inv_gen(cusolverDnDpotri_bufferSize, cusolverDnDpotrf,
                            cusolverDnDpotri, N, N, mat, device);
}

template <class T>
int carma_potr_inv(CarmaObj<T> *mat) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_potr_inv<float>(CarmaObjS *mat) {
  long N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  int device_id = mat->get_device();
  CarmaDevice *device = mat->get_context()->get_device(device_id);
  return carma_potr_inv_gen(cusolverDnSpotri_bufferSize, cusolverDnSpotrf,
                            cusolverDnSpotri, N, N, mat->get_data(), device);
}

template <>
int carma_potr_inv<double>(CarmaObjD *mat) {
  long N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  int device_id = mat->get_device();
  CarmaDevice *device = mat->get_context()->get_device(device_id);
  return carma_potr_inv_gen(cusolverDnDpotri_bufferSize, cusolverDnDpotrf,
                            cusolverDnDpotri, N, N, mat->get_data(), device);
}

struct CarmaCusolverInterfacer {
  template <typename T_data>
  static void call() {
    force_keep((int (*)(char, long, T_data *, T_data *, CarmaDevice *)) &
               carma_syevd<T_data>);
    force_keep((int (*)(char, CarmaObj<T_data> *, CarmaObj<T_data> *)) &
               carma_syevd<T_data>);
    force_keep((int (*)(long, T_data *, CarmaDevice *)) &
               carma_potr_inv<T_data>);
    force_keep((int (*)(CarmaObj<T_data> *)) & carma_potr_inv<T_data>);
  }
};

void declare_carma() { apply<CarmaCusolverInterfacer, TypeListObj>(); }
