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

//! \file      carma_cusolver.cpp
//! \ingroup   libcarma
//! \brief     this file provides wrappers to the cuSolver functions
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.1
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

#define CHECK_CUSOLVER(fct, dealloc)                                    \
  do {                                                                  \
    cusolverStatus_t info = fct;                                        \
    if (info != CUSOLVER_STATUS_SUCCESS) {                              \
      printf("%s@%d cuSolver returned error %d.\n", __FILE__, __LINE__, \
             (int)info);                                                \
      dealloc;                                                          \
      return EXIT_FAILURE;                                              \
    }                                                                   \
  } while (0)

template <class T, typename Fn_bufferSize, typename Fn>
int carma_syevd_gen(Fn_bufferSize const &ptr_syevd_gpu_bufferSize,
                    Fn const &ptr_syevd_gpu, cusolverEigMode_t jobz, int N,
                    int lda, T *d_mat, T *d_eigenvals) {
  cusolverDnHandle_t cusolverH = NULL;

  int *devInfo = NULL;
  T *d_work = NULL;
  int lwork = 0;

  int info_gpu = 0;

  // step 1: create cusolver/cublas handle
  CHECK_CUSOLVER(cusolverDnCreate(&cusolverH), nullptr);

  // step 2: query working space of syevd
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  CHECK_CUSOLVER(ptr_syevd_gpu_bufferSize(cusolverH, jobz, uplo, N, d_mat, lda,
                                          d_eigenvals, &lwork),
                 nullptr);
  cudaMalloc((void **)&d_work, sizeof(T) * lwork);
  cudaMalloc((void **)&devInfo, sizeof(int));

  // step 3: compute spectrum
  CHECK_CUSOLVER(ptr_syevd_gpu(cusolverH, jobz, uplo, N, d_mat, lda,
                               d_eigenvals, d_work, lwork, devInfo),
                 nullptr);

  if (devInfo) cudaFree(devInfo);
  if (d_work) cudaFree(d_work);
  cudaDeviceSynchronize();
  return EXIT_SUCCESS;
}

template <class T>
int carma_syevd(cusolverEigMode_t jobz, long N, T *mat, T *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_syevd<float>(cusolverEigMode_t jobz, long N, float *mat,
                       float *eigenvals) {
  return carma_syevd_gen(cusolverDnSsyevd_bufferSize, cusolverDnSsyevd, jobz, N,
                         N, mat, eigenvals);
}

template <>
int carma_syevd<double>(cusolverEigMode_t jobz, long N, double *mat,
                        double *eigenvals) {
  return carma_syevd_gen(cusolverDnDsyevd_bufferSize, cusolverDnDsyevd, jobz, N,
                         N, mat, eigenvals);
}

template <class T>
int carma_syevd(cusolverEigMode_t jobz, carma_obj<T> *mat,
                carma_obj<T> *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_syevd<float>(cusolverEigMode_t jobz, caObjS *mat, caObjS *eigenvals) {
  long N = mat->getDims(1);

  if (N != mat->getDims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  return carma_syevd_gen(cusolverDnSsyevd_bufferSize, cusolverDnSsyevd, jobz, N,
                         N, mat->getData(), eigenvals->getData());
}

template <>
int carma_syevd<double>(cusolverEigMode_t jobz, caObjD *mat,
                        caObjD *eigenvals) {
  long N = mat->getDims(1);

  if (N != mat->getDims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  return carma_syevd_gen(cusolverDnDsyevd_bufferSize, cusolverDnDsyevd, jobz, N,
                         N, mat->getData(), eigenvals->getData());
}

struct CarmaCusolverInterfacer {
  template <typename T_data>
  static void call() {
    force_keep((int (*)(cusolverEigMode_t, long, T_data *, T_data *)) &
               carma_syevd<T_data>);
    force_keep(
        (int (*)(cusolverEigMode_t, carma_obj<T_data> *, carma_obj<T_data> *)) &
        carma_syevd<T_data>);
  }
};

void declare_carma() { apply<CarmaCusolverInterfacer, TypeListObj>(); }
