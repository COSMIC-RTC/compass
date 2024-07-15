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

//! \file      carma_utils.cu
//! \ingroup   libcarma
//! \brief     this file provides utilities CUDA kernels
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include "carma_utils.cuh"
#include "carma_utils.hpp"

template <class T_data>
__global__ void find_nnz_krnl(T_data *d_data, int32_t *colind, int32_t *d_nnz, int32_t N) {
  int32_t *sdata = SharedMemory<int32_t>();
  int32_t tid = threadIdx.x + blockDim.x * blockIdx.x;
  int32_t sid = threadIdx.x;
  if (tid == 0) d_nnz[0] = 0;

  // Load shared memory with 1 if d_data[tid]!= 0, with 0 else
  if (tid < N) {
    sdata[sid] = (d_data[tid] != 0);
    colind[tid] = (sdata[sid]) ? tid : N + tid;  // Init colind for further sort
  } else {
    sdata[sid] = 0;
  }
  __syncthreads();
  reduce_krnl(sdata, blockDim.x, sid);
  __syncthreads();

  if (threadIdx.x == 0)
    //		intensities[blockIdx.x] = sdata[0];
    atomicAdd(d_nnz, sdata[0]);
}
template <class T_data>
int32_t find_nnz(T_data *d_data, int32_t *colind, int32_t N, int32_t *d_nnz, int32_t &h_nnz,
             CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  int32_t smemSize = nb_threads * sizeof(int32_t);

  find_nnz_krnl<<<grid, threads, smemSize>>>(d_data, colind, d_nnz, N);
  carma_check_msg("find_nnz_krnl<<<>>> execution failed\n");

  // wrap raw pointer with a device_ptr
  thrust::device_ptr<int32_t> dev_ptr(colind);

  thrust::sort(dev_ptr, dev_ptr + N);
  carma_safe_call(cudaMemcpy(&h_nnz, d_nnz, sizeof(int32_t), cudaMemcpyDeviceToHost));

  return EXIT_SUCCESS;
}

template int32_t find_nnz<float>(float *d_data, int32_t *colind, int32_t N, int32_t *d_nnz,
                             int32_t &h_nnz, CarmaDevice *device);
template int32_t find_nnz<double>(double *d_data, int32_t *colind, int32_t N, int32_t *d_nnz,
                              int32_t &h_nnz, CarmaDevice *device);

template <class T_data>
__global__ void fill_sparse_vect_krnl(T_data *dense_data, int32_t *colind_sorted,
                                      T_data *values, int32_t *colind, int32_t *rowind,
                                      int32_t nnz) {
  int32_t tid = threadIdx.x + blockDim.x * blockIdx.x;
  if (tid == 0) rowind[0] = 0;
  if (tid == 1) rowind[1] = nnz;

  // Load shared memory with 1 if d_data[tid]!= 0, with 0 else
  if (tid < nnz) {
    values[tid] = dense_data[colind_sorted[tid]];
    colind[tid] = colind_sorted[tid];
  }
  __syncthreads();
}
template <class T_data>
int32_t fill_sparse_vect(T_data *dense_data, int32_t *colind_sorted, T_data *values,
                     int32_t *colind, int32_t *rowind, int32_t nnz, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, nnz, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_sparse_vect_krnl<<<grid, threads>>>(dense_data, colind_sorted, values,
                                           colind, rowind, nnz);
  carma_check_msg("fill_sparse_vect_krnl<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}
template int32_t fill_sparse_vect<float>(float *dense_data, int32_t *colind_sorted,
                                     float *values, int32_t *colind, int32_t *rowind,
                                     int32_t nnz, CarmaDevice *device);
template int32_t fill_sparse_vect<double>(double *dense_data, int32_t *colind_sorted,
                                      double *values, int32_t *colind, int32_t *rowind,
                                      int32_t nnz, CarmaDevice *device);

__global__ void float_to_double_krnl(float *i_data, double *o_data, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    o_data[tid] = (double)i_data[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void double_to_float_krnl(double *i_data, float *o_data, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    o_data[tid] = (float)i_data[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int32_t float_to_double(float *i_data, double *o_data, int32_t N, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  float_to_double_krnl<<<grid, threads>>>(i_data, o_data, N);
  carma_check_msg("float_to_double_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int32_t double_to_float(double *i_data, float *o_data, int32_t N, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  double_to_float_krnl<<<grid, threads>>>(i_data, o_data, N);
  carma_check_msg("float_to_double_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template <typename T>
__global__ void fill_array_krnl(T *d_data, T value, int32_t N) {
  int32_t tid = threadIdx.x + blockDim.x * blockIdx.x;
  while (tid < N) {
    d_data[tid] = value;
    tid += blockDim.x * gridDim.x;
  }
}

template <typename T_data>
int32_t fill_array_with_value(T_data *d_data, T_data value, int32_t N,
                          CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  fill_array_krnl<<<grid, threads>>>(d_data, value, N);
  carma_check_msg("fill_array_with_value\n");

  return EXIT_SUCCESS;
}

template int32_t fill_array_with_value<float>(float *d_data, float value, int32_t N,
                                          CarmaDevice *device);
template int32_t fill_array_with_value<double>(double *d_data, double value, int32_t N,
                                           CarmaDevice *device);
template int32_t fill_array_with_value<int32_t>(int32_t *d_data, int32_t value, int32_t N,
                                        CarmaDevice *device);
template int32_t fill_array_with_value<uint32_t>(uint32_t *d_data,
                                                 uint32_t value, int32_t N,
                                                 CarmaDevice *device);
template int32_t fill_array_with_value<uint16_t>(uint16_t *d_data, uint16_t value,
                                             int32_t N, CarmaDevice *device);
template int32_t fill_array_with_value<cuFloatComplex>(cuFloatComplex *d_data,
                                                   cuFloatComplex value, int32_t N,
                                                   CarmaDevice *device);
template int32_t fill_array_with_value<cuDoubleComplex>(cuDoubleComplex *d_data,
                                                    cuDoubleComplex value,
                                                    int32_t N,
                                                    CarmaDevice *device);