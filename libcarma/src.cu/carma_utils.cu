// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_utils.cu
//! \ingroup   libcarma
//! \brief     this file provides utilities CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24

#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include "carma_utils.cuh"
#include "carma_utils.h"

template <class T_data>
__global__ void find_nnz_krnl(T_data *d_data, int *colind, int *d_nnz, int N) {
  int *sdata = SharedMemory<int>();
  int tid = threadIdx.x + blockDim.x * blockIdx.x;
  int sid = threadIdx.x;
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
int find_nnz(T_data *d_data, int *colind, int N, int *d_nnz, int &h_nnz,
             CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  int smemSize = nb_threads * sizeof(int);

  find_nnz_krnl<<<grid, threads, smemSize>>>(d_data, colind, d_nnz, N);
  carma_check_msg("find_nnz_krnl<<<>>> execution failed\n");

  // wrap raw pointer with a device_ptr
  thrust::device_ptr<int> dev_ptr(colind);

  thrust::sort(dev_ptr, dev_ptr + N);
  carma_safe_call(cudaMemcpy(&h_nnz, d_nnz, sizeof(int), cudaMemcpyDeviceToHost));

  return EXIT_SUCCESS;
}

template int find_nnz<float>(float *d_data, int *colind, int N, int *d_nnz,
                             int &h_nnz, CarmaDevice *device);
template int find_nnz<double>(double *d_data, int *colind, int N, int *d_nnz,
                              int &h_nnz, CarmaDevice *device);

template <class T_data>
__global__ void fill_sparse_vect_krnl(T_data *dense_data, int *colind_sorted,
                                      T_data *values, int *colind, int *rowind,
                                      int nnz) {
  int tid = threadIdx.x + blockDim.x * blockIdx.x;
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
int fill_sparse_vect(T_data *dense_data, int *colind_sorted, T_data *values,
                     int *colind, int *rowind, int nnz, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, nnz, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_sparse_vect_krnl<<<grid, threads>>>(dense_data, colind_sorted, values,
                                           colind, rowind, nnz);
  carma_check_msg("fill_sparse_vect_krnl<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}
template int fill_sparse_vect<float>(float *dense_data, int *colind_sorted,
                                     float *values, int *colind, int *rowind,
                                     int nnz, CarmaDevice *device);
template int fill_sparse_vect<double>(double *dense_data, int *colind_sorted,
                                      double *values, int *colind, int *rowind,
                                      int nnz, CarmaDevice *device);

__global__ void float_to_double_krnl(float *i_data, double *o_data, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    o_data[tid] = (double)i_data[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void double_to_float_krnl(double *i_data, float *o_data, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    o_data[tid] = (float)i_data[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int float_to_double(float *i_data, double *o_data, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  float_to_double_krnl<<<grid, threads>>>(i_data, o_data, N);
  carma_check_msg("float_to_double_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int double_to_float(double *i_data, float *o_data, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  double_to_float_krnl<<<grid, threads>>>(i_data, o_data, N);
  carma_check_msg("float_to_double_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
__global__ void float_to_half_array_krnl(float *source, half *dest, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    dest[tid] = __float2half(source[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

half *float_to_half_array(float *source, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  half *dest;
  carma_safe_call(cudaMalloc((void **)&(dest), sizeof(half) * N));
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  float_to_half_array_krnl<<<grid, threads>>>(source, dest, N);
  carma_check_msg("float_to_half_array_krnl\n");

  return dest;
}

__global__ void half_to_float_array_krnl(half *source, float *dest, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    dest[tid] = __half2float(source[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

float *half_to_float_array(half *source, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  float *dest;
  carma_safe_call(cudaMalloc((void **)&(dest), sizeof(float) * N));
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  half_to_float_array_krnl<<<grid, threads>>>(source, dest, N);
  carma_check_msg("half_to_float_array_krnl\n");

  return dest;
}

__global__ void copy_from_float_to_half_krnl(const float *data, half *dest, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    dest[tid] = __float2half(data[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

int copy_from_float_to_half(const float *h_data, half *d_dest, int N,
                        CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  float *d_data;
  carma_safe_call(cudaMalloc((void **)&d_data, sizeof(float) * N));
  carma_safe_call(
      cudaMemcpy(d_data, h_data, sizeof(float) * N, cudaMemcpyHostToDevice));

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  copy_from_float_to_half_krnl<<<grid, threads>>>(d_data, d_dest, N);
  carma_check_msg("copy_from_float_to_half_krnl\n");
  cudaFree(d_data);

  return EXIT_SUCCESS;
}

__global__ void copy_from_half_to_float_krnl(const half *data, float *dest, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    dest[tid] = __half2float(data[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

int copy_from_half_to_float(const half *d_data, float *h_dest, int N,
                        CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  float *d_dest;
  carma_safe_call(cudaMalloc((void **)&d_dest, sizeof(float) * N));
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  copy_from_half_to_float_krnl<<<grid, threads>>>(d_data, d_dest, N);
  carma_check_msg("copy_from_half_to_float_krnl\n");
  carma_safe_call(
      cudaMemcpy(h_dest, d_dest, sizeof(float) * N, cudaMemcpyDeviceToHost));

  cudaFree(d_dest);

  return EXIT_SUCCESS;
}
#endif

template <typename T>
__global__ void fill_array_krnl(T *d_data, T value, int N) {
  int tid = threadIdx.x + blockDim.x * blockIdx.x;
  while (tid < N) {
    d_data[tid] = value;
    tid += blockDim.x * gridDim.x;
  }
}

template <typename T_data>
int fill_array_with_value(T_data *d_data, T_data value, int N,
                          CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  fill_array_krnl<<<grid, threads>>>(d_data, value, N);
  carma_check_msg("fill_array_with_value\n");

  return EXIT_SUCCESS;
}

template int fill_array_with_value<float>(float *d_data, float value, int N,
                                          CarmaDevice *device);
template int fill_array_with_value<double>(double *d_data, double value, int N,
                                           CarmaDevice *device);
template int fill_array_with_value<int>(int *d_data, int value, int N,
                                        CarmaDevice *device);
template int fill_array_with_value<unsigned int>(unsigned int *d_data,
                                                 unsigned int value, int N,
                                                 CarmaDevice *device);
template int fill_array_with_value<uint16_t>(uint16_t *d_data, uint16_t value,
                                             int N, CarmaDevice *device);
template int fill_array_with_value<cuFloatComplex>(cuFloatComplex *d_data,
                                                   cuFloatComplex value, int N,
                                                   CarmaDevice *device);
template int fill_array_with_value<cuDoubleComplex>(cuDoubleComplex *d_data,
                                                    cuDoubleComplex value,
                                                    int N,
                                                    CarmaDevice *device);

#ifdef CAN_DO_HALF
template int fill_array_with_value<half>(half *d_data, half value, int N,
                                         CarmaDevice *device);
#endif
