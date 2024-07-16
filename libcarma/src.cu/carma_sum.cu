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

//! \file      carma_sum.cu
//! \ingroup   libcarma
//! \brief     this file provides summation CUDA kernels
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24


#include <carma_obj.hpp>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/complex.h>
#include "carma_utils.cuh"

/*
 Parallel sum reduction using shared memory
 - takes log(n) steps for n input elements
 - uses n threads
 - only works for power-of-2 arrays
 */

/*
 This version adds multiple elements per thread sequentially.  This reduces the
 overall cost of the algorithm while keeping the work complexity O(n) and the
 step complexity O(log n). (Brent's Theorem optimization)

 Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
 In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
 If blockSize > 32, allocate blockSize*sizeof(T) bytes.
 */
// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template <class T, uint32_t blockSize, bool nIsPow2>
__global__ void reduce6(T *g_idata, T *g_odata, uint32_t n) {
  T *sdata = SharedMemory<T>();

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockSize * 2 + threadIdx.x;
  uint32_t gridSize = blockSize * 2 * gridDim.x;

  T mySum = 0;

  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  while (i < n) {
    mySum += g_idata[i];
    // ensure we don't read out of bounds -- this is optimized away for powerOf2
    // sized arrays
    if (nIsPow2 || i + blockSize < n) mySum += g_idata[i + blockSize];
    i += gridSize;
  }

  // each thread puts its local sum into shared memory
  sdata[tid] = mySum;
  __syncthreads();

  // do reduction in shared mem
  if (blockSize >= 1024) {
    if (tid < 512) {
      sdata[tid] = mySum = mySum + sdata[tid + 512];
    }
    __syncthreads();
  }
  if (blockSize >= 512) {
    if (tid < 256) {
      sdata[tid] = mySum = mySum + sdata[tid + 256];
    }
    __syncthreads();
  }
  if (blockSize >= 256) {
    if (tid < 128) {
      sdata[tid] = mySum = mySum + sdata[tid + 128];
    }
    __syncthreads();
  }
  if (blockSize >= 128) {
    if (tid < 64) {
      sdata[tid] = mySum = mySum + sdata[tid + 64];
    }
    __syncthreads();
  }

#ifndef __DEVICE_EMULATION__
  if (tid < 32)
#endif
  {
    // now that we are using warp-synchronous programming (below)
    // we need to declare our shared memory volatile so that the compiler
    // doesn't reorder stores to it and induce incorrect behavior.
    volatile T *smem = sdata;
    if (blockSize >= 64) {
      smem[tid] = mySum = mySum + smem[tid + 32];
      __syncthreads();
    }
    if (blockSize >= 32) {
      smem[tid] = mySum = mySum + smem[tid + 16];
      __syncthreads();
    }
    if (blockSize >= 16) {
      smem[tid] = mySum = mySum + smem[tid + 8];
      __syncthreads();
    }
    if (blockSize >= 8) {
      smem[tid] = mySum = mySum + smem[tid + 4];
      __syncthreads();
    }
    if (blockSize >= 4) {
      smem[tid] = mySum = mySum + smem[tid + 2];
      __syncthreads();
    }
    if (blockSize >= 2) {
      smem[tid] = mySum = mySum + smem[tid + 1];
      __syncthreads();
    }
  }

  // write result for this block to global mem
  if (tid == 0) atomicAdd(g_odata, sdata[0]);
}

template <class T>
void reduce(int32_t size, int32_t threads, int32_t blocks, T *d_idata, T *d_odata) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

  if (is_pow2(size)) {
    switch (threads) {
      case 1024:
        reduce6<T, 1024, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 512:
        reduce6<T, 512, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 256:
        reduce6<T, 256, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 128:
        reduce6<T, 128, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 64:
        reduce6<T, 64, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 32:
        reduce6<T, 32, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 16:
        reduce6<T, 16, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 8:
        reduce6<T, 8, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 4:
        reduce6<T, 4, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 2:
        reduce6<T, 2, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 1:
        reduce6<T, 1, true>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    }
  } else {
    switch (threads) {
      case 1024:
        reduce6<T, 1024, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 512:
        reduce6<T, 512, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 256:
        reduce6<T, 256, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 128:
        reduce6<T, 128, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 64:
        reduce6<T, 64, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 32:
        reduce6<T, 32, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 16:
        reduce6<T, 16, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 8:
        reduce6<T, 8, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 4:
        reduce6<T, 4, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 2:
        reduce6<T, 2, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 1:
        reduce6<T, 1, false>
            <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    }
  }
}

template void reduce<int32_t>(int32_t size, int32_t threads, int32_t blocks, int32_t *d_idata,
                          int32_t *d_odata);

template void reduce<float>(int32_t size, int32_t threads, int32_t blocks, float *d_idata,
                            float *d_odata);

template void reduce<uint32_t>(int32_t size, int32_t threads, int32_t blocks,
                                   uint32_t *d_idata,
                                   uint32_t *d_odata);

#if (__CUDA_ARCH__ < 600)
template <>
void reduce<double>(int32_t size, int32_t threads, int32_t blocks, double *d_idata,
                    double *d_odata) {
  DEBUG_TRACE(
      "Not implemented, only supported by devices of compute capability 6.x "
      "and higher.");
}
#else
template void reduce<double>(int32_t size, int32_t threads, int32_t blocks, double *d_idata,
                             double *d_odata);
#endif

template <>
void reduce<cuFloatComplex>(int32_t size, int32_t threads, int32_t blocks,
                            cuFloatComplex *d_idata, cuFloatComplex *d_odata) {
  DEBUG_TRACE("Not implemented");
}
// template <>
// void reduce<tuple_t<float>>(int32_t size, int32_t threads, int32_t blocks,
//                             tuple_t<float> *d_idata, tuple_t<float> *d_odata)
//                             {
//   DEBUG_TRACE("Not implemented");
// }
template <>
void reduce<cuDoubleComplex>(int32_t size, int32_t threads, int32_t blocks,
                             cuDoubleComplex *d_idata,
                             cuDoubleComplex *d_odata) {
  DEBUG_TRACE("Not implemented");
}

template <class T>
T reduce(T *data, int32_t N) {
  thrust::device_ptr<T> dev_ptr(data);
  return thrust::reduce(dev_ptr, dev_ptr + N);
}

template float reduce<float>(float *data, int32_t N);

template double reduce<double>(double *data, int32_t N);

template int32_t reduce<int32_t>(int32_t *data, int32_t N);

template <>
uint32_t reduce<uint32_t>(uint32_t *data, int32_t N) {
  DEBUG_TRACE("Not implemented for this data type");
  return 0;
}

template <>
uint16_t reduce<uint16_t>(uint16_t *data, int32_t N) {
  DEBUG_TRACE("Not implemented for this data type");
  return 0;
}
template <>
cuFloatComplex reduce<cuFloatComplex>(cuFloatComplex *data, int32_t N) {
  thrust::device_ptr<thrust::complex<float>> dev_ptr(reinterpret_cast<thrust::complex<float>*>(data));
  thrust::complex<float> result = thrust::reduce(dev_ptr, dev_ptr + N);
  cuFloatComplex* sum = reinterpret_cast<cuFloatComplex*>(&result);
  return *sum;
}

template <>
cuDoubleComplex reduce<cuDoubleComplex>(cuDoubleComplex *data, int32_t N) {
  DEBUG_TRACE("Not implemented for this data type");
  return make_cuDoubleComplex(0, 0);
}

// template <>
// tuple_t<float> reduce<tuple_t<float>>(tuple_t<float> *data, int32_t N) {
//   DEBUG_TRACE("Not implemented for this data type");
//   return {0, 0.f};
// }

template <class T>
void init_reduceCubCU(T *&cub_data, size_t &cub_data_size, T *data, T *&o_data,
                      int32_t N) {
  // Determine temporary device storage requirements
  cudaMalloc(&o_data, sizeof(T));
  cub_data = NULL;
  cub_data_size = 0;
  cub::DeviceReduce::Sum(cub_data, cub_data_size, data, o_data, N);
  // Allocate temporary storage
  cudaMalloc(&cub_data, cub_data_size);
}

template void init_reduceCubCU<int32_t>(int32_t *&cub_data, size_t &cub_data_size,
                                    int32_t *data, int32_t *&o_data, int32_t N);
template void init_reduceCubCU<uint16_t>(uint16_t *&cub_data,
                                         size_t &cub_data_size, uint16_t *data,
                                         uint16_t *&o_data, int32_t N);
template void init_reduceCubCU<uint32_t>(uint32_t *&cub_data,
                                             size_t &cub_data_size,
                                             uint32_t *data,
                                             uint32_t *&o_data, int32_t N);
template void init_reduceCubCU<float>(float *&cub_data, size_t &cub_data_size,
                                      float *data, float *&o_data, int32_t N);
template void init_reduceCubCU<double>(double *&cub_data, size_t &cub_data_size,
                                       double *data, double *&o_data, int32_t N);
template <>
void init_reduceCubCU<cuFloatComplex>(cuFloatComplex *&cub_data,
                                      size_t &cub_data_size,
                                      cuFloatComplex *data,
                                      cuFloatComplex *&o_data, int32_t N) {
  DEBUG_TRACE("Not implemented");
}
// template <>
// void init_reduceCubCU<tuple_t<float>>(tuple_t<float> *&cub_data,
//                                       size_t &cub_data_size,
//                                       tuple_t<float> *data,
//                                       tuple_t<float> *&o_data, int32_t N) {
//   DEBUG_TRACE("Not implemented");
// }
template <>
void init_reduceCubCU<cuDoubleComplex>(cuDoubleComplex *&cub_data,
                                       size_t &cub_data_size,
                                       cuDoubleComplex *data,
                                       cuDoubleComplex *&o_data, int32_t N) {
  DEBUG_TRACE("Not implemented");
}

template <class T>
void reduceCubCU(T *cub_data, size_t cub_data_size, T *data, T *o_data, int32_t N, cudaStream_t stream) {
  cub::DeviceReduce::Sum(cub_data, cub_data_size, data, o_data, N, stream);
}

template void reduceCubCU<int32_t>(int32_t *cub_data, size_t cub_data_size, int32_t *data,
                               int32_t *o_data, int32_t N, cudaStream_t stream);
template void reduceCubCU<uint32_t>(uint32_t *cub_data,
                                        size_t cub_data_size,
                                        uint32_t *data,
                                        uint32_t *o_data, int32_t N, cudaStream_t stream);
template void reduceCubCU<uint16_t>(uint16_t *cub_data, size_t cub_data_size,
                                    uint16_t *data, uint16_t *o_data, int32_t , cudaStream_t stream);
template void reduceCubCU<float>(float *cub_data, size_t cub_data_size,
                                 float *data, float *o_data, int32_t N, cudaStream_t stream);
template void reduceCubCU<double>(double *cub_data, size_t cub_data_size,
                                  double *data, double *o_data, int32_t N, cudaStream_t stream);
template <>
void reduceCubCU<cuFloatComplex>(cuFloatComplex *cub_data, size_t cub_data_size,
                                 cuFloatComplex *data, cuFloatComplex *o_data,
                                 int32_t N, cudaStream_t stream) {
  DEBUG_TRACE("Not implemented");
}
// template <>
// void reduceCubCU<tuple_t<float>>(tuple_t<float> *cub_data, size_t
// cub_data_size,
//                                  tuple_t<float> *data, tuple_t<float>
//                                  *o_data, int32_t N) {
//   DEBUG_TRACE("Not implemented");
// }
template <>
void reduceCubCU<cuDoubleComplex>(cuDoubleComplex *cub_data,
                                  size_t cub_data_size, cuDoubleComplex *data,
                                  cuDoubleComplex *o_data, int32_t N, cudaStream_t stream) {
  DEBUG_TRACE("Not implemented");
}
