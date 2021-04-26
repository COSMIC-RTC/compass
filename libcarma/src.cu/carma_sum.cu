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

//! \file      carma_sum.cu
//! \ingroup   libcarma
//! \brief     this file provides summation CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#include <carma_obj.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
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
template <class T, unsigned int blockSize, bool nIsPow2>
__global__ void reduce6(T *g_idata, T *g_odata, unsigned int n) {
  T *sdata = SharedMemory<T>();

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockSize * 2 + threadIdx.x;
  unsigned int gridSize = blockSize * 2 * gridDim.x;

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
void reduce(int size, int threads, int blocks, T *d_idata, T *d_odata) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
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

template void reduce<int>(int size, int threads, int blocks, int *d_idata,
                          int *d_odata);

template void reduce<float>(int size, int threads, int blocks, float *d_idata,
                            float *d_odata);

template void reduce<unsigned int>(int size, int threads, int blocks,
                                   unsigned int *d_idata,
                                   unsigned int *d_odata);

#if (__CUDA_ARCH__ < 600)
template <>
void reduce<double>(int size, int threads, int blocks, double *d_idata,
                    double *d_odata) {
  DEBUG_TRACE(
      "Not implemented, only supported by devices of compute capability 6.x "
      "and higher.");
}
#else
template void reduce<double>(int size, int threads, int blocks, double *d_idata,
                             double *d_odata);
#endif

template <>
void reduce<cuFloatComplex>(int size, int threads, int blocks,
                            cuFloatComplex *d_idata, cuFloatComplex *d_odata) {
  DEBUG_TRACE("Not implemented");
}
// template <>
// void reduce<tuple_t<float>>(int size, int threads, int blocks,
//                             tuple_t<float> *d_idata, tuple_t<float> *d_odata)
//                             {
//   DEBUG_TRACE("Not implemented");
// }
template <>
void reduce<cuDoubleComplex>(int size, int threads, int blocks,
                             cuDoubleComplex *d_idata,
                             cuDoubleComplex *d_odata) {
  DEBUG_TRACE("Not implemented");
}
#ifdef CAN_DO_HALF
template <>
void reduce<half>(int size, int threads, int blocks, half *d_idata,
                  half *d_odata) {
  DEBUG_TRACE("Not implemented");
}
#endif

template <class T>
T reduce(T *data, int N) {
  thrust::device_ptr<T> dev_ptr(data);
  return thrust::reduce(dev_ptr, dev_ptr + N);
}

template float reduce<float>(float *data, int N);

template double reduce<double>(double *data, int N);

template int reduce<int>(int *data, int N);

template <>
unsigned int reduce<unsigned int>(unsigned int *data, int N) {
  DEBUG_TRACE("Not implemented for this data type");
  return 0;
}

template <>
uint16_t reduce<uint16_t>(uint16_t *data, int N) {
  DEBUG_TRACE("Not implemented for this data type");
  return 0;
}
template <>
cuFloatComplex reduce<cuFloatComplex>(cuFloatComplex *data, int N) {
  DEBUG_TRACE("Not implemented for this data type");
  return make_cuComplex(0, 0);
}

template <>
cuDoubleComplex reduce<cuDoubleComplex>(cuDoubleComplex *data, int N) {
  DEBUG_TRACE("Not implemented for this data type");
  return make_cuDoubleComplex(0, 0);
}

#ifdef CAN_DO_HALF
template <>
half reduce<half>(half *data, int N) {
  DEBUG_TRACE("Not implemented for thhis data type");
  return 0;
}
#endif

// template <>
// tuple_t<float> reduce<tuple_t<float>>(tuple_t<float> *data, int N) {
//   DEBUG_TRACE("Not implemented for this data type");
//   return {0, 0.f};
// }

template <class T>
void init_reduceCubCU(T *&cub_data, size_t &cub_data_size, T *data, T *&o_data,
                      int N) {
  // Determine temporary device storage requirements
  cudaMalloc(&o_data, sizeof(T));
  cub_data = NULL;
  cub_data_size = 0;
  cub::DeviceReduce::Sum(cub_data, cub_data_size, data, o_data, N);
  // Allocate temporary storage
  cudaMalloc(&cub_data, cub_data_size);
}

template void init_reduceCubCU<int>(int *&cub_data, size_t &cub_data_size,
                                    int *data, int *&o_data, int N);
template void init_reduceCubCU<uint16_t>(uint16_t *&cub_data,
                                         size_t &cub_data_size, uint16_t *data,
                                         uint16_t *&o_data, int N);
template void init_reduceCubCU<unsigned int>(unsigned int *&cub_data,
                                             size_t &cub_data_size,
                                             unsigned int *data,
                                             unsigned int *&o_data, int N);
template void init_reduceCubCU<float>(float *&cub_data, size_t &cub_data_size,
                                      float *data, float *&o_data, int N);
template void init_reduceCubCU<double>(double *&cub_data, size_t &cub_data_size,
                                       double *data, double *&o_data, int N);
template <>
void init_reduceCubCU<cuFloatComplex>(cuFloatComplex *&cub_data,
                                      size_t &cub_data_size,
                                      cuFloatComplex *data,
                                      cuFloatComplex *&o_data, int N) {
  DEBUG_TRACE("Not implemented");
}
// template <>
// void init_reduceCubCU<tuple_t<float>>(tuple_t<float> *&cub_data,
//                                       size_t &cub_data_size,
//                                       tuple_t<float> *data,
//                                       tuple_t<float> *&o_data, int N) {
//   DEBUG_TRACE("Not implemented");
// }
template <>
void init_reduceCubCU<cuDoubleComplex>(cuDoubleComplex *&cub_data,
                                       size_t &cub_data_size,
                                       cuDoubleComplex *data,
                                       cuDoubleComplex *&o_data, int N) {
  DEBUG_TRACE("Not implemented");
}
#ifdef CAN_DO_HALF
template void init_reduceCubCU<half>(half *&cub_data, size_t &cub_data_size,
                                     half *data, half *&o_data, int N);
#endif

template <class T>
void reduceCubCU(T *cub_data, size_t cub_data_size, T *data, T *o_data, int N, cudaStream_t stream) {
  cub::DeviceReduce::Sum(cub_data, cub_data_size, data, o_data, N, stream);
}

template void reduceCubCU<int>(int *cub_data, size_t cub_data_size, int *data,
                               int *o_data, int N, cudaStream_t stream);
template void reduceCubCU<unsigned int>(unsigned int *cub_data,
                                        size_t cub_data_size,
                                        unsigned int *data,
                                        unsigned int *o_data, int N, cudaStream_t stream);
template void reduceCubCU<uint16_t>(uint16_t *cub_data, size_t cub_data_size,
                                    uint16_t *data, uint16_t *o_data, int , cudaStream_t stream);
template void reduceCubCU<float>(float *cub_data, size_t cub_data_size,
                                 float *data, float *o_data, int N, cudaStream_t stream);
template void reduceCubCU<double>(double *cub_data, size_t cub_data_size,
                                  double *data, double *o_data, int N, cudaStream_t stream);
template <>
void reduceCubCU<cuFloatComplex>(cuFloatComplex *cub_data, size_t cub_data_size,
                                 cuFloatComplex *data, cuFloatComplex *o_data,
                                 int N, cudaStream_t stream) {
  DEBUG_TRACE("Not implemented");
}
// template <>
// void reduceCubCU<tuple_t<float>>(tuple_t<float> *cub_data, size_t
// cub_data_size,
//                                  tuple_t<float> *data, tuple_t<float>
//                                  *o_data, int N) {
//   DEBUG_TRACE("Not implemented");
// }
template <>
void reduceCubCU<cuDoubleComplex>(cuDoubleComplex *cub_data,
                                  size_t cub_data_size, cuDoubleComplex *data,
                                  cuDoubleComplex *o_data, int N, cudaStream_t stream) {
  DEBUG_TRACE("Not implemented");
}
#ifdef CAN_DO_HALF
template void reduceCubCU<half>(half *cub_data, size_t cub_data_size,
                                half *data, half *o_data, int N, cudaStream_t stream);
#endif
