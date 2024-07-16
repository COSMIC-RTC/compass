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

//! \file      carma_utils.cuh
//! \ingroup   libcarma
//! \brief     this file provides tools to CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _CARMA_UTILS_CUH_
#define _CARMA_UTILS_CUH_

#include <cufft.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <driver_types.h>
#include <vector_types.h>

#include <cub/cub.cuh>
#include <cuda.h>
#include <cuda_runtime_api.h>

template <class T> struct SharedMemory {
  __device__ inline operator T *() {
    extern __shared__ int32_t __smem[];
    return (T *)__smem;
  }

  __device__ inline operator const T *() const {
    extern __shared__ int32_t __smem[];
    return (T *)__smem;
  }
};

// specialize for double to avoid unaligned memory
// access compile errors
template <> struct SharedMemory<double> {
  __device__ inline operator double *() {
    extern __shared__ double __smem_d[];
    return (double *)__smem_d;
  }

  __device__ inline operator const double *() const {
    extern __shared__ double __smem_d[];
    return (double *)__smem_d;
  }
};
template <> struct SharedMemory<float> {
  __device__ inline operator float *() {
    extern __shared__ float __smem_f[];
    return (float *)__smem_f;
  }

  __device__ inline operator const float *() const {
    extern __shared__ float __smem_f[];
    return (float *)__smem_f;
  }
};

template <class T> __device__ inline void mswap(T &a, T &b) {
  T tmp = a;
  a = b;
  b = tmp;
}

template <class T>
__inline__ __device__ void reduce_krnl(T *sdata, int32_t size, int32_t n) {
  if (size & (size - 1)) { // test if size is not a power of 2
    uint32_t s;
    if ((size & 1) != 0)
      s = size / 2 + 1; //(size&1)==size%2
    else
      s = size / 2;
    uint32_t s_old = size;
    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        sdata[n] += sdata[n + s];
      }
      //__threadfence_block();
      //__threadfence();
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2 * s < s_old) && (s != 0))
        s += 1;
    }
  } else {
    // do reduction in shared mem
    for (uint32_t s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        sdata[n] += sdata[n + s];
      }
      //__threadfence_block();
      //__threadfence();
      __syncthreads();
    }
  }
}

template <class T_data>
__inline__ __device__ T_data carma_clip(T_data n, T_data min, T_data max) {
  return n > max ? max : (n < min ? min : n);
}

#endif //_CARMA_UTILS_CUH_
