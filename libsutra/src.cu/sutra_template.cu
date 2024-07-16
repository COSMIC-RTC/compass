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

//! \file      sutra_template.cu
//! \ingroup   libsutra
//! \class     SutraTemplate
//! \brief     this class provides a class template to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_template.hpp>
#include "carma_utils.cuh"

template <class T>
__global__ void comp_aotemplate_krnl(T *g_idata, T *g_odata, int32_t sh_size,
                                     int32_t N) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    // fill shared mem with data
    sdata[tid] = g_idata[i];
  }

  __syncthreads();

  if (i < N) {
    // write result for this block to global mem
    g_odata[i] =
        sin((sdata[tid] - sdata[(tid + 1) % sh_size]) * 2.0f * CARMA_PI);
  }
}

template <class T>
void comp_aotemplate(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  comp_aotemplate_krnl<T>
      <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, smemSize, N);

  carma_check_msg("comp_aotemplate_kernel<<<>>> execution failed\n");
}

template void comp_aotemplate<float>(int32_t threads, int32_t blocks, float *d_idata,
                                     float *d_odata, int32_t N);

template void comp_aotemplate<double>(int32_t threads, int32_t blocks, double *d_idata,
                                      double *d_odata, int32_t N);
