// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_template.cu
//! \ingroup   libsutra
//! \class     SutraTemplate
//! \brief     this class provides a class template to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#include <sutra_template.h>
#include "carma_utils.cuh"

template <class T>
__global__ void comp_aotemplate_krnl(T *g_idata, T *g_odata, int sh_size,
                                     int N) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

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
void comp_aotemplate(int threads, int blocks, T *d_idata, T *d_odata, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  comp_aotemplate_krnl<T>
      <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, smemSize, N);

  carma_check_msg("comp_aotemplate_kernel<<<>>> execution failed\n");
}

template void comp_aotemplate<float>(int threads, int blocks, float *d_idata,
                                     float *d_odata, int N);

template void comp_aotemplate<double>(int threads, int blocks, double *d_idata,
                                      double *d_odata, int N);
