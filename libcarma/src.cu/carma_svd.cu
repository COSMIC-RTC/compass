// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_svd.cu
//! \ingroup   libcarma
//! \brief     this file provides SVD CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#if 0
#include <carma_svd.hpp>

__global__ void kernel_setidd(double *d,int32_t N) {
  int32_t idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) {
    d[idx+N*idx] = 1.0;
  }
}

int32_t carma_setidd(double *d,int32_t n) {
  int32_t blockSize = 8;
  int32_t nb_blocks = n / blockSize + (n % blockSize == 0?0:1);

  kernel_setidd <<< nb_blocks, blockSize >>> ((double *)d,n);

  return EXIT_SUCCESS;
}

__global__ void kernel_setidf(float *d,int32_t N) {
  int32_t idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) {
    d[idx+N*idx] = 1.0;
  }
}

int32_t carma_setidf(float *d,int32_t n) {
  int32_t blockSize = 8;
  int32_t nb_blocks = n / blockSize + (n % blockSize == 0?0:1);

  kernel_setidf <<< nb_blocks, blockSize >>> ((float *)d,n);

  return EXIT_SUCCESS;
}
#endif
