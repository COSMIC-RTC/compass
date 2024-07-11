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

//! \file      carma_svd.cu
//! \ingroup   libcarma
//! \brief     this file provides SVD CUDA kernels
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
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
