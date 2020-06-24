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

//! \file      carma_svd.cu
//! \ingroup   libcarma
//! \brief     this file provides SVD CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#if 0
#include <carma_svd.h>

__global__ void kernel_setidd(double *d,int N) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) {
    d[idx+N*idx] = 1.0;
  }
}

int carma_setidd(double *d,int n) {
  int blockSize = 8;
  int nBlocks = n / blockSize + (n % blockSize == 0?0:1);

  kernel_setidd <<< nBlocks, blockSize >>> ((double *)d,n);

  return EXIT_SUCCESS;
}

__global__ void kernel_setidf(float *d,int N) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) {
    d[idx+N*idx] = 1.0;
  }
}

int carma_setidf(float *d,int n) {
  int blockSize = 8;
  int nBlocks = n / blockSize + (n % blockSize == 0?0:1);

  kernel_setidf <<< nBlocks, blockSize >>> ((float *)d,n);

  return EXIT_SUCCESS;
}
#endif
