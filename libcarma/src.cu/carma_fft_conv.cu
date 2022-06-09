// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      carma_fft_conv.cu
//! \ingroup   libcarma
//! \brief     this file provides fft convolution CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <carma_obj.h>
#include <convolutionFFT2D_common.h>

__global__ void fftconv_upadkrnl(float *odata, float *idata, int fftW,
                                 int dataW, int N, int n) {
  __shared__ float cache[BLOCK_SZ][BLOCK_SZ];

  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  // int tid = x + y *blockDim.x * gridDim.x;

  if (y * fftW + x < N)
    cache[BLOCK_SZ - 1 - threadIdx.x][BLOCK_SZ - 1 - threadIdx.y] =
        idata[y * fftW + x];

  __syncthreads();

  if (y * dataW + x < n)
    odata[y * dataW + x] =
        cache[BLOCK_SZ - 1 - threadIdx.x][BLOCK_SZ - 1 - threadIdx.y];
}

int fftconv_unpad_old(float *d_odata, float *d_idata, int fftW, int dataH,
                      int dataW, int N, int n) {
  dim3 blocks(dataH / BLOCK_SZ, dataW / BLOCK_SZ), threads(BLOCK_SZ, BLOCK_SZ);

  fftconv_upadkrnl<<<blocks, threads>>>(d_odata, d_idata, fftW, dataW, N, n);

  return EXIT_SUCCESS;
}

__global__ void unpad_krnl(float *odata, float *idata, int fftW, int dataW,
                           int N, int n, int nim) {
  const int y = blockDim.y * blockIdx.y + threadIdx.y;
  const int x = blockDim.x * blockIdx.x + threadIdx.x;
  const int z = blockDim.z * blockIdx.z + threadIdx.z;

  int kz_src = z * N;
  int kz_dst = z * n;

  if ((y * fftW + x < N) && (z < nim)) {
    odata[y * dataW + x + kz_dst] = idata[y * fftW + x + kz_src];
  }
}

int fftconv_unpad(float *d_odata, float *d_idata, int fftW, int dataH,
                  int dataW, int N, int n, int nim) {
  dim3 threads(16, 8, 8);
  dim3 grid(i_div_up(dataW, threads.x), i_div_up(dataH, threads.y),
            i_div_up(nim, threads.z));

  unpad_krnl<<<grid, threads>>>(d_odata, d_idata, fftW, dataW, N, n, nim);
  carma_check_msg("unpad_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
