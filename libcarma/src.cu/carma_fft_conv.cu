// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_fft_conv.cu
//! \ingroup   libcarma
//! \brief     this file provides fft convolution CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <carma_obj.hpp>
#include <convolutionFFT2D_common.h>

__global__ void fftconv_upadkrnl(float *odata, float *idata, int32_t fftW,
                                 int32_t dataW, int32_t N, int32_t n) {
  __shared__ float cache[BLOCK_SZ][BLOCK_SZ];

  int32_t x = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t y = threadIdx.y + blockIdx.y * blockDim.y;
  // int32_t tid = x + y *blockDim.x * gridDim.x;

  if (y * fftW + x < N)
    cache[BLOCK_SZ - 1 - threadIdx.x][BLOCK_SZ - 1 - threadIdx.y] =
        idata[y * fftW + x];

  __syncthreads();

  if (y * dataW + x < n)
    odata[y * dataW + x] =
        cache[BLOCK_SZ - 1 - threadIdx.x][BLOCK_SZ - 1 - threadIdx.y];
}

int32_t fftconv_unpad_old(float *d_odata, float *d_idata, int32_t fftW, int32_t dataH,
                      int32_t dataW, int32_t N, int32_t n) {
  dim3 blocks(dataH / BLOCK_SZ, dataW / BLOCK_SZ), threads(BLOCK_SZ, BLOCK_SZ);

  fftconv_upadkrnl<<<blocks, threads>>>(d_odata, d_idata, fftW, dataW, N, n);

  return EXIT_SUCCESS;
}

__global__ void unpad_krnl(float *odata, float *idata, int32_t fftW, int32_t dataW,
                           int32_t N, int32_t n, int32_t nim) {
  const int32_t y = blockDim.y * blockIdx.y + threadIdx.y;
  const int32_t x = blockDim.x * blockIdx.x + threadIdx.x;
  const int32_t z = blockDim.z * blockIdx.z + threadIdx.z;

  int32_t kz_src = z * N;
  int32_t kz_dst = z * n;

  if ((y * fftW + x < N) && (z < nim)) {
    odata[y * dataW + x + kz_dst] = idata[y * fftW + x + kz_src];
  }
}

int32_t fftconv_unpad(float *d_odata, float *d_idata, int32_t fftW, int32_t dataH,
                  int32_t dataW, int32_t N, int32_t n, int32_t nim) {
  dim3 threads(16, 8, 8);
  dim3 grid(i_div_up(dataW, threads.x), i_div_up(dataH, threads.y),
            i_div_up(nim, threads.z));

  unpad_krnl<<<grid, threads>>>(d_odata, d_idata, fftW, dataW, N, n, nim);
  carma_check_msg("unpad_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
