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

//! \file      sutra_gamora.cu
//! \ingroup   libsutra
//! \class     SutraGamora
//! \brief     this class provides the gamora features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_gamora.hpp>

__global__ void fillamplikrnl(cuFloatComplex *amplipup, float *phase,
                              int32_t *wherephase, float scale, int32_t Npts, int32_t nx,
                              int32_t Nx, int32_t puponly) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t nim;
  int32_t nline;
  int32_t ncol;

  while (tid < Npts) {
    nim = wherephase[tid];
    nline = nim / nx;
    ncol = nim - nline * nx;
    nim = ncol + nline * Nx;
    if (puponly == 1) {
      amplipup[nim].x = 1.0f;
      amplipup[nim].y = 0.0f;
    } else if (puponly == 0) {
      amplipup[nim].x = cosf(scale * phase[tid]);
      amplipup[nim].y = sinf(scale * phase[tid]);
    } else if (puponly == 2) {
      amplipup[nim].x = phase[tid];
      amplipup[nim].y = 0.0f;
    }
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fill_amplipup(cuFloatComplex *amplipup, float *phase, int32_t *wherephase,
                  float scale, int32_t Npts, int32_t nx, int32_t Nx, int32_t puponly,
                  CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, Npts, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fillamplikrnl<<<grid, threads>>>(amplipup, phase, wherephase, scale, Npts, nx,
                                   Nx, puponly);
  carma_check_msg("fillamplikrnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void cumulpsf_krnl(float *odata, cuFloatComplex *idata, int32_t N) {
  cuFloatComplex cache;

  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    cache = idata[tid];
    odata[tid] += (cache.x * cache.x + cache.y * cache.y);
    tid += blockDim.x * gridDim.x;
  }
}

int32_t cumulpsf(float *d_odata, cuFloatComplex *d_idata, int32_t N,
             CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  cumulpsf_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carma_check_msg("cumulpsf_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void abs2complex_krnl(cuFloatComplex *d_odata,
                                 cuFloatComplex *d_idata, int32_t N) {
  cuFloatComplex idata;
  cuFloatComplex odata;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    idata = d_idata[tid];
    odata.x = idata.x * idata.x + idata.y * idata.y;
    odata.y = 0.0f;
    d_odata[tid] = odata;
    tid += blockDim.x * gridDim.x;
  }
}
int32_t abs2complex(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t N,
                CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  abs2complex_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carma_check_msg("abs2complex_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void modulus2_krnl(float *d_odata, cuFloatComplex *d_idata, int32_t N) {
  cuFloatComplex idata;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    idata = d_idata[tid];
    d_odata[tid] = idata.x * idata.x + idata.y * idata.y;
    tid += blockDim.x * gridDim.x;
  }
}
int32_t modulus2(float *d_odata, cuFloatComplex *d_idata, int32_t N,
             CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  modulus2_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carma_check_msg("modulus2_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void real_krnl(float *d_odata, cuFloatComplex *d_idata, int32_t N) {
  cuFloatComplex cache;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    cache = d_idata[tid];
    d_odata[tid] = cache.x;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t real(float *d_odata, cuFloatComplex *d_idata, int32_t N, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  real_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carma_check_msg("real_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillmask_krnl(float *d_odata, float *d_idata, int32_t N, int32_t norm) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    d_odata[tid] = d_idata[tid] < (1e-5 * norm) ? 0.0f : 1.0f;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fill_mask(float *d_odata, float *d_idata, int32_t N, int32_t norm,
              CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fillmask_krnl<<<grid, threads>>>(d_odata, d_idata, N, norm);
  carma_check_msg("fillmask_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void pow2_krnl(cuFloatComplex *d_odata, cuFloatComplex *d_idata,
                          int32_t N) {
  cuFloatComplex idata;
  cuFloatComplex odata;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    idata = d_idata[tid];
    odata.x = idata.x * idata.x;
    odata.y = 0.0f;
    d_odata[tid] = odata;
    tid += blockDim.x * gridDim.x;
  }
}
int32_t pow2(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t N,
         CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  pow2_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carma_check_msg("pow2_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillterm1_krnl(float *d_odata, cuFloatComplex *d_idata,
                               cuFloatComplex *d_pupfft, int32_t N) {
  cuFloatComplex idata;
  cuFloatComplex pupfft;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    idata = d_idata[tid];
    pupfft = d_pupfft[tid];
    d_odata[tid] = idata.x * pupfft.x + idata.y * pupfft.y;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fill_term1(float *d_odata, cuFloatComplex *d_idata,
               cuFloatComplex *d_pupfft, int32_t N, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fillterm1_krnl<<<grid, threads>>>(d_odata, d_idata, d_pupfft, N);
  carma_check_msg("fillterm1_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void add2Dphi_krnl(cuFloatComplex *d_odata, float *d_term1,
                              float *d_term2, float e, int32_t N) {
  cuFloatComplex cache;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    cache = d_odata[tid];
    cache.x += 2 * ((d_term1[tid] - d_term2[tid]) * e);
    cache.y = 0.0f;
    d_odata[tid] = cache;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t add2Dphi(cuFloatComplex *d_odata, float *d_term1, float *d_term2, float e,
             int32_t N, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  add2Dphi_krnl<<<grid, threads>>>(d_odata, d_term1, d_term2, e, N);
  carma_check_msg("add2Dphi_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void computeOTFvii_krnl(float *d_otfVii, cuFloatComplex *d_Dphi,
                                   float *d_otftel, float *d_mask, float scale,
                                   int32_t N) {
  cuFloatComplex Dphi;
  float tmp;
  float den;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    Dphi = d_Dphi[tid];
    den = d_otftel[tid] > 1e-9 ? 1.0f / d_otftel[tid] : 0.0f;
    tmp = Dphi.x * den * d_mask[tid] * scale * scale;
    d_otfVii[tid] = expf(-0.5 * tmp) * d_mask[tid];

    tid += blockDim.x * gridDim.x;
  }
}
int32_t computeOTFvii(float *d_otfVii, cuFloatComplex *d_Dphi, float *d_otftel,
                  float *d_mask, float scale, int32_t N, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  computeOTFvii_krnl<<<grid, threads>>>(d_otfVii, d_Dphi, d_otftel, d_mask,
                                        scale, N);
  carma_check_msg("computeOTFvii_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void ifftscale_krnl(cuFloatComplex *d_odata, float scale, int32_t N) {
  cuFloatComplex cache;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    cache = d_odata[tid];
    cache.x = cache.x * scale;
    cache.y = cache.y * scale;
    d_odata[tid] = cache;
    tid += blockDim.x * gridDim.x;
  }
}
int32_t ifftscale(cuFloatComplex *d_odata, float scale, int32_t N,
              CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  ifftscale_krnl<<<grid, threads>>>(d_odata, scale, N);
  carma_check_msg("iffscale_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
