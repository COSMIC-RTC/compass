// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_coronagraph.cu
//! \ingroup   libsutra
//! \class     SutraCoronagraph
//! \brief     this class provides the coronagraph features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_coronagraph.hpp>
#include <carma_utils.cuh>

__global__ void compute_electric_field_krnl(cuFloatComplex *ef, float* opd, float scale,
                            float* amplitude, float* mask, int32_t N) {
    int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {
        float A = amplitude[tid];
        float _opd = opd[tid];
        float _mask = mask[tid];
        ef[tid].x = A * cosf(scale * _opd) * _mask;
        ef[tid].y = A * sinf(scale * _opd) * _mask;
        tid += blockDim.x * gridDim.x;
    }
}

int32_t compute_electric_field(cuFloatComplex *electric_field, float* phase_opd, float scale,
                            float* amplitude, float* mask, int32_t dimx, int32_t dimy, CarmaDevice *device) {
    int32_t nBlocks, nThreads;
    get_num_blocks_and_threads(device, dimx*dimy, nBlocks, nThreads);
    dim3 grid(nBlocks), threads(nThreads);

    compute_electric_field_krnl<<<grid, threads>>>(electric_field, phase_opd, scale,
                                                    amplitude, mask, dimx*dimy);
    return EXIT_SUCCESS;
}

__global__ void remove_complex_avg_krnl(cuFloatComplex *ef, cuFloatComplex sum,
                                    float* mask, int32_t Nvalid, int32_t N) {
    int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {
        float _mask = mask[tid];
        ef[tid].x -= (sum.x / Nvalid * _mask);
        ef[tid].y -= (sum.y / Nvalid * _mask);
        tid += blockDim.x * gridDim.x;
    }
}
int32_t remove_complex_avg(cuFloatComplex *electric_field, cuFloatComplex sum, float* mask, int32_t Nvalid,
                        int32_t dimx, int32_t dimy, CarmaDevice *device) {

    int32_t nBlocks, nThreads;
    get_num_blocks_and_threads(device, dimx*dimy, nBlocks, nThreads);
    dim3 grid(nBlocks), threads(nThreads);

    remove_complex_avg_krnl<<<grid, threads>>>(electric_field, sum, mask, Nvalid, dimx*dimy);
    return EXIT_SUCCESS;
}

__global__ void accumulate_abs2_krnl(cuFloatComplex *img, float* abs2img, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    cuFloatComplex cache = img[tid];
    abs2img[tid] += (cache.x * cache.x + cache.y * cache.y);
    tid += blockDim.x * gridDim.x;
  }
}

int32_t accumulate_abs2(cuFloatComplex *img, float* abs2img, int32_t N, CarmaDevice *device) {
  int32_t nBlocks, nThreads;
  get_num_blocks_and_threads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  accumulate_abs2_krnl<<<grid, threads>>>(img, abs2img, N);

  return EXIT_SUCCESS;
}

__global__ void apply_mask_krnl(cuFloatComplex *electric_field, float* mask, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    float cache = mask[tid];
    electric_field[tid].x *= cache;
    electric_field[tid].y *= cache;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t apply_mask(cuFloatComplex *electric_field, float* mask, int32_t N, CarmaDevice *device) {
  int32_t nBlocks, nThreads;
  get_num_blocks_and_threads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  apply_mask_krnl<<<grid, threads>>>(electric_field, mask, N);

  return EXIT_SUCCESS;

}
