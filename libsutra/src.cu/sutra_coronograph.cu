// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_coronograph.cu
//! \ingroup   libsutra
//! \class     SutraCoronograph
//! \brief     this class provides the coronograph features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#include <sutra_coronograph.h>
#include <carma_utils.cuh>

__global__ void compute_electric_field_krnl(cuFloatComplex *ef, float* opd, float scale, 
                            float* amplitude, float* mask, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    cuFloatComplex cache;
    while (tid < N) {
        cache = ef[tid];
        cache.x = amplitude[tid] * cosf(scale * opd[tid]) * mask[tid];
        cache.y = amplitude[tid] * sinf(scale * opd[tid]) * mask[tid];
        tid += blockDim.x * gridDim.x;
    } 
}

int compute_electric_field(cuFloatComplex *electric_field, float* phase_opd, float scale, 
                            float* amplitude, float* mask, int dimx, int dimy, CarmaDevice *device) {
    int nBlocks, nThreads;
    get_num_blocks_and_threads(device, dimx*dimy, nBlocks, nThreads);
    dim3 grid(nb_blocks), threads(nb_threads);

    compute_electric_field_krnl<<<grid, threads>>>(electric_field, phase_opd, scale,
                                                    amplitude, mask, dimx*dimy);
    return EXIT_SUCCESS;
}

__global__ void remove_complex_avg(cuFloatComplex *ef, cuFloatComplex sum, 
                                    float* mask, int Nvalid, int N) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    cuFloatComplex cache;
    while (tid < N) {
        cache = ef[tid];
        cache.x -= (sum / Nvalid * mask[tid]);
        cache.y = (sum / Nvalid * mask[tid]);
        tid += blockDim.x * gridDim.x;
    } 
}
int remove_complex_avg(cuFloatComplex *electric_field, cuFloatComplex sum, float* mask, int Nvalid, 
                        int dimx, int dimy, CarmaDevice *device) {
    
    int nBlocks, nThreads;
    get_num_blocks_and_threads(device, dimx*dimy, nBlocks, nThreads);
    dim3 grid(nb_blocks), threads(nb_threads);

    remove_complex_avg_krnl<<<grid, threads>>>(electric_field, sum, mask, Nvalid, dimx*dimy);
    return EXIT_SUCCESS;
}

__global__ void accumulate_abs2_krnl(cuFloatComplex *img, float* abs2img, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  cuFloatComplex cache;
  while (tid < N) {
    cache = img[tid];
    abs2img[tid] += (cache.x * cache.x + cache.y * cache.y);
    tid += blockDim.x * gridDim.x;
  }
}

int accumulate_abs2(cuFloatComplex *img, float* abs2img, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  accumulate_abs2_krnl<<<grid, threads>>>(img, abs2img, N);

  return EXIT_SUCCESS;

}