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

//! \file      sutra_tscreen.cu
//! \ingroup   libsutra
//! \class     sutra_tscreen
//! \brief     this class provides the tscreen features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <cublas_v2.h>
#include <sutra_tscreen.h>

extern __shared__ float cache_shm[];

__global__ void vonkarman_krnl(cuFloatComplex *odata, float *idata, float k0,
                               int nalias, int nx, int ny, int blockSize) {
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  cache_shm[threadIdx.x + threadIdx.y * blockSize] = 0.0f;

  if ((x < nx) && (y < ny)) {
    // generate von karman spectrum
    for (int i = -nalias; i <= nalias; i++) {
      for (int j = -nalias; j <= nalias; j++) {
        if ((i == 0) && (j == 0)) {
          float xc = nx / 2;
          float yc = ny / 2;
          float tmp =
              sqrtf((xc - x) * (xc - x) + (yc - y) * (yc - y) + k0 * k0);
          if (tmp > 1.)
            cache_shm[threadIdx.x + threadIdx.y * blockSize] =
                (6.88f * 0.00969f) * pow(tmp, -1.83333f);
          else
            cache_shm[threadIdx.x + threadIdx.y * blockSize] =
                (6.88f * 0.00969f);
        } else {
          float xc = x * nx + nx / 2;
          float yc = y * ny + ny / 2;
          cache_shm[threadIdx.x + threadIdx.y * blockSize] +=
              (6.88f * 0.00969f) *
              pow(sqrtf((xc - x) * (xc - x) + (yc - y) * (yc - y) + k0 * k0),
                  -1.83333f);
        }
      }
    }

    odata[x + y * nx].x = cache_shm[threadIdx.x + threadIdx.y * blockSize] *
                          cosf(2.0f * CARMA_PI * idata[x + y * nx]);
    odata[x + y * nx].y = cache_shm[threadIdx.x + threadIdx.y * blockSize] *
                          sinf(2.0f * CARMA_PI * idata[x + y * nx]);
  }

  if ((x == 0) && (y == 0)) {
    odata[x + y * nx].x = 0.0f;
    odata[x + y * nx].y = 0.0f;
  }
}

int gene_vonkarman(cuFloatComplex *d_odata, float *d_idata, float k0,
                   int nalias, int nx, int ny, int block_size) {
  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  int nny = ny + block_size - ny % block_size;
  dim3 blocks(nnx / block_size, nny / block_size),
      threads(block_size, block_size);

  int smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  vonkarman_krnl<<<blocks, threads, smemSize>>>(d_odata, d_idata, k0, nalias,
                                                nx, ny, block_size);

  carmaCheckMsg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void dphix_krnl(float *odata, float *idata, int N, int iter,
                           int nx) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid + iter < N) {
    if (tid % nx < nx - iter)
      odata[tid] =
          (idata[tid] - idata[tid + iter]) * (idata[tid] - idata[tid + iter]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void dphiy_krnl(float *odata, float *idata, int N, int iter,
                           int nx) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid + iter * nx < N) {
    odata[tid] = (idata[tid] - idata[tid + iter * nx]) *
                 (idata[tid] - idata[tid + iter * nx]);
    tid += blockDim.x * gridDim.x;
  }
}

int norm_pscreen(float *d_odata, float *d_idata, int nx, int ny,
                 float norm_fact, carma_device *device) {
  float sfx, sfy, norm = 0;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, nx * ny, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  cublasHandle_t cublas_handle;
  cublasCreate(&cublas_handle);

  int npts = 5;
  for (int i = 1; i < npts + 1; i++) {
    carmaSafeCall(cudaMemset(d_odata, 0, sizeof(float) * nx * ny));
    dphix_krnl<<<grid, threads>>>(d_odata, d_idata, nx * ny, i, nx);
    carmaCheckMsg("dphix_kernel<<<>>> execution failed\n");
    // sfx  = cublasSasum(nx*ny,d_odata,1)/((nx-i)*ny);
    // here we can use asum because the initial array is positive (result of a
    // square)
    cublasSasum(cublas_handle, nx * ny, d_odata, 1, &sfx);

    carmaSafeCall(cudaMemset(d_odata, 0, sizeof(float) * nx * ny));
    dphiy_krnl<<<grid, threads>>>(d_odata, d_idata, nx * ny, i, nx);
    carmaCheckMsg("dphiy_kernel<<<>>> execution failed\n");
    // sfy  = cublasSasum(nx*ny,d_odata,1)/((ny-i)*nx);
    cublasSasum(cublas_handle, nx * ny, d_odata, 1, &sfy);

    // norm += sqrtf((sfx/((nx-i)*ny) +
    // sfy/((ny-i)*nx))/2.0f)/sqrtf(6.88*pow(i,1.66));
    norm += sqrtf(sfx / ((nx - i) * ny)) / sqrtf(6.88 * pow(i, 1.66));
  }
  norm /= npts;

  carmaSafeCall(cudaMemset(d_odata, 0, sizeof(float) * nx * ny));
  // cublasSaxpy(nx*ny,1.0f/norm*norm_fact, d_idata, 1, d_odata, 1);
  norm = (1.0f / norm) * norm_fact;
  cublasSaxpy(cublas_handle, nx * ny, &norm, d_idata, 1, d_odata, 1);

  cublasDestroy(cublas_handle);
  return EXIT_SUCCESS;
}
