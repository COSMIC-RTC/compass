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

//! \file      sutra_target.cu
//! \ingroup   libsutra
//! \class     SutraTarget
//! \brief     this class provides the target features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_target.hpp>

__device__ void generic_raytrace(float *odata, float *idata, int32_t nx, int32_t ny,
                                 float xoff, float yoff, float G, float thetaML,
                                 float dx, float dy, int32_t Nx, int32_t blockSize,
                                 int32_t istart, float delta) {
  int32_t x = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t y = threadIdx.y + blockIdx.y * blockDim.y;
  y += istart;

  int32_t tido;

  int32_t iref, jref;  // tidi;
  float x1 = ((x - nx / 2.0f) * G);
  float y1 = ((y - ny / 2.0f) * G);

  float x2 = (cosf(thetaML) * x1 - sinf(thetaML) * y1) + nx / 2.0f + xoff / delta;
  float y2 = (sinf(thetaML) * x1 + cosf(thetaML) * y1) + ny / 2.0f + yoff / delta;

  float xref = x2 * delta + dx;
  float yref = y2 * delta + dy;

  float xshift, yshift, wx1, wx2, wy1, wy2;

  iref = (int32_t)xref;
  jref = (int32_t)yref;

  if ((x < nx) && (y < ny)) {
    tido = x + y * nx;

    xshift = xref - iref;
    yshift = yref - jref;

    wx1 = (1.0f - xshift);
    wx2 = xshift;
    wy1 = (1.0f - yshift);
    wy2 = yshift;

    if ((iref + 1 < Nx) && (jref + 1 < Nx)) {
      odata[tido] += (wx1 * wy1 * idata[iref + jref * Nx] +
                      wx2 * wy1 * idata[iref + 1 + jref * Nx] +
                      wx1 * wy2 * idata[iref + (jref + 1) * Nx] +
                      wx2 * wy2 * idata[iref + 1 + (jref + 1) * Nx]);
    } else {
      odata[tido] += 0.0f;
    }
  }
}

__global__ void raytrace_krnl(float *odata, float *idata, int32_t nx, int32_t ny,
                              float xoff, float yoff, float G, float thetaML,
                              float dx, float dy, int32_t Nx, int32_t blockSize, float delta) {
  generic_raytrace(odata, idata, nx, ny, xoff, yoff, G, thetaML, dx, dy, Nx,
                   blockSize, 0, delta);
}


__global__ void raytrace_krnl(float *odata, float *idata, int32_t nx, int32_t ny,
                              float xoff, float yoff, int32_t Nx, int32_t blockSize,
                              int32_t istart, float delta) {
  generic_raytrace(odata, idata, nx, ny, xoff, yoff, 1.0f, 0.0f, 0.0f, 0.0f, Nx,
                   blockSize, istart, delta);
}

int32_t target_raytrace(float *d_odata, float *d_idata, int32_t nx, int32_t ny, int32_t Nx,
                    float xoff, float yoff, float G, float thetaML, float dx,
                    float dy, int32_t block_size, float delta) {
  int32_t nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  int32_t nny = ny + block_size - ny % block_size;
  dim3 blocks(nnx / block_size, nny / block_size),
      threads(block_size, block_size);

  int32_t smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  raytrace_krnl<<<blocks, threads, smemSize>>>(
      d_odata, d_idata, nx, ny, xoff, yoff, G, thetaML, dx, dy, Nx, block_size, delta);

  carma_check_msg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}


int32_t target_raytrace_async(CarmaStreams *streams, float *d_odata,
                          float *d_idata, int32_t nx, int32_t ny, int32_t Nx, float xoff,
                          float yoff, int32_t block_size) {
  int32_t nstreams = streams->get_nb_streams();

  int32_t nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  dim3 blocks(nnx / block_size, 1), threads(block_size, block_size);
  int32_t smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  for (int32_t i = 0; i < nstreams; i++)
    raytrace_krnl<<<blocks, threads, smemSize, streams->get_stream(i)>>>(
        d_odata, d_idata, nx, ny, xoff, yoff, Nx, block_size, i * block_size, 1.0f);

  carma_check_msg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int32_t target_raytrace_async(CarmaHostObj<float> *phase_telemetry,
                          float *d_odata, float *d_idata, int32_t nx, int32_t ny,
                          int32_t Nx, float xoff, float yoff, int32_t block_size) {
  float *hdata = phase_telemetry->get_data();
  int32_t nstreams = phase_telemetry->get_nb_streams();

  int32_t nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  dim3 blocks(nnx / block_size, 1), threads(block_size, block_size);
  int32_t smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  for (int32_t i = 0; i < nstreams; i++)
    raytrace_krnl<<<blocks, threads, smemSize,
                    phase_telemetry->get_cuda_stream(i)>>>(
        d_odata, d_idata, nx, ny, xoff, yoff, Nx, block_size, i * block_size, 1.0f);

  carma_check_msg("raytrace_kernel<<<>>> execution failed\n");
  // asynchronously launch nstreams memcopies.  Note that memcopy in stream x
  // will only
  //   commence executing when all previous CUDA calls in stream x have
  //   completed
  int32_t delta = block_size * nx;
  for (int32_t i = 0; i < nstreams; i++) {
    int32_t nbcopy = nx * ny - (i + 1) * block_size * nx;
    if (nbcopy > 0) {
      nbcopy = delta;
    } else
      nbcopy = delta + nbcopy;
    cudaMemcpyAsync(&(hdata[i * block_size * nx]),
                    &(d_odata[i * block_size * nx]), sizeof(float) * nbcopy,
                    cudaMemcpyDeviceToHost,
                    phase_telemetry->get_cuda_stream(i));
  }

  return EXIT_SUCCESS;
}

__global__ void fillamplikrnl(cuFloatComplex *amplipup, float *phase,
                              float *mask, float scale, int32_t puponly, int32_t nx,
                              int32_t Np, int32_t Nx) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Np) {
    int32_t nlinep = tid / nx;
    int32_t ncol = tid - nlinep * nx;
    int32_t nim = ncol + nlinep * Nx;

    if (puponly == 1) {
      amplipup[nim].x = mask[tid] > 0 ? 1.0f : 0.0f;
      amplipup[nim].y = 0.0f;
    } else {
      if (mask[tid] > 0) {
          amplipup[nim].x = cosf(scale * phase[tid]);
          amplipup[nim].y = sinf(scale * phase[tid]);
      }
      else {
        amplipup[nim].x = 0.0f;
        amplipup[nim].y = 0.0f;
      }
    }
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fill_amplipup(cuFloatComplex *amplipup, float *phase, float *mask,
                  float scale, int32_t puponly, int32_t nx, int32_t ny, int32_t Nx,
                  CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nx * ny, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fillamplikrnl<<<grid, threads>>>(amplipup, phase, mask, scale, puponly, nx,
                                   nx * ny, Nx);
  carma_check_msg("fillamplikrnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
