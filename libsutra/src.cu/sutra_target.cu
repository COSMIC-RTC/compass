// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_target.cu
//! \ingroup   libsutra
//! \class     SutraTarget
//! \brief     this class provides the target features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24

#include <sutra_target.h>

__device__ void generic_raytrace(float *odata, float *idata, int nx, int ny,
                                 float xoff, float yoff, float G, float thetaML,
                                 float dx, float dy, int Nx, int blockSize,
                                 int istart, float delta) {
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  y += istart;

  int tido;

  int iref, jref;  // tidi;
  float x1 = ((x - nx / 2.0f) * G);
  float y1 = ((y - ny / 2.0f) * G);

  float x2 = (cosf(thetaML) * x1 - sinf(thetaML) * y1) + nx / 2.0f + xoff / delta;
  float y2 = (sinf(thetaML) * x1 + cosf(thetaML) * y1) + ny / 2.0f + yoff / delta;

  float xref = x2 * delta + dx;
  float yref = y2 * delta + dy;

  float xshift, yshift, wx1, wx2, wy1, wy2;

  iref = (int)xref;
  jref = (int)yref;

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

__global__ void raytrace_krnl(float *odata, float *idata, int nx, int ny,
                              float xoff, float yoff, float G, float thetaML,
                              float dx, float dy, int Nx, int blockSize, float delta) {
  generic_raytrace(odata, idata, nx, ny, xoff, yoff, G, thetaML, dx, dy, Nx,
                   blockSize, 0, delta);
}


__global__ void raytrace_krnl(float *odata, float *idata, int nx, int ny,
                              float xoff, float yoff, int Nx, int blockSize,
                              int istart, float delta) {
  generic_raytrace(odata, idata, nx, ny, xoff, yoff, 1.0f, 0.0f, 0.0f, 0.0f, Nx,
                   blockSize, istart, delta);
}

int target_raytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
                    float xoff, float yoff, float G, float thetaML, float dx,
                    float dy, int block_size, float delta) {
  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  int nny = ny + block_size - ny % block_size;
  dim3 blocks(nnx / block_size, nny / block_size),
      threads(block_size, block_size);

  int smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  raytrace_krnl<<<blocks, threads, smemSize>>>(
      d_odata, d_idata, nx, ny, xoff, yoff, G, thetaML, dx, dy, Nx, block_size, delta);

  carma_check_msg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}


int target_raytrace_async(CarmaStreams *streams, float *d_odata,
                          float *d_idata, int nx, int ny, int Nx, float xoff,
                          float yoff, int block_size) {
  int nstreams = streams->get_nb_streams();

  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  dim3 blocks(nnx / block_size, 1), threads(block_size, block_size);
  int smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  for (int i = 0; i < nstreams; i++)
    raytrace_krnl<<<blocks, threads, smemSize, streams->get_stream(i)>>>(
        d_odata, d_idata, nx, ny, xoff, yoff, Nx, block_size, i * block_size, 1.0f);

  carma_check_msg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int target_raytrace_async(CarmaHostObj<float> *phase_telemetry,
                          float *d_odata, float *d_idata, int nx, int ny,
                          int Nx, float xoff, float yoff, int block_size) {
  float *hdata = phase_telemetry->get_data();
  int nstreams = phase_telemetry->get_nb_streams();

  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  dim3 blocks(nnx / block_size, 1), threads(block_size, block_size);
  int smemSize = (block_size + 1) * (block_size + 1) * sizeof(float);

  for (int i = 0; i < nstreams; i++)
    raytrace_krnl<<<blocks, threads, smemSize,
                    phase_telemetry->get_cuda_stream(i)>>>(
        d_odata, d_idata, nx, ny, xoff, yoff, Nx, block_size, i * block_size, 1.0f);

  carma_check_msg("raytrace_kernel<<<>>> execution failed\n");
  // asynchronously launch nstreams memcopies.  Note that memcopy in stream x
  // will only
  //   commence executing when all previous CUDA calls in stream x have
  //   completed
  int delta = block_size * nx;
  for (int i = 0; i < nstreams; i++) {
    int nbcopy = nx * ny - (i + 1) * block_size * nx;
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
                              float *mask, float scale, int puponly, int nx,
                              int Np, int Nx) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Np) {
    int nlinep = tid / nx;
    int ncol = tid - nlinep * nx;
    int nim = ncol + nlinep * Nx;

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

int fill_amplipup(cuFloatComplex *amplipup, float *phase, float *mask,
                  float scale, int puponly, int nx, int ny, int Nx,
                  CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nx * ny, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fillamplikrnl<<<grid, threads>>>(amplipup, phase, mask, scale, puponly, nx,
                                   nx * ny, Nx);
  carma_check_msg("fillamplikrnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
