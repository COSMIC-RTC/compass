/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <cutil.h>
#include <carma_obj.hpp>
#include "convolutionFFT2D.cuh"
#include "convolutionFFT2D_common.h"

////////////////////////////////////////////////////////////////////////////////
/// Position convolution kernel center at (0, 0) in the image
////////////////////////////////////////////////////////////////////////////////
extern "C" void pad_kernel(float *d_Dst, float *d_Src, int32_t fftH, int32_t fftW,
                          int32_t kernelH, int32_t kernelW, int32_t kernelY, int32_t kernelX) {
  assert(d_Src != d_Dst);
  dim3 threads(32, 8);
  dim3 grid(i_div_up(kernelW, threads.x), i_div_up(kernelH, threads.y));

  SET_FLOAT_BASE;
  padKernel_kernel<<<grid, threads>>>(d_Dst, d_Src, fftH, fftW, kernelH,
                                      kernelW, kernelY, kernelX);
  carma_check_msg("padKernel_kernel<<<>>> execution failed\n");
}

extern "C" void pad_kernel_3d(float *d_Dst, float *d_Src, int32_t fftH, int32_t fftW,
                            int32_t kernelH, int32_t kernelW, int32_t kernelY, int32_t kernelX,
                            int32_t nim) {
  assert(d_Src != d_Dst);
  dim3 threads(16, 8, 8);
  dim3 grid(i_div_up(kernelW, threads.x), i_div_up(kernelH, threads.y),
            i_div_up(nim, threads.z));
  // dim3 grid(i_div_up(kernelW, threads.x), i_div_up(kernelH, threads.y),nim);

  SET_FLOAT_BASE;
  pad_kernel_3d_kernel<<<grid, threads>>>(d_Dst, d_Src, fftH, fftW, kernelH,
                                        kernelW, kernelY, kernelX, nim);

  carma_check_msg("pad_kernel_3d_kernel<<<>>> execution failed\n");
}

////////////////////////////////////////////////////////////////////////////////
// Prepare data for "pad to border" addressing mode
////////////////////////////////////////////////////////////////////////////////
extern "C" void pad_data_clamp_to_border(float *d_Dst, float *d_Src, int32_t fftH,
                                     int32_t fftW, int32_t dataH, int32_t dataW,
                                     int32_t kernelW, int32_t kernelH, int32_t kernelY,
                                     int32_t kernelX) {
  assert(d_Src != d_Dst);
  dim3 threads(32, 8);
  dim3 grid(i_div_up(fftW, threads.x), i_div_up(fftH, threads.y));

  SET_FLOAT_BASE;
  pad_data_clamp_to_border_kernel<<<grid, threads>>>(d_Dst, d_Src, fftH, fftW,
                                                 dataH, dataW, kernelH, kernelW,
                                                 kernelY, kernelX);
  carma_check_msg("pad_data_clamp_to_border_kernel<<<>>> execution failed\n");
}

extern "C" void pad_data_clamp_to_border_3d(float *d_Dst, float *d_Src, int32_t fftH,
                                       int32_t fftW, int32_t dataH, int32_t dataW,
                                       int32_t kernelW, int32_t kernelH, int32_t kernelY,
                                       int32_t kernelX, int32_t nim) {
  assert(d_Src != d_Dst);
  dim3 threads(16, 8, 8);
  dim3 grid(i_div_up(fftW, threads.x), i_div_up(fftH, threads.y),
            i_div_up(nim, threads.z));

  SET_FLOAT_BASE;
  pad_data_clamp_to_border_3d_kernel<<<grid, threads>>>(
      d_Dst, d_Src, fftH, fftW, dataH, dataW, kernelH, kernelW, kernelY,
      kernelX, nim);
  carma_check_msg("pad_data_clamp_to_border_3d_kernel<<<>>> execution failed\n");
}

////////////////////////////////////////////////////////////////////////////////
// Modulate Fourier image of padded data by Fourier image of padded kernel
// and normalize by FFT size
////////////////////////////////////////////////////////////////////////////////
extern "C" void modulate_and_normalize(fComplex *d_Dst, fComplex *d_Src, int32_t fftH,
                                     int32_t fftW, int32_t padding, int32_t nim) {
  assert(fftW % 2 == 0);
  const int32_t dataSize = fftH * (fftW / 2 + padding) * nim;

  modulate_and_normalize_kernel<<<i_div_up(dataSize, 256), 256>>>(
      d_Dst, d_Src, dataSize, 1.0f / (float)(fftW * fftH));
  carma_check_msg("modulate_and_normalize() execution failed\n");
}

////////////////////////////////////////////////////////////////////////////////
// 2D R2C / C2R post/preprocessing kernels
////////////////////////////////////////////////////////////////////////////////
static const double PI = 3.1415926535897932384626433832795;
static const uint BLOCKDIM = 256;

extern "C" void sp_postprocess_2d(void *d_Dst, void *d_Src, uint DY, uint DX,
                                uint padding, int32_t dir) {
  assert(d_Src != d_Dst);
  assert(DX % 2 == 0);

#if (POWER_OF_TWO)
  uint log2DX, log2DY;
  uint factorizationRemX = factorRadix2(log2DX, DX);
  uint factorizationRemY = factorRadix2(log2DY, DY);
  assert(factorizationRemX == 1 && factorizationRemY == 1);
#endif

  const uint threadCount = DY * (DX / 2);
  const double phaseBase = dir * PI / (double)DX;

  SET_FCOMPLEX_BASE;
  sp_postprocess_2d_kernel<<<i_div_up(threadCount, BLOCKDIM), BLOCKDIM>>>(
      (fComplex *)d_Dst, (fComplex *)d_Src, DY, DX, threadCount, padding,
      (float)phaseBase);
  carma_check_msg("sp_postprocess_2d_kernel<<<>>> execution failed\n");
}

extern "C" void sp_preprocess_2d(void *d_Dst, void *d_Src, uint DY, uint DX,
                               uint padding, int32_t dir) {
  assert(d_Src != d_Dst);
  assert(DX % 2 == 0);

#if (POWER_OF_TWO)
  uint log2DX, log2DY;
  uint factorizationRemX = factorRadix2(log2DX, DX);
  uint factorizationRemY = factorRadix2(log2DY, DY);
  assert(factorizationRemX == 1 && factorizationRemY == 1);
#endif

  const uint threadCount = DY * (DX / 2);
  const double phaseBase = -dir * PI / (double)DX;

  SET_FCOMPLEX_BASE;
  sp_preprocess_2d_kernel<<<i_div_up(threadCount, BLOCKDIM), BLOCKDIM>>>(
      (fComplex *)d_Dst, (fComplex *)d_Src, DY, DX, threadCount, padding,
      (float)phaseBase);
  carma_check_msg("sp_preprocess_2d_kernel<<<>>> execution failed\n");
}

////////////////////////////////////////////////////////////////////////////////
// Combined sp_postprocess_2d + modulate_and_normalize + sp_preprocess_2d
////////////////////////////////////////////////////////////////////////////////
extern "C" void sp_process_2d(void *d_Dst, void *d_SrcA, void *d_SrcB, uint DY,
                            uint DX, int32_t dir) {
  assert(DY % 2 == 0);

#if (POWER_OF_TWO)
  uint log2DX, log2DY;
  uint factorizationRemX = factorRadix2(log2DX, DX);
  uint factorizationRemY = factorRadix2(log2DY, DY);
  assert(factorizationRemX == 1 && factorizationRemY == 1);
#endif

  const uint threadCount = (DY / 2) * DX;
  const double phaseBase = dir * PI / (double)DX;

  SET_FCOMPLEX_BASE_A;
  SET_FCOMPLEX_BASE_B;
  sp_process_2d_kernel<<<i_div_up(threadCount, BLOCKDIM), BLOCKDIM>>>(
      (fComplex *)d_Dst, (fComplex *)d_SrcA, (fComplex *)d_SrcB, DY, DX,
      threadCount, (float)phaseBase, 0.5f / (float)(DY * DX));
  carma_check_msg("sp_process_2d_kernel<<<>>> execution failed\n");
}
