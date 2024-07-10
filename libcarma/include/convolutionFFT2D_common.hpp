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

#ifndef CONVOLUTIONFFT2D_COMMON_H
#define CONVOLUTIONFFT2D_COMMON_H

typedef uint32_t uint;

#ifdef __CUDACC__
typedef float2 fComplex;
#else
typedef struct {
  float x;
  float y;
} fComplex;
#endif

////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////
// Round a / b to nearest higher integer value
inline int32_t i_div_up(int32_t a, int32_t b) { return (a % b != 0) ? (a / b + 1) : (a / b); }

// Align a to nearest higher multiple of b
inline int32_t i_align_up(int32_t a, int32_t b) { return (a % b != 0) ? (a - a % b + b) : a; }

extern "C" {

void pad_kernel(float *d_PaddedKernel, float *d_Kernel, int32_t fftH, int32_t fftW,
               int32_t kernelH, int32_t kernelW, int32_t kernelY, int32_t kernelX);

void pad_kernel_3d(float *d_PaddedKernel, float *d_Kernel, int32_t fftH, int32_t fftW,
                 int32_t kernelH, int32_t kernelW, int32_t kernelY, int32_t kernelX, int32_t nim);

void pad_data_clamp_to_border(float *d_PaddedData, float *d_Data, int32_t fftH,
                          int32_t fftW, int32_t dataH, int32_t dataW, int32_t kernelH,
                          int32_t kernelW, int32_t kernelY, int32_t kernelX);

void pad_data_clamp_to_border_3d(float *d_PaddedData, float *d_Data, int32_t fftH,
                            int32_t fftW, int32_t dataH, int32_t dataW, int32_t kernelH,
                            int32_t kernelW, int32_t kernelY, int32_t kernelX, int32_t nim);

void modulate_and_normalize(fComplex *d_Dst, fComplex *d_Src, int32_t fftH, int32_t fftW,
                          int32_t padding, int32_t nim);

void sp_postprocess_2d(void *d_Dst, void *d_Src, uint DY, uint DX, uint padding,
                     int32_t dir);

void sp_preprocess_2d(void *d_Dst, void *d_Src, uint DY, uint DX, uint padding,
                    int32_t dir);

void sp_process_2d(void *d_Data, void *d_Data0, void *d_Kernel0, uint DY, uint DX,
                 int32_t dir);
}
#endif  // CONVOLUTIONFFT2D_COMMON_H
