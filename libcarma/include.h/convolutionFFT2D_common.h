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

typedef unsigned int uint;

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
inline int i_div_up(int a, int b) { return (a % b != 0) ? (a / b + 1) : (a / b); }

// Align a to nearest higher multiple of b
inline int i_align_up(int a, int b) { return (a % b != 0) ? (a - a % b + b) : a; }

extern "C" {

void pad_kernel(float *d_PaddedKernel, float *d_Kernel, int fftH, int fftW,
               int kernelH, int kernelW, int kernelY, int kernelX);

void pad_kernel_3d(float *d_PaddedKernel, float *d_Kernel, int fftH, int fftW,
                 int kernelH, int kernelW, int kernelY, int kernelX, int nim);

void pad_data_clamp_to_border(float *d_PaddedData, float *d_Data, int fftH,
                          int fftW, int dataH, int dataW, int kernelH,
                          int kernelW, int kernelY, int kernelX);

void pad_data_clamp_to_border_3d(float *d_PaddedData, float *d_Data, int fftH,
                            int fftW, int dataH, int dataW, int kernelH,
                            int kernelW, int kernelY, int kernelX, int nim);

void modulate_and_normalize(fComplex *d_Dst, fComplex *d_Src, int fftH, int fftW,
                          int padding, int nim);

void sp_postprocess_2d(void *d_Dst, void *d_Src, uint DY, uint DX, uint padding,
                     int dir);

void sp_preprocess_2d(void *d_Dst, void *d_Src, uint DY, uint DX, uint padding,
                    int dir);

void sp_process_2d(void *d_Data, void *d_Data0, void *d_Kernel0, uint DY, uint DX,
                 int dir);
}
#endif  // CONVOLUTIONFFT2D_COMMON_H
