// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_lgs.cu
//! \ingroup   libsutra
//! \class     SutraLGS
//! \brief     this class provides the lgs features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_lgs.h>

// texture<float, 3, cudaReadModeElementType> tex3;  // 3D texture

__device__ float interp_transf(float *idata, float tidim, int32_t npix,
                               float pixsize, float norm, float shiftx,
                               float deltah, int32_t hmax) {
  int32_t xref;
  float weightx;

  // determine the interpolated position
  if (npix % 2 == 0)
    tidim -= (npix / 2.0f + 0.5);  // centered pix pos
  else
    tidim -= (npix / 2.0f);  // centered pix pos

  tidim *= pixsize;  // centered pix size
  tidim /= norm;     // in meter
  tidim += shiftx;   // in meter on profile centered on cog
  tidim /= deltah;   // in pixels on profile array

  if (tidim < 0) {
    xref = 0;
    weightx = 0.;
  } else {
    if (tidim < hmax - 2) {
      xref = (int32_t)tidim;
      weightx = tidim - xref;
    } else {
      xref = hmax - 2;
      weightx = 1.;
    }
  }

  return (idata[xref] * (1 - weightx) + idata[xref + 1] * weightx);
}

__global__ void iprof_krnl(cuFloatComplex *profout, float *profin,
                           float *profinc, int32_t npix, float *doffaxis, float hg,
                           float pixsize, float h0, float deltah, int32_t hmax,
                           int32_t Ntot) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Ntot) {
    int32_t nim = tid / npix;
    float tidim = tid - nim * npix;
    float norm;

    norm = 206265.0f * doffaxis[nim] / (hg * hg);

    if (pixsize > deltah * norm) {
      profout[tid].x = interp_transf(profinc, tidim + 1, npix, pixsize, norm,
                                     hg - h0, deltah, hmax) -
                       interp_transf(profinc, tidim, npix, pixsize, norm,
                                     hg - h0, deltah, hmax);
      profout[tid].y = 0.0f;
    } else {
      profout[tid].x = interp_transf(profin, tidim, npix, pixsize, norm,
                                     hg - h0, deltah, hmax);
      profout[tid].y = 0.0f;
    }

    tid += blockDim.x * gridDim.x;
  }
}

int32_t interp_prof(cuFloatComplex *profout, float *prof1d, float *profcum,
                int32_t npix, float *doffaxis, float hg, float pixsize, float h0,
                float deltah, int32_t hmax, int32_t Ntot, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, Ntot, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  iprof_krnl<<<grid, threads>>>(profout, prof1d, profcum, npix, doffaxis, hg,
                                pixsize, h0, deltah, hmax, Ntot);
  carma_check_msg("iprof_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void tftbeam_krnl(cuFloatComplex *profout, cuFloatComplex *fbeam,
                             int32_t N, int32_t Ntot) {
  int32_t idim;
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Ntot) {
    idim = tid % N;
    cuFloatComplex tmp = profout[tid];
    profout[tid].x = tmp.x * fbeam[idim].x - tmp.y * fbeam[idim].y;
    profout[tid].y = tmp.y * fbeam[idim].x + tmp.x * fbeam[idim].y;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t times_ftbeam(cuFloatComplex *profout, cuFloatComplex *ftbeam, int32_t N,
                 int32_t Ntot, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, Ntot, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  tftbeam_krnl<<<grid, threads>>>(profout, ftbeam, N, Ntot);
  carma_check_msg("tftbeam_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void rollbeamexp_krnl(float *imout, cuFloatComplex *iprof,
                                 float *beam, int32_t N, int32_t Ntot) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Ntot) {
    int32_t nim = tid / (N * N);
    int32_t tidim = tid % (N * N);
    int32_t tidprof = tidim % N;
    int32_t tidbeam = tidim / N;
    int32_t x = (tidprof + (N / 2)) % N + nim * N;  // for roll
    imout[tid] = abs(iprof[x].x) * beam[tidbeam];
    tid += blockDim.x * gridDim.x;
  }
}

int32_t roll_beam_exp(float *imout, cuFloatComplex *iprof, float *beam, int32_t N,
                int32_t Ntot, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, Ntot, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  rollbeamexp_krnl<<<grid, threads>>>(imout, iprof, beam, N, Ntot);
  carma_check_msg("beamexp_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void rotate_krnl(cuFloatComplex *odata, float *idata, int32_t width,
                            int32_t height, float *theta, float center, int32_t N,
                            int32_t Ntot)
// rotate a cube manually
{
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < Ntot) {
    int32_t nim = tid / N;
    int32_t tidim = tid - nim * N;

    int32_t y = tidim / width;
    int32_t x = tidim - y * width;

    float ucent = width / 2 - center;
    float vcent = height / 2 - center;

    float u = x - ucent;
    float v = y - vcent;
    // transform coordinates
    float tu = u * cos(theta[nim]) + v * sin(theta[nim]) + ucent;
    float tv = -u * sin(theta[nim]) + v * cos(theta[nim]) + ucent;

    int32_t uref, vref;
    float weightu, weightv;

    if (tu < 0) {
      uref = 0;
      weightu = 0.;
    } else {
      if (tu < width - 2) {
        uref = (int32_t)tu;
        weightu = tu - uref;
      } else {
        uref = width - 2;
        weightu = 1.;
      }
    }

    if (tv < 0) {
      vref = 0;
      weightv = 0.;
    } else {
      if (tv < height - 2) {
        vref = (int32_t)tv;
        weightv = tv - vref;
      } else {
        vref = height - 2;
        weightv = 1.;
      }
    }

    int32_t ind1, ind2, ind3, ind4;
    ind1 = ind2 = ind3 = ind4 = nim * N;

    ind1 += uref + vref * width;
    ind2 += uref + 1 + vref * width;
    ind3 += uref + (vref + 1) * width;
    ind4 += uref + 1 + (vref + 1) * width;

    odata[tid].x =
        (idata[ind1] * (1 - weightu) + idata[ind2] * weightu) * (1 - weightv) +
        (idata[ind3] * (1 - weightu) + idata[ind4] * weightu) * weightv;
    odata[tid].y = 0.0f;

    tid += blockDim.x * gridDim.x;
  }
}

int32_t lgs_rotate(cuFloatComplex *odata, float *idata, int32_t width, int32_t height,
               float *theta, float center, int32_t Ntot, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, Ntot, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  int32_t N = width * height;

  rotate_krnl<<<grid, threads>>>(odata, idata, width, height, theta, center, N,
                                 Ntot);
  carma_check_msg("rotate_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

// __global__ void rotate3d_krnl(cuFloatComplex *g_odata, int32_t width, int32_t height,
//                               int32_t N, float *theta, float center, int32_t Ntot)
// // rotate a cube using 3d texture fetch
// {
//   int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

//   while (tid < Ntot) {
//     int32_t nim = tid / N;
//     int32_t tidim = tid - nim * N;

//     int32_t y = tidim / width;
//     int32_t x = tidim - y * width;

//     float ucent = width / 2 - center;
//     float vcent = height / 2 - center;

//     float u = x - ucent;
//     float v = y - vcent;
//     // transform coordinates
//     float tu = u * cos(theta[nim]) + v * sin(theta[nim]) + ucent;
//     float tv = -u * sin(theta[nim]) + v * cos(theta[nim]) + ucent;

//     g_odata[tid].x = tex3D(tex3, tu + 0.5, tv + 0.5, nim + 0.5);
//     g_odata[tid].y = 0.0f;
//     tid += blockDim.x * gridDim.x;
//   }
// }

// int32_t rotate3d(cuFloatComplex *d_odata, cudaMemcpy3DParms copyParams,
//              cudaArray *d_array, cudaChannelFormatDesc channel_desc, int32_t width,
//              int32_t height, float *theta, float center, int32_t Ntot,
//              CarmaDevice *device) {
//   tex3.normalized = false;
//   tex3.filterMode = cudaFilterModeLinear;      // linear interpolation
//   tex3.addressMode[0] = cudaAddressModeClamp;  // wrap texture coordinates
//   tex3.addressMode[1] = cudaAddressModeClamp;
//   tex3.addressMode[2] = cudaAddressModeClamp;

//   // copy the data into the array
//   carma_safe_call(cudaMemcpy3D(&copyParams));
//   // bind array to 3D texture
//   carma_safe_call(cudaBindTextureToArray(tex3, d_array, channel_desc));

//   int32_t N = width * height;

//   int32_t nb_blocks, nb_threads;
//   get_num_blocks_and_threads(device, Ntot, nb_blocks, nb_threads);
//   dim3 grid(nb_blocks), threads(nb_threads);

//   // cout << nb_blocks << " " << nx / nb_blocks << " " << ny / nb_blocks << " " <<  nz
//   // / nb_blocks <<endl;
//   rotate3d_krnl<<<grid, threads, 0>>>(d_odata, width, height, N, theta, center,
//                                       Ntot);

//   carma_check_msg("rotate3d_krnl <<<>>> execution failed\n");

//   return EXIT_SUCCESS;
// }
