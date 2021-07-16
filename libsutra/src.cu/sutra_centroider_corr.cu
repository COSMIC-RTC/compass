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

//! \file      sutra_centroider_corr.h
//! \ingroup   libsutra
//! \class     SutraCentroiderCorr
//! \brief     this class provides the centroider_corr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#include <sutra_centroider_corr.h>
#include <carma_utils.cuh>

template <class Tcu, class T>
__global__ void fillcorrcube_krnl(Tcu *d_out, T *d_in, int npix_in, int Npix_in,
                                  int npix_out, int Npix_out, int N) {
  int nim, npix, nlig, ncol, idx;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / Npix_in;
    npix = tid - nim * Npix_in;
    nlig = npix / npix_in;
    ncol = npix - nlig * npix_in;
    idx = nlig * npix_out + ncol + nim * Npix_out;
    d_out[idx].x = d_in[tid];
    d_out[idx].y = 0.0;
    tid += blockDim.x * gridDim.x;
  }
}

template <class Tcu, class T>
__global__ void fillcorrim_krnl(Tcu *d_out, T *d_in, int npix_in, int Npix_in,
                                int npix_out, int Npix_out, int N) {
  int nim, npix, nlig, ncol, idx;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / Npix_in;
    npix = tid - nim * Npix_in;
    nlig = npix / npix_in;
    ncol = npix - nlig * npix_in;
    idx = nlig * npix_out + ncol + nim * Npix_out;
    d_out[idx].x = d_in[npix];
    d_out[idx].y = 0.0;
    tid += blockDim.x * gridDim.x;
  }
}

template <class Tcu, class T>
int fillcorr(Tcu *d_out, T *d_in, int npix_in, int npix_out, int N, int nvalid,
             CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  if (nvalid == 1) {
    // cout << "3d" << endl;
    fillcorrcube_krnl<<<grid, threads>>>(d_out, d_in, npix_in,
                                         npix_in * npix_in, npix_out,
                                         npix_out * npix_out, N);
  } else {
    fillcorrim_krnl<<<grid, threads>>>(d_out, d_in, npix_in, npix_in * npix_in,
                                       npix_out, npix_out * npix_out, N);
    // cout << "2d" << endl;
  }

  carma_check_msg("fillcorr_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int fillcorr<cuFloatComplex, float>(cuFloatComplex *d_out, float *d_in,
                                             int npix_in, int npix_out, int N,
                                             int nvalid, CarmaDevice *device);
template int fillcorr<cuDoubleComplex, double>(cuDoubleComplex *d_out,
                                               double *d_in, int npix_in,
                                               int npix_out, int N, int nvalid,
                                               CarmaDevice *device);

template <class T>
__global__ void corr_krnl(T *odata, T *idata, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ T tmp;

  while (tid < N) {
    tmp.x = idata[tid].x * odata[tid].x + idata[tid].y * odata[tid].y;
    tmp.y = -1.0f * idata[tid].y * odata[tid].x + idata[tid].x * odata[tid].y;
    odata[tid].x = tmp.x;
    odata[tid].y = tmp.y;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int correl(T *d_odata, T *d_idata, int N, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  corr_krnl<<<grid, threads>>>(d_odata, d_idata, N);

  carma_check_msg("corr_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int correl<cuFloatComplex>(cuFloatComplex *d_odata,
                                    cuFloatComplex *d_idata, int N,
                                    CarmaDevice *device);
template int correl<cuDoubleComplex>(cuDoubleComplex *d_odata,
                                     cuDoubleComplex *d_idata, int N,
                                     CarmaDevice *device);

template <class Tcu, class T>
__global__ void roll2real_krnl(T *odata, Tcu *idata, int n, int Npix, int N) {
  // here we need to roll and keep only the (2:,2:,) elements
  // startegy is to go through all elements of input array
  // if final index > 0 (in x & y) keep it in output with (idx-1,idy-1)
  // n is always odd because it is 2 x npix
  int nim, idpix, idx, idy, idt;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / Npix;
    idpix = tid - nim * Npix;
    idy = idpix / n;
    idx = idpix - n * idy;

    idx = (idx + n / 2) % n;
    idy = (idy + n / 2) % n;

    if ((idx > 0) && (idy > 0)) {
      idt = idx - 1 + (idy - 1) * (n - 1) + nim * (n - 1) * (n - 1);
      odata[idt] = idata[tid].x;
    }
    tid += blockDim.x * gridDim.x;
  }
}

template <class Tcu, class T>
int roll2real(T *d_odata, Tcu *d_idata, int n, int Npix, int N,
              CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  roll2real_krnl<<<grid, threads>>>(d_odata, d_idata, n, Npix, N);

  carma_check_msg("roll2real_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int roll2real<cuFloatComplex, float>(float *d_odata,
                                              cuFloatComplex *d_idata, int n,
                                              int Npix, int N,
                                              CarmaDevice *device);
template int roll2real<cuDoubleComplex, double>(double *d_odata,
                                                cuDoubleComplex *d_idata, int n,
                                                int Npix, int N,
                                                CarmaDevice *device);

template <class T>
__global__ void corrnorm_krnl(T *odata, T *idata, int Npix, int N) {
  int nim, idpix;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ T tmp;

  while (tid < N) {
    nim = tid / Npix;
    idpix = tid - nim * Npix;
    tmp = odata[tid] / idata[idpix];
    odata[tid] = tmp;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int corr_norm(T *d_odata, T *d_idata, int Npix, int N, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  corrnorm_krnl<<<grid, threads>>>(d_odata, d_idata, Npix, N);

  carma_check_msg("corrnorm_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int corr_norm<float>(float *d_odata, float *d_idata, int Npix, int N,
                              CarmaDevice *device);
template int corr_norm<double>(double *d_odata, double *d_idata, int Npix,
                               int N, CarmaDevice *device);

template <class T>
__device__ inline void sortmaxi_krnl(T *sdata, unsigned int *values, int size,
                                     int n) {
  if (!((size & (size - 1)) == 0)) {
    unsigned int s;
    if ((size & 1) != 0)
      s = size / 2 + 1;  //(size&1)==size%2
    else
      s = size / 2;
    unsigned int s_old = size;
    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2 * s < s_old) && (s != 0)) s += 1;
    }
  } else {
    // do reduction in shared mem
    for (unsigned int s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
    }
  }
}

template <class T>
__global__ void sortmaxi(T *g_idata, int *values, int nmax, int offx, int offy,
                         int npix, int Npix, int Npix2) {
  extern __shared__ uint svalues[];
  T *sdata = (T *)&svalues[blockDim.x];

  // load shared mem
  unsigned int tid = threadIdx.x;

  svalues[tid] = tid;
  int nlig = tid / npix;
  int ncol = tid - nlig * npix;
  int idx = offx + ncol + (nlig + offy) * Npix + blockIdx.x * Npix2;
  sdata[tid] = g_idata[idx];

  __syncthreads();

  for (int cc = 0; cc < nmax; cc++) {
    if (tid >= cc)
      sortmaxi_krnl(&(sdata[cc]), &(svalues[cc]), blockDim.x - cc, tid - cc);

    __syncthreads();
  }

  if (tid < nmax) {
    nlig = svalues[tid] / npix;
    ncol = svalues[tid] - nlig * npix;
    idx = offx + ncol + (nlig + offy) * Npix;
    values[nmax * blockIdx.x + tid] = idx;
  }
  __syncthreads();
}

template <class T>
void subap_sortmaxi(int threads, int blocks, T *d_idata, int *values, int nmax,
                    int offx, int offy, int npix, int Npix)
// here idata is a [Npix,Npix,nvalid] array
// we want to get the [nvalid] max into subregions of [npix,npix] starting at
// [xoff,yoff] number of threads is npix * npix and number of blocks : nvalid
{
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = (threads <= 32) ? 2 * threads * (sizeof(T) + sizeof(uint))
                                 : threads * (sizeof(T) + sizeof(uint));
  sortmaxi<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, values, nmax, offx,
                                               offy, npix, Npix, Npix * Npix);

  carma_check_msg("sortmaxi_kernel<<<>>> execution failed\n");
}
template void subap_sortmaxi<float>(int threads, int blocks, float *d_idata,
                                    int *values, int nmax, int offx, int offy,
                                    int npix, int Npix);

template void subap_sortmaxi<double>(int threads, int blocks, double *d_idata,
                                     int *values, int nmax, int offx, int offy,
                                     int npix, int Npix);

template <class T>
__global__ void interp_parab(T *g_idata, T *g_centroids, int *g_values,
                             T *g_matinterp, int sizex, int sizey, int nvalid,
                             int Npix, int Npix2, float scale, float offset) {
  T *pidata = SharedMemory<T>();
  T *scoeff = (T *)&pidata[blockDim.x];
  T *m_interp = (T *)&scoeff[6];

  int offy = g_values[blockIdx.x] / Npix;
  int offx = g_values[blockIdx.x] - offy * Npix;
  offx -= sizex / 2;
  offy -= sizey / 2;
  int nlig = threadIdx.x / sizex;
  int ncol = threadIdx.x - nlig * sizex;
  int idx = offx + ncol + (nlig + offy) * Npix + blockIdx.x * Npix2;

  // load shared mem
  pidata[threadIdx.x] = g_idata[idx];

  __syncthreads();

  for (int cc = 0; cc < 6; cc++)
    m_interp[cc * sizex * sizey + threadIdx.x] =
        g_matinterp[cc * sizex * sizey + threadIdx.x] * pidata[threadIdx.x];

  __syncthreads();

  // do reduction for each 6 coeffs
  if (threadIdx.x < 6) {
    scoeff[threadIdx.x] = 0.0f;
    for (int cc = 0; cc < sizex * sizey; cc++)
      scoeff[threadIdx.x] += m_interp[cc + threadIdx.x * sizex * sizey];
  }

  __syncthreads();

  // now retreive x0 and y0 from scoeff
  if (threadIdx.x < 2) {
    T denom = scoeff[2] * scoeff[2] - 4.0f * scoeff[1] * scoeff[0];
    if (denom == 0) {
      if (threadIdx.x == 0) g_centroids[blockIdx.x] = 0.0f;
      if (threadIdx.x == 1) g_centroids[blockIdx.x + nvalid] = 0.0f;
    } else {
      if (threadIdx.x == 0) {
        g_centroids[blockIdx.x] =
            (2.0f * scoeff[1] * scoeff[3] - scoeff[4] * scoeff[2]) / denom;
        int xm = (2 * offx + sizex);
        xm = ((xm & 1) == 0) ? xm / 2 : xm / 2 + 1;  //(xm&1)==xm%2
        g_centroids[blockIdx.x] += (xm - 0.5 - (Npix + 1) / 4);
        g_centroids[blockIdx.x] = (g_centroids[blockIdx.x] - offset) * scale;
      }
      if (threadIdx.x == 1) {
        g_centroids[blockIdx.x + nvalid] =
            (2.0f * scoeff[0] * scoeff[4] - scoeff[3] * scoeff[2]) / denom;
        int ym = (2 * offy + sizey);
        ym = ((ym & 1) == 0) ? ym / 2 : ym / 2 + 1;  //(ym&1)==ym%2
        g_centroids[blockIdx.x + nvalid] += (ym - 0.5 - (Npix + 1) / 4);
        g_centroids[blockIdx.x + nvalid] =
            (g_centroids[blockIdx.x + nvalid] - offset) * scale;
      }
    }
  }

  __syncthreads();
  /*
   if (threadIdx.x == 0) g_centroids[blockIdx.x] = (g_centroids[blockIdx.x] -
   offset) * scale; if (threadIdx.x == 1) g_centroids[blockIdx.x+nvalid] =
   (g_centroids[blockIdx.x+nvalid] - offset) * scale;
   */
}

/*
 algorithm for parabolic interpolation
 we do the interpolation on nx x ny points
 we use a (nx * ny, 6) interp matrix
 we thus need nx * ny arrays in shared mem this is the number of threads
 we have nvalid blocks
 */

template <class T>
void subap_pinterp(int threads, int blocks, T *d_idata, int *values,
                   T *d_centroids, T *d_matinterp, int sizex, int sizey,
                   int nvalid, int Npix, float scale, float offset)
// here idata is a [Npix,Npix,nvalid] array
// we want to get the [nvalid] (x0,y0) into subregions of [sizex,sizey] around
// gvalue number of threads is sizex * sizey  and number of blocks : nvalid
{
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  int shsize = (threads + threads * 6 + 6);
  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds

  int smemSize = (shsize <= 32) ? 2 * shsize * sizeof(T) : shsize * sizeof(T);
  interp_parab<T><<<dimGrid, dimBlock, smemSize>>>(
      d_idata, d_centroids, values, d_matinterp, sizex, sizey, nvalid, Npix,
      Npix * Npix, scale, offset);

  carma_check_msg("sortmaxi_kernel<<<>>> execution failed\n");
}

template void subap_pinterp<float>(int threads, int blocks, float *d_idata,
                                   int *values, float *d_centroids,
                                   float *d_matinterp, int sizex, int sizey,
                                   int nvalid, int Npix, float scale,
                                   float offset);

template void subap_pinterp<double>(int threads, int blocks, double *d_idata,
                                    int *values, double *d_centroids,
                                    double *d_matinterp, int sizex, int sizey,
                                    int nvalid, int Npix, float scale,
                                    float offset);
