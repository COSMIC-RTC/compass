// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_centroider_tcog.cu
//! \ingroup   libsutra
//! \class     SutraCentroiderTcog
//! \brief     this class provides the centroider_tcog features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_centroider_tcog.hpp>
#include <sutra_centroider_utils.cuh>
#include <carma_utils.cuh>

template <int32_t nb_threads, typename T>
__global__ void centroids(float *d_img, T *d_centroids, T *ref, int32_t *validx,
                          int32_t *validy, float *d_intensities, float threshold,
                          uint32_t npix, sutra::SlopesIndex si, uint32_t size, T scale,
                          T offset, uint32_t nelem_thread) {
  if (blockDim.x > nb_threads) {
    if (threadIdx.x == 0) printf("Wrong size argument\n");
    return;
  }
  // Specialize BlockReduce for a 1D block of 128 threads on type int32_t
  typedef cub::BlockReduce<float, nb_threads> BlockReduce;
  // Allocate shared memory for BlockReduce
  __shared__ typename BlockReduce::TempStorage temp_storage;

  float idata = 0;
  float xdata = 0;
  float ydata = 0;
  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t xvalid = validx[blockIdx.x];
  uint32_t yvalid = validy[blockIdx.x];
  uint32_t x, y;
  int32_t idim;

  for (int32_t cc = 0; cc < nelem_thread; cc++) {
    x = ((tid * nelem_thread + cc) % npix);
    y = ((tid * nelem_thread + cc) / npix);
    // idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) *
    // blockIdx.x;
    idim = (x + xvalid) + (y + yvalid) * size;
    if (idim < size * size) {
      float data_thresh =
          (d_img[idim] > threshold) ? d_img[idim] - threshold : 1.e-6;
      idata += data_thresh;
      xdata += data_thresh * x;
      ydata += data_thresh * y;
      d_img[idim] = data_thresh;
    }
  }

  // sdata[tid] = (i < N) ? g_idata[i] * x : 0;
  __syncthreads();
  float intensity = BlockReduce(temp_storage).Sum(idata, npix * npix);
  __syncthreads();
  float slopex = BlockReduce(temp_storage).Sum(xdata, npix * npix);
  __syncthreads();
  float slopey = BlockReduce(temp_storage).Sum(ydata, npix * npix);

  // write result for this block to global mem
  if (tid == 0) {
    d_centroids[si.x(blockIdx.x)] = (T(slopex / (intensity + 1.e-6)) - offset) * scale - ref[si.x(blockIdx.x)];
    d_centroids[si.y(blockIdx.x)] = (T(slopey / (intensity + 1.e-6)) - offset) * scale - ref[si.y(blockIdx.x)];
    d_intensities[blockIdx.x] = intensity;
  }
}

template <class T>
void get_centroids(int32_t size, int32_t threads, int32_t blocks, int32_t npix, float *d_img,
                   T *d_centroids, T *ref, int32_t *validx, int32_t *validy,
                   float *intensities, float threshold, float scale,
                   float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device, cudaStream_t stream) {
  int32_t maxThreads = device->get_properties().maxThreadsPerBlock;
  uint32_t nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  sutra::SlopesIndex si{blocks, slope_order};

  threads /= nelem_thread;
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  if (threads <= 16)
    centroids<  16><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, threshold, npix, si, size,
                                             T(scale), T(offset), nelem_thread);
  else if (threads <= 36)
    centroids<  36><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, threshold, npix, si, size,
                                             T(scale), T(offset), nelem_thread);

  else if (threads <= 64)
    centroids<  64><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, threshold, npix, si, size,
                                             T(scale), T(offset), nelem_thread);
  else if (threads <= 100)
    centroids< 100><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, threshold, npix, si, size,
                                             T(scale), T(offset), nelem_thread);
  else if (threads <= 144)
    centroids< 144><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, threshold, npix, si, size,
                                             T(scale), T(offset), nelem_thread);
  else if (threads <= 256)
    centroids< 256><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, threshold, npix, si, size,
                                             T(scale), T(offset), nelem_thread);
  else if (threads <= 512)
    centroids< 512><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, threshold, npix, si, size,
                                             T(scale), T(offset), nelem_thread);
  else if (threads <= 1024)
    centroids<1024><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, threshold, npix, si, size,
                                             T(scale), T(offset), nelem_thread);
  else
    printf("SH way too big !!!\n");

  carma_check_msg("centroids_kernel<<<>>> execution failed\n");
}

template void get_centroids<float>(int32_t size, int32_t threads, int32_t blocks, int32_t npix,
                                   float *d_img, float *d_centroids, float *ref,
                                   int32_t *validx, int32_t *validy, float *intensities,
                                   float threshold, float scale, float offset,
                                   SlopeOrder slope_order,
                                   CarmaDevice *device, cudaStream_t stream);

template void get_centroids<double>(int32_t size, int32_t threads, int32_t blocks, int32_t npix,
                                    float *d_img, double *d_centroids,
                                    double *ref, int32_t *validx, int32_t *validy,
                                    float *intensities, float threshold,
                                    float scale, float offset,
                                    SlopeOrder slope_order,
                                    CarmaDevice *device, cudaStream_t stream);