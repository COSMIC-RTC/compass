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

//! \file      sutra_centroider_pbcog.cu
//! \ingroup   libsutra
//! \class     sutra_centroider_pbcog
//! \brief     this class provides the centroider_pbcog features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_centroider_bpcog.hpp>
#include <sutra_centroider_utils.cuh>
#include <carma_utils.cuh>

template <int32_t BLOCK_THREADS, typename T>
__launch_bounds__(BLOCK_THREADS) __global__
    void centroids(float *d_img, T *d_centroids, T *ref, int32_t *validx,
                   int32_t *validy, float *d_intensities, int32_t nbpix,
                   uint32_t npix, sutra::SlopesIndex si, uint32_t size, T scale, T offset,
                   uint32_t nelem_thread) {
  // Specialize BlockRadixSort for a 1D block of BLOCK_THREADS threads owning 1
  // item each
  typedef cub::BlockRadixSort<float, BLOCK_THREADS, 1> BlockRadixSortT;
  typedef cub::BlockReduce<float, BLOCK_THREADS> BlockReduce;

  // Allocate shared memory for BlockRadixSort
  __shared__ typename BlockRadixSortT::TempStorage temp_storageSort;
  __shared__ typename BlockReduce::TempStorage temp_storageSum;
  __shared__ float threshold;

  float idata = 0;
  float xdata = 0;
  float ydata = 0;

  uint32_t tid = threadIdx.x;
  uint32_t xvalid = validx[blockIdx.x];
  uint32_t yvalid = validy[blockIdx.x];
  uint32_t x = tid % npix;
  uint32_t y = tid / npix;
  int32_t idim = (x + xvalid) + (y + yvalid) * size;

  float items[1];
  items[0] = ((idim < size * size) && (tid < npix * npix)) ? d_img[idim] : 0.f;

  __syncthreads();
  BlockRadixSortT(temp_storageSort).SortDescending(items);

  if (tid == nbpix) threshold = items[0];

  __syncthreads();
  if ((idim < size * size) && (tid < npix * npix)) {
    float data_thresh =
        (d_img[idim] > threshold) ? d_img[idim] - threshold : 0.f;
    idata += data_thresh;
    xdata += data_thresh * x;
    ydata += data_thresh * y;
    d_img[idim] = data_thresh;
  }

  __syncthreads();
  float intensity = BlockReduce(temp_storageSum).Sum(idata, npix * npix);
  __syncthreads();
  float slopex = BlockReduce(temp_storageSum).Sum(xdata, npix * npix);
  __syncthreads();
  float slopey = BlockReduce(temp_storageSum).Sum(ydata, npix * npix);

  if (tid == 0) {
    d_centroids[si.x(blockIdx.x)] = (T(slopex / (intensity + 1.e-6)) - offset) * scale - ref[si.x(blockIdx.x)];
    d_centroids[si.y(blockIdx.x)] = (T(slopey / (intensity + 1.e-6)) - offset) * scale - ref[si.y(blockIdx.x)];
    d_intensities[blockIdx.x] = intensity;
  }
}

template <class T>
void get_centroids(int32_t size, int32_t threads, int32_t blocks, int32_t npix, float *d_img,
                   T *d_centroids, T *ref, int32_t *validx, int32_t *validy,
                   float *intensities, int32_t nbpix, float scale, float offset,
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
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 36)
    centroids<  36><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 64)
    centroids<  64><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 100)
    centroids< 100><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 144)
    centroids< 144><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 256)
    centroids< 256><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 512)
    centroids< 512><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 1024)
    centroids<1024><<<dimGrid, threads, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else
    printf("SH way too big !!!\n");

  carma_check_msg("centroids_kernel<<<>>> execution failed\n");
}

template void get_centroids<float>(int32_t size, int32_t threads, int32_t blocks, int32_t npix,
                                   float *d_img, float *d_centroids, float *ref,
                                   int32_t *validx, int32_t *validy, float *intensities,
                                   int32_t nbpix, float scale, float offset,
                                   SlopeOrder slope_order,
                                   CarmaDevice *device, cudaStream_t stream);

template void get_centroids<double>(int32_t size, int32_t threads, int32_t blocks, int32_t npix,
                                    float *d_img, double *d_centroids,
                                    double *ref, int32_t *validx, int32_t *validy,
                                    float *intensities, int32_t nbpix, float scale,
                                    float offset, SlopeOrder slope_order,
                                    CarmaDevice *device, cudaStream_t stream);