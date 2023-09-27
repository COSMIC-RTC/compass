// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_pbcog.cu
//! \ingroup   libsutra
//! \class     sutra_centroider_pbcog
//! \brief     this class provides the centroider_pbcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_centroider_bpcog.h>
#include <sutra_centroider_utils.cuh>
#include <carma_utils.cuh>

template <int BLOCK_THREADS, typename T>
__launch_bounds__(BLOCK_THREADS) __global__
    void centroids(float *d_img, T *d_centroids, T *ref, int *validx,
                   int *validy, float *d_intensities, int nbpix,
                   unsigned int npix, sutra::SlopesIndex si, unsigned int size, T scale, T offset,
                   unsigned int nelem_thread) {
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

  unsigned int tid = threadIdx.x;
  unsigned int xvalid = validx[blockIdx.x];
  unsigned int yvalid = validy[blockIdx.x];
  unsigned int x = tid % npix;
  unsigned int y = tid / npix;
  int idim = (x + xvalid) + (y + yvalid) * size;

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
void get_centroids(int size, int threads, int blocks, int npix, float *d_img,
                   T *d_centroids, T *ref, int *validx, int *validy,
                   float *intensities, int nbpix, float scale, float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device, cudaStream_t stream) {
  int maxThreads = device->get_properties().maxThreadsPerBlock;
  unsigned int nelem_thread = 1;
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

template void get_centroids<float>(int size, int threads, int blocks, int npix,
                                   float *d_img, float *d_centroids, float *ref,
                                   int *validx, int *validy, float *intensities,
                                   int nbpix, float scale, float offset,
                                   SlopeOrder slope_order,
                                   CarmaDevice *device, cudaStream_t stream);

template void get_centroids<double>(int size, int threads, int blocks, int npix,
                                    float *d_img, double *d_centroids,
                                    double *ref, int *validx, int *validy,
                                    float *intensities, int nbpix, float scale,
                                    float offset, SlopeOrder slope_order,
                                    CarmaDevice *device, cudaStream_t stream);
#ifdef CAN_DO_HALF
template void get_centroids<half>(int size, int threads, int blocks, int npix,
                                  float *d_img, half *d_centroids, half *ref,
                                  int *validx, int *validy, float *intensities,
                                  int nbpix, float scale, float offset,
                                  SlopeOrder slope_order,
                                  CarmaDevice *device, cudaStream_t stream);
#endif