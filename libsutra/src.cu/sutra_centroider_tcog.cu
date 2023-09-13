// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_tcog.cu
//! \ingroup   libsutra
//! \class     SutraCentroiderTcog
//! \brief     this class provides the centroider_tcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

#include <sutra_centroider_tcog.h>
#include <sutra_centroider_utils.cuh>
#include <carma_utils.cuh>

template <int nb_threads, typename T>
__global__ void centroids(float *d_img, T *d_centroids, T *ref, int *validx,
                          int *validy, float *d_intensities, float threshold,
                          unsigned int npix, sutra::SlopesIndex si, unsigned int size, T scale,
                          T offset, unsigned int nelem_thread) {
  if (blockDim.x > nb_threads) {
    if (threadIdx.x == 0) printf("Wrong size argument\n");
    return;
  }
  // Specialize BlockReduce for a 1D block of 128 threads on type int
  typedef cub::BlockReduce<float, nb_threads> BlockReduce;
  // Allocate shared memory for BlockReduce
  __shared__ typename BlockReduce::TempStorage temp_storage;

  float idata = 0;
  float xdata = 0;
  float ydata = 0;
  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int xvalid = validx[blockIdx.x];
  unsigned int yvalid = validy[blockIdx.x];
  unsigned int x, y;
  int idim;

  for (int cc = 0; cc < nelem_thread; cc++) {
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
void get_centroids(int size, int threads, int blocks, int npix, float *d_img,
                   T *d_centroids, T *ref, int *validx, int *validy,
                   float *intensities, float threshold, float scale,
                   float offset,
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

template void get_centroids<float>(int size, int threads, int blocks, int npix,
                                   float *d_img, float *d_centroids, float *ref,
                                   int *validx, int *validy, float *intensities,
                                   float threshold, float scale, float offset,
                                   SlopeOrder slope_order,
                                   CarmaDevice *device, cudaStream_t stream);

template void get_centroids<double>(int size, int threads, int blocks, int npix,
                                    float *d_img, double *d_centroids,
                                    double *ref, int *validx, int *validy,
                                    float *intensities, float threshold,
                                    float scale, float offset,
                                    SlopeOrder slope_order,
                                    CarmaDevice *device, cudaStream_t stream);
#ifdef CAN_DO_HALF
template void get_centroids<half>(int size, int threads, int blocks, int npix,
                                  float *d_img, half *d_centroids, half *ref,
                                  int *validx, int *validy, float *intensities,
                                  float threshold, float scale, float offset,
                                  SlopeOrder slope_order,
                                  CarmaDevice *device, cudaStream_t stream);
#endif