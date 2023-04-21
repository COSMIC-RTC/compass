// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_wcog.cu
//! \ingroup   libsutra
//! \class     SutraCentroiderWcog
//! \brief     this class provides the centroider_wcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24

#include <sutra_centroider_wcog.h>
#include <sutra_centroider_utils.cuh>
#include <carma_utils.cuh>

template <class T>
__global__ void fillweights_krnl(T *d_out, T *weights, int Npix, int N) {
  int nim, idx;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / Npix;
    idx = tid - nim * Npix;
    d_out[tid] = weights[idx];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int fill_weights(T *d_out, T *d_in, int npix, int N, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fillweights_krnl<<<grid, threads>>>(d_out, d_in, npix * npix, N);
  carma_check_msg("<<<fillweights_krnl>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int fill_weights<float>(float *d_out, float *d_in, int npix, int N,
                                CarmaDevice *device);

template int fill_weights<double>(double *d_out, double *d_in, int npix, int N,
                                 CarmaDevice *device);

template <int nb_threads, typename T>
__global__ void centroids(float *d_img, T *d_centroids, T *ref, int *validx,
                          int *validy, float *d_intensities, float *weights, float threshold, 
                          unsigned int npix, sutra::SlopesIndex si, unsigned int size, float scale,
                          float offset, unsigned int nelem_thread) {
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
  int idim, wdim;

  for (int cc = 0; cc < nelem_thread; cc++) {
    x = ((tid * nelem_thread + cc) % npix);
    y = ((tid * nelem_thread + cc) / npix);
    // idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) *
    // blockIdx.x;
    idim = (x + xvalid) + (y + yvalid) * size;
    wdim = x + y * npix + blockIdx.x * npix * npix;
    if (idim < size * size) {
      float data_thresh =
          (d_img[idim] > threshold) ? d_img[idim] - threshold : 1.e-6;
      idata += data_thresh * weights[wdim];
      xdata += data_thresh * x * weights[wdim];
      ydata += data_thresh * y * weights[wdim];
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
    d_centroids[si.x(blockIdx.x)] = (T(slopex * 1.0 / (intensity + 1.e-6)) - offset) * scale - ref[si.x(blockIdx.x)];
    d_centroids[si.y(blockIdx.x)] = (T(slopey * 1.0 / (intensity + 1.e-6)) - offset) * scale - ref[si.y(blockIdx.x)];
    d_intensities[blockIdx.x] = intensity;
  }
}

template <class T>
void get_centroids(int size, int threads, int blocks, int npix, float *d_img,
                   T *d_centroids, T *ref, int *validx, int *validy,
                   float *intensities, float *weights, float threshold, float scale,
                   float offset,
                   SlopeOrder slope_order, CarmaDevice *device) {
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
    centroids<  16><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, weights, threshold, npix, si, size, scale,
                                             offset, nelem_thread);
  else if (threads <= 36)
    centroids<  36><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, weights, threshold, npix, si, size, scale,
                                             offset, nelem_thread);
  else if (threads <= 64)
    centroids<  64><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, weights, threshold, npix, si, size, scale,
                                             offset, nelem_thread);
  else if (threads <= 100)
    centroids< 100><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, weights, threshold, npix, si, size, scale,
                                             offset, nelem_thread);
  else if (threads <= 144)
    centroids< 144><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, weights, threshold, npix, si, size, scale,
                                             offset, nelem_thread);
  else if (threads <= 256)
    centroids< 256><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, weights, threshold, npix, si, size, scale,
                                             offset, nelem_thread);
  else if (threads <= 512)
    centroids< 512><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                             intensities, weights, threshold, npix, si, size, scale,
                                             offset, nelem_thread);
  else if (threads <= 1024)
    centroids<1024><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx,
                                             validy, intensities, weights, threshold, npix, si,
                                             size, scale, offset, nelem_thread);
  else
    printf("SH way too big !!!\n");

  carma_check_msg("centroids_kernel<<<>>> execution failed\n");
}

template void get_centroids<float>(int size, int threads, int blocks, int npix,
                                   float *d_img, float *d_centroids, float *ref,
                                   int *validx, int *validy, float *intensities,
                                   float *weights, float threshold, float scale, float offset,
                                   SlopeOrder slope_order,
                                   CarmaDevice *device);

template void get_centroids<double>(int size, int threads, int blocks, int npix,
                                    float *d_img, double *d_centroids,
                                    double *ref, int *validx, int *validy,
                                    float *intensities, float *weights, float threshold,
                                    float scale, float offset,
                                    SlopeOrder slope_order,
                                    CarmaDevice *device);

// template <class T>
// __global__ void centroidx(T *g_idata, T *g_odata, T *alpha, T *weights,
//                           unsigned int n, unsigned int N, float scale, float
//                           offset, unsigned int nelem_thread) {
//   T *sdata = SharedMemory<T>();

//   // load shared mem
//   unsigned int tid = threadIdx.x;
//   // unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
//   // unsigned int x = (tid % n) + 1;
//   unsigned int x;
//   int idim;
//   sdata[tid] = 0.0f;
//   for (int cc = 0; cc < nelem_thread; cc++) {
//     x = ((tid * nelem_thread + cc) % n);
//     idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) *
//     blockIdx.x; if (idim < N)
//       sdata[tid] += g_idata[idim] * x * weights[idim];
//     else
//       sdata[tid] += 0.0f;
//   }
//   __syncthreads();

//   // sdata[tid] = (i < N) ? g_idata[i] * x * weights[i] : 0;

//   __syncthreads();

//   reduce_krnl(sdata, blockDim.x, tid);
//   // if(tid == 0)
//   //	printf("blockIdx %d sdata %f \n",blockIdx.x,sdata[tid]);
//   // write result for this block to global mem
//   if (tid == 0)
//     g_odata[blockIdx.x] =
//         ((sdata[tid] * 1.0 / (alpha[blockIdx.x] + 1.e-6)) - offset) * scale;
// }

// template <class T>
// __global__ void centroidy(T *g_idata, T *g_odata, T *alpha, T *weights,
//                           unsigned int n, unsigned int N, float scale, float
//                           offset, unsigned int nelem_thread) {
//   T *sdata = SharedMemory<T>();

//   // load shared mem
//   unsigned int tid = threadIdx.x;
//   // unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
//   // unsigned int y = (tid / n) + 1;
//   unsigned int y;
//   int idim;
//   sdata[tid] = 0.0f;
//   for (int cc = 0; cc < nelem_thread; cc++) {
//     y = ((tid * nelem_thread + cc) / n);
//     idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) *
//     blockIdx.x; if (idim < N)
//       sdata[tid] += g_idata[idim] * y * weights[idim];
//     else
//       sdata[tid] += 0.0f;
//   }

//   // sdata[tid] = (i < N) ? g_idata[i] * y * weights[i] : 0;

//   __syncthreads();

//   reduce_krnl(sdata, blockDim.x, tid);

//   // write result for this block to global mem
//   if (tid == 0)
//     g_odata[blockIdx.x] =
//         ((sdata[tid] * 1.0 / (alpha[blockIdx.x] + 1.e-6)) - offset) * scale;
// }

// template <class T>
// void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
//                    T *d_odata, T *alpha, T *weights, float scale, float
//                    offset, CarmaDevice *device) {
//   int maxThreads = device->get_properties().maxThreadsPerBlock;
//   unsigned int nelem_thread = 1;
//   while ((threads / nelem_thread > maxThreads) ||
//          (threads % nelem_thread != 0)) {
//     nelem_thread++;
//   }

//   threads /= nelem_thread;
//   dim3 dimBlock(threads, 1, 1);
//   dim3 dimGrid(blocks, 1, 1);

//   // when there is only one warp per block, we need to allocate two warps
//   // worth of shared memory so that we don't index shared memory out of
//   bounds int smemSize =
//       (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
//   centroidx<T><<<dimGrid, dimBlock, smemSize>>>(
//       d_idata, d_odata, alpha, weights, n, size, scale, offset,
//       nelem_thread);

//   carma_check_msg("centroidx_kernel<<<>>> execution failed\n");

//   centroidy<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, &(d_odata[blocks]),
//                                                 alpha, weights, n, size,
//                                                 scale, offset, nelem_thread);

//   carma_check_msg("centroidy_kernel<<<>>> execution failed\n");
// }

// template void get_centroids<float>(int size, int threads, int blocks, int n,
//                                    float *d_idata, float *d_odata, float
//                                    *alpha, float *weights, float scale, float
//                                    offset, CarmaDevice *device);

// template void get_centroids<double>(int size, int threads, int blocks, int n,
//                                     double *d_idata, double *d_odata,
//                                     double *alpha, double *weights,
//                                     double scale, double offset,
//                                     CarmaDevice *device);
