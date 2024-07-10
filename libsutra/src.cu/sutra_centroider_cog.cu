// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_cog.cu
//! \ingroup   libsutra
//! \class     SutraCentroiderCog
//! \brief     this class provides the centroider_cog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24


#include <sutra_centroider_cog.hpp>
#include <sutra_centroider_utils.cuh>
#include <carma_utils.cuh>

template <int32_t nb_threads, typename T>
__global__ void centroids(float *d_img, T *d_centroids, T *ref, int32_t *validx,
                          int32_t *validy, float *d_intensities,
                          uint32_t npix, sutra::SlopesIndex si,
                          uint32_t size, T scale, T offset,
                          uint32_t nelem_thread) {
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
      idata += d_img[idim];
      xdata += d_img[idim] * x;
      ydata += d_img[idim] * y;
    }
  }

  // sdata[tid] = (i < N) ? g_idata[i] * x : 0;
  __syncthreads();

  float intensity = BlockReduce(temp_storage).Sum(idata, npix * npix);
  __syncthreads();
  float slopex    = BlockReduce(temp_storage).Sum(xdata, npix * npix);
  __syncthreads();
  float slopey    = BlockReduce(temp_storage).Sum(ydata, npix * npix);

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
                   float *intensities, float scale, float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device, cudaStream_t stream) {
  int32_t maxThreads = device->get_properties().maxThreadsPerBlock;
  uint32_t nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  threads /= nelem_thread;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  sutra::SlopesIndex si{blocks, slope_order};

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  if (threads <= 16)
    centroids<  16><<<dimGrid, dimBlock, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                                         intensities, npix, si, size,
                                                         T(scale), T(offset), nelem_thread);
  else if (threads <= 36)
    centroids<  36><<<dimGrid, dimBlock, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                                         intensities, npix, si, size,
                                                         T(scale), T(offset), nelem_thread);
  else if (threads <= 64)
    centroids<  64><<<dimGrid, dimBlock, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                                         intensities, npix, si, size,
                                                         T(scale), T(offset), nelem_thread);
  else if (threads <= 100)
    centroids< 100><<<dimGrid, dimBlock, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                                         intensities, npix, si, size,
                                                         T(scale), T(offset), nelem_thread);
  else if (threads <= 144)
    centroids< 144><<<dimGrid, dimBlock, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                                         intensities, npix, si, size,
                                                         T(scale), T(offset), nelem_thread);
  // else if (threads <= 196)
  // centroids< 196><<<dimGrid, 196, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
  //                                     intensities, npix, size, T(scale),
  //                                     T(offset), nelem_thread);
  else if (threads <= 256)
    centroids< 256><<<dimGrid, dimBlock, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                                         intensities, npix, si, size,
                                                         T(scale), T(offset), nelem_thread);
  else if (threads <= 512)
    centroids< 512><<<dimGrid, dimBlock, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                                         intensities, npix, si, size,
                                                         T(scale), T(offset), nelem_thread);
  else if (threads <= 1024)
    centroids<1024><<<dimGrid, dimBlock, 0, stream>>>(d_img, d_centroids, ref, validx, validy,
                                                         intensities, npix, si, size,
                                                         T(scale), T(offset), nelem_thread);
  else
    printf("SH way too big !!!\n");

  carma_check_msg("centroids_kernel<<<>>> execution failed\n");

  //   centroidy<T><<<dimGrid, dimBlock, smemSize>>>(
  //       d_idata, &(d_odata[blocks]), alpha, n, size, scale, offset,
  //       nelem_thread);

  //   carma_check_msg("centroidy_kernel<<<>>> execution failed\n");
}

template void get_centroids<float>(int32_t size, int32_t threads, int32_t blocks, int32_t npix,
                                   float *d_img, float *d_centroids, float *ref,
                                   int32_t *validx, int32_t *validy, float *intensities,
                                   float scale, float offset,
                                   SlopeOrder slope_order,
                                   CarmaDevice *device, cudaStream_t stream);

template void get_centroids<double>(int32_t size, int32_t threads, int32_t blocks, int32_t npix,
                                    float *d_img, double *d_centroids, double *ref,
                                    int32_t *validx, int32_t *validy, float *intensities,
                                    float scale, float offset,
                                    SlopeOrder slope_order,
                                    CarmaDevice *device, cudaStream_t stream);

#ifdef CAN_DO_HALF
template void get_centroids<half>(int32_t size, int32_t threads, int32_t blocks, int32_t npix,
                                  float *d_img, half *d_centroids, half *ref,
                                  int32_t *validx, int32_t *validy, float *intensities,
                                  float scale, float offset,
                                  SlopeOrder slope_order,
                                  CarmaDevice *device, cudaStream_t stream);
#endif

// template <class T>
// __global__ void centroidx(T *g_idata, T *g_odata, T *alpha, uint32_t n,
//                           uint32_t N, float scale, float offset,
//                           uint32_t nelem_thread) {
//   T *sdata = SharedMemory<T>();

//   // load shared mem
//   uint32_t tid = threadIdx.x;
//   // uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//   // uint32_t x = (tid % n) + 1;
//   uint32_t x;
//   int32_t idim;
//   sdata[tid] = 0;
//   for (int32_t cc = 0; cc < nelem_thread; cc++) {
//     x = ((tid * nelem_thread + cc) % n);
//     idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) *
//     blockIdx.x; if (idim < N)
//       sdata[tid] += g_idata[idim] * x;
//     else
//       sdata[tid] += 0;
//   }

//   // sdata[tid] = (i < N) ? g_idata[i] * x : 0;

//   __syncthreads();

//   reduce_krnl(sdata, blockDim.x, tid);

//   // write result for this block to global mem
//   if (tid == 0)
//     g_odata[blockIdx.x] =
//         ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6)) - offset) * scale;
// }

// template <class T>
// __global__ void centroidy(T *g_idata, T *g_odata, T *alpha, uint32_t n,
//                           uint32_t N, float scale, float offset,
//                           uint32_t nelem_thread) {
//   T *sdata = SharedMemory<T>();

//   // load shared mem
//   uint32_t tid = threadIdx.x;
//   // uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//   // uint32_t y = (tid / n) + 1;
//   uint32_t y;
//   int32_t idim;
//   sdata[tid] = 0;
//   for (int32_t cc = 0; cc < nelem_thread; cc++) {
//     y = ((tid * nelem_thread + cc) / n);
//     idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) *
//     blockIdx.x; if (idim < N)
//       sdata[tid] += g_idata[idim] * y;
//     else
//       sdata[tid] += 0;
//   }

//   // sdata[tid] = (i < N) ? g_idata[i] * y : 0;

//   __syncthreads();

//   reduce_krnl(sdata, blockDim.x, tid);

//   // write result for this block to global mem
//   if (tid == 0)
//     g_odata[blockIdx.x] =
//         ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6)) - offset) * scale;
// }

// template <class T>
// __global__ void centroidx_async(T *g_idata, T *g_odata, T *alpha,
//                                 uint32_t n, uint32_t N, float scale,
//                                 float offset, int32_t stream_offset) {
//   T *sdata = SharedMemory<T>();

//   // load shared mem
//   uint32_t tid = threadIdx.x;
//   uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//   uint32_t x = (tid % n) + 1;
//   i += stream_offset * blockDim.x;

//   sdata[tid] = (i < N) ? g_idata[i] * x : 0;

//   __syncthreads();

//   reduce_krnl(sdata, blockDim.x, tid);

//   // write result for this block to global mem
//   if (tid == 0)
//     g_odata[blockIdx.x + stream_offset] =
//         ((sdata[0] * 1.0 / (alpha[blockIdx.x + stream_offset] + 1.e-6)) -
//          offset) *
//         scale;
// }

// template <class T>
// __global__ void centroidy_async(T *g_idata, T *g_odata, T *alpha,
//                                 uint32_t n, uint32_t N, float scale,
//                                 float offset, int32_t stream_offset) {
//   T *sdata = SharedMemory<T>();

//   // load shared mem
//   uint32_t tid = threadIdx.x;
//   uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//   uint32_t y = (tid / n) + 1;
//   i += stream_offset * blockDim.x;

//   sdata[tid] = (i < N) ? g_idata[i] * y : 0;

//   __syncthreads();

//   reduce_krnl(sdata, blockDim.x, tid);

//   // write result for this block to global mem
//   if (tid == 0)
//     g_odata[blockIdx.x + stream_offset] =
//         ((sdata[0] * 1.0 / (alpha[blockIdx.x + stream_offset] + 1.e-6)) -
//          offset) *
//         scale;
// }

// template <class T>
// void get_centroids_async(int32_t threads, int32_t blocks, int32_t n, CarmaStreams
// *streams,
//                          T *d_idata, T *d_odata, T *alpha, float scale, float
//                          offset)
//                          {
//   int32_t nstreams = streams->get_nb_streams();
//   int32_t nbelem = threads * blocks;

//   dim3 dimBlock(threads);
//   dim3 dimGrid(blocks / nstreams);

//   // when there is only one warp per block, we need to allocate two warps
//   // worth of shared memory so that we don't index shared memory out of
//   bounds int32_t smemSize =
//       (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
//   for (int32_t i = 0; i < nstreams; i++) {
//     centroidx_async<T><<<dimGrid, dimBlock, smemSize,
//     streams->get_stream(i)>>>(
//         d_idata, d_odata, alpha, n, nbelem, scale, offset,
//         i * blocks / nstreams);

//     carma_check_msg("centroidx_kernel<<<>>> execution failed\n");

//     centroidy_async<T><<<dimGrid, dimBlock, smemSize,
//     streams->get_stream(i)>>>(
//         d_idata, &(d_odata[blocks]), alpha, n, nbelem, scale, offset,
//         i * blocks / nstreams);

//     carma_check_msg("centroidy_kernel<<<>>> execution failed\n");
//   }
// }

// template void get_centroids_async<float>(int32_t threads, int32_t blocks, int32_t n,
//                                          CarmaStreams *streams, float
//                                          *d_idata, float *d_odata, float
//                                          *alpha, float scale, float offset);
// template void get_centroids_async<double>(int32_t threads, int32_t blocks, int32_t n,
//                                           CarmaStreams *streams,
//                                           double *d_idata, double *d_odata,
//                                           double *alpha, double scale,
//                                           double offset);
