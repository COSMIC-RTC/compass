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
//! \version   5.4.2
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
                   CarmaDevice *device) {
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
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 36)
    centroids<  36><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 64)
    centroids<  64><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 100)
    centroids< 100><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 144)
    centroids< 144><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 256)
    centroids< 256><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 512)
    centroids< 512><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
                                          intensities, nbpix, npix, si, size,
                                          T(scale), T(offset), nelem_thread);
  else if (threads <= 1024)
    centroids<1024><<<dimGrid, threads>>>(d_img, d_centroids, ref, validx, validy,
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
                                   CarmaDevice *device);

template void get_centroids<double>(int size, int threads, int blocks, int npix,
                                    float *d_img, double *d_centroids,
                                    double *ref, int *validx, int *validy,
                                    float *intensities, int nbpix, float scale,
                                    float offset, SlopeOrder slope_order,
                                    CarmaDevice *device);
#ifdef CAN_DO_HALF
template void get_centroids<half>(int size, int threads, int blocks, int npix,
                                  float *d_img, half *d_centroids, half *ref,
                                  int *validx, int *validy, float *intensities,
                                  int nbpix, float scale, float offset,
                                  SlopeOrder slope_order,
                                  CarmaDevice *device);
#endif

template <class T>
__device__ inline void sortmax_krnl(T *sdata, unsigned int *values, int size,
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
__global__ void sortmax(T *g_idata, T *g_odata, unsigned int *values, int nmax,
                        int Npix, int size, int nelem_thread) {
  extern __shared__ uint svalues[];
  T *sdata = (T *)&svalues[Npix];
  /*
    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    svalues[tid] = tid;
    sdata[tid] = g_idata[i];
  */
  unsigned int tid = threadIdx.x;
  // unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  // unsigned int y = (tid / n) + 1;
  int idim;
  int sdim;

  for (int cc = 0; cc < nelem_thread; cc++) {
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    sdim = tid * nelem_thread + cc;
    if (idim < size) {
      sdata[sdim] = g_idata[idim];
      svalues[sdim] = sdim;
    }
  }

  __syncthreads();

  for (int cc = 0; cc < nmax; cc++) {
    for (int cpt = 0; cpt < nelem_thread; cpt++) {
      sdim = tid * nelem_thread + cpt;
      if (sdim >= cc)
        sortmax_krnl(&(sdata[cc]), &(svalues[cc]), Npix - cc, sdim - cc);
      __syncthreads();
    }
  }
  for (int cpt = 0; cpt < nelem_thread; cpt++) {
    sdim = tid * nelem_thread + cpt;
    if (sdim < nmax) {
      g_odata[nmax * blockIdx.x + sdim] = sdata[sdim];  // - sdata[nmax - 1];
      values[nmax * blockIdx.x + sdim] = svalues[sdim];
    }
  }
  /*
   __syncthreads();
   if ((blockIdx.x == 0) && (tid < nmax))
   printf("tid %d sdata %f \n",tid,g_odata[tid]);
   */
}

template <class T>
void subap_sortmax(int threads, int blocks, T *d_idata, T *d_odata,
                   unsigned int *values, int nmax, CarmaDevice *device) {
  int maxThreads = device->get_properties().maxThreadsPerBlock;
  unsigned int nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  dim3 dimBlock(threads / nelem_thread, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  size_t smemSize = threads * (sizeof(T) + sizeof(uint));
  sortmax<T><<<dimGrid, dimBlock, smemSize>>>(
      d_idata, d_odata, values, nmax, threads, threads * blocks, nelem_thread);

  carma_check_msg("sortmax_kernel<<<>>> execution failed\n");
}
template void subap_sortmax<float>(int threads, int blocks, float *d_idata,
                                   float *d_odata, unsigned int *values,
                                   int nmax, CarmaDevice *device);
template void subap_sortmax<double>(int threads, int blocks, double *d_idata,
                                    double *d_odata, unsigned int *values,
                                    int nmax, CarmaDevice *device);

template <class T>
__global__ void centroid_bpix(int nsub, int n, T *g_idata, unsigned int *values,
                              T *g_odata, float scale, float offset) {
  extern __shared__ uint svalues[];
  T *sdata = (T *)&svalues[blockDim.x];
  T intensities;
  // T minimum;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  svalues[tid] = values[i];
  sdata[tid] = g_idata[i] - g_idata[blockIdx.x * blockDim.x + blockDim.x - 1];

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  __syncthreads();
  // get the sum per subap
  if (tid == 0) intensities = (abs(sdata[tid]) > 1.e-6 ? sdata[tid] : 0.0f);

  __syncthreads();

  // Reload sdata
  sdata[tid] = g_idata[i] - g_idata[blockIdx.x * blockDim.x + blockDim.x - 1];

  __syncthreads();

  // compute the centroid on the first part of the array
  sdata[tid] *= ((svalues[tid] % n));
  // x centroid
  __syncthreads();
  reduce_krnl(sdata, blockDim.x, tid);
  //__syncthreads();
  if (tid == 0)
    g_odata[blockIdx.x] =
        (intensities != 0.0f ? ((sdata[tid] / intensities) - offset) * scale
                             : 0.0f);
  __syncthreads();
  sdata[tid] = g_idata[i] - g_idata[blockIdx.x * blockDim.x + blockDim.x - 1];

  __syncthreads();

  // compute the centroid on the first part of the array
  sdata[tid] *= (svalues[tid] / n);
  // y centroid
  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);
  //__syncthreads();
  if (tid == 0)
    g_odata[blockIdx.x + nsub] =
        (intensities != 0.0f ? ((sdata[tid] / intensities) - offset) * scale
                             : 0.0f);
}

template <class T>
void subap_bpcentro(int threads, int blocks, int npix, T *d_idata,
                    unsigned int *values, T *d_odata, float scale,
                    float offset) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = threads * (sizeof(T) + sizeof(uint));

  centroid_bpix<T><<<dimGrid, dimBlock, smemSize>>>(
      blocks, npix, d_idata, values, d_odata, scale, offset);

  carma_check_msg("centroid_bpix<<<>>> execution failed\n");
}
template void subap_bpcentro<float>(int threads, int blocks, int npix,
                                    float *d_idata, unsigned int *values,
                                    float *d_odata, float scale, float offset);
template void subap_bpcentro<double>(int threads, int blocks, int npix,
                                     double *d_idata, unsigned int *values,
                                     double *d_odata, float scale,
                                     float offset);
