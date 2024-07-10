// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_dm.cu
//! \ingroup   libsutra
//! \class     SutraDm
//! \brief     this class provides the dm features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_dm.hpp>
#include "carma_utils.cuh"

#include <cuda.h>

#include <stdio.h>

__device__ inline struct doubleint getInterval(const int32_t pos,
                                               const int32_t *iStart_t) {
  const int32_t start = iStart_t[pos];
  return {start, iStart_t[pos + 1] - start};
}

template <class T>
__device__ inline T get_data(const int32_t pos, const T *cmdVector,
                            const T *influData, const int32_t *iPos) {
  return cmdVector[iPos[pos]] * influData[pos];
}

template <class T>
__global__ void compShape2(T *outData, const T *cmdVector, const T *influData,
                           const int32_t *iPos, const int32_t *iStart_t, const int32_t N) {
  const int32_t id = threadIdx.x + blockIdx.x * blockDim.x;

  if (id < N) {
    const struct doubleint interval = getInterval(id, iStart_t);

    T sum = 0;

    for (int32_t pos = interval.start; pos < interval.start + interval.nbInflu;
         ++pos)
      sum += get_data(pos, cmdVector, influData, iPos);

    outData[id] = sum;
  }
}

template <class T>
void comp_dmshape2(T *outData, const T *cmdVector, const T *influData,
                   const int32_t *iStart_t, const int32_t *iPos, const int32_t roiLength,
                   const dim3 threads, const dim3 blocks, const int32_t shared) {
  compShape2<T><<<blocks, threads>>>(outData, cmdVector, influData, iPos,
                                     iStart_t, roiLength);

  carma_check_msg("comp_dmshape2<<<>>> execution failed\n");
}

template void comp_dmshape2<float>(float *outData, const float *cmdVector,
                                   const float *influData, const int32_t *iStart_t,
                                   const int32_t *iPos, const int32_t roiLength,
                                   const dim3 threads, const dim3 blocks,
                                   const int32_t shared);

template void comp_dmshape2<double>(double *outData, const double *cmdVector,
                                    const double *influData,
                                    const int32_t *iStart_t, const int32_t *iPos,
                                    const int32_t roiLength, const dim3 threads,
                                    const dim3 blocks, const int32_t shared);

template <class T>
__global__ void dmshape_krnl(T *g_idata, T *g_odata, int32_t *pos, int32_t *istart,
                             int32_t *npts, T *comm, uint32_t n, int32_t N) {
  // T *sdata = SharedMemory<T>();
  // load shared mem
  // uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int32_t local_istart = istart[i];
    int32_t local_npts = npts[i];

    // sdata[tid] = 0;
    T value = 0;

    if (local_npts > 0) {
      for (int32_t cc = 0; cc < local_npts; cc++) {
        int32_t lpos = pos[local_istart + cc];
        int32_t ninflu = lpos / n;
        // sdata[tid] += comm[ninflu] * g_idata[lpos];
        value += comm[ninflu] * g_idata[lpos];
      }
    }
    g_odata[i] = value;
    i += blockDim.x * gridDim.x;
  }
  /*
  __syncthreads();

  if (i < N) {
  // write result for this block to global mem
  g_odata[i] = sdata[tid];
  }
   */
}

template <class T>
void comp_dmshape(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t *pos,
                  int32_t *istart, int32_t *npts, T *comm, uint32_t n, int32_t N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  // int32_t smemSize =
  //    (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

  dmshape_krnl<T><<<dimGrid, dimBlock /*, smemSize*/>>>(
      d_idata, d_odata, pos, istart, npts, comm, n, N);

  carma_check_msg("dmshape_kernel<<<>>> execution failed\n");
}

template void comp_dmshape<float>(int32_t threads, int32_t blocks, float *d_idata,
                                  float *d_odata, int32_t *pos, int32_t *istart,
                                  int32_t *npts, float *comm, uint32_t n,
                                  int32_t N);

template void comp_dmshape<double>(int32_t threads, int32_t blocks, double *d_idata,
                                   double *d_odata, int32_t *pos, int32_t *istart,
                                   int32_t *npts, double *comm, uint32_t n,
                                   int32_t N);

template <class T>
__global__ void oneactu_krnl_fast(T *g_idata, T *g_odata, int32_t nactu, T ampli,
                                  int32_t *xoff, int32_t *yoff, int32_t dim_im,
                                  int32_t dim_influ, int32_t N) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int32_t iy = i / dim_influ;
    int32_t ix = i - iy * dim_influ;
    int32_t ixactu = ix + xoff[nactu];
    int32_t iyactu = iy + yoff[nactu];

    // write result for this block to global mem
    if ((ixactu > -1) && (ixactu < dim_im) && (iyactu > -1) &&
        (iyactu < dim_im)) {
      int32_t tid = ixactu + iyactu * dim_im;
      g_odata[tid] = ampli * g_idata[i + nactu * N];
    }
    i += blockDim.x * gridDim.x;
  }
}

template <class T>
void oneactu(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t nactu,
             T ampli, int32_t *xoff, int32_t *yoff, int32_t dim_im, int32_t dim_influ, int32_t N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  oneactu_krnl_fast<T><<<dimGrid, dimBlock>>>(d_idata, d_odata, nactu, ampli,
                                              xoff, yoff, dim_im, dim_influ, N);

  carma_check_msg("oneactu_kernel<<<>>> execution failed\n");
}

template void oneactu<float>(int32_t threads, int32_t blocks, float *d_idata,
                             float *d_odata, int32_t nactu, float ampli, int32_t *xoff,
                             int32_t *yoff, int32_t dim_im, int32_t dim_influ, int32_t N);

template void oneactu<double>(int32_t threads, int32_t blocks, double *d_idata,
                              double *d_odata, int32_t nactu, double ampli,
                              int32_t *xoff, int32_t *yoff, int32_t dim_im, int32_t dim_influ,
                              int32_t N);

template <class T>
__global__ void oneactu_krnl_fast(T *g_idata, T *g_odata, int32_t nactu, T ampli,
                                  int32_t dim_im, int32_t dim_influ, int32_t N) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int32_t iy = i / dim_influ;
    int32_t ix = i - iy * dim_influ;

    // write result for this block to global mem
    if ((ix > -1) && (ix < dim_im) && (iy > -1) && (iy < dim_im)) {
      int32_t tid = ix + iy * dim_im;
      g_odata[tid] = ampli * g_idata[i + nactu * dim_influ * dim_influ];
    }
    i += blockDim.x * gridDim.x;
  }
}

template <class T>
void oneactu(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t nactu,
             T ampli, int32_t dim_im, int32_t dim_influ, int32_t N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  oneactu_krnl_fast<T><<<dimGrid, dimBlock>>>(d_idata, d_odata, nactu, ampli,
                                              dim_im, dim_influ, N);

  carma_check_msg("oneactu_kernel<<<>>> execution failed\n");
}

template void oneactu<float>(int32_t threads, int32_t blocks, float *d_idata,
                             float *d_odata, int32_t nactu, float ampli, int32_t dim_im,
                             int32_t dim_influ, int32_t N);

template void oneactu<double>(int32_t threads, int32_t blocks, double *d_idata,
                              double *d_odata, int32_t nactu, double ampli,
                              int32_t dim_im, int32_t dim_influ, int32_t N);

template <class T>
__global__ void fulldmshape_krnl(T *g_idata, T *g_odata, int32_t ninflu,
                                 int32_t diminflu, T *comm, int32_t N) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  uint32_t tid = threadIdx.x;
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  // TODO : replace if by while smartly
  if (i < N) {
    sdata[tid] = 0;

    for (int32_t cc = 0; cc < ninflu; cc++) {
      sdata[tid] += comm[cc] * g_idata[i + cc * diminflu];
    }
  }
  __syncthreads();

  if (i < N) {
    // write result for this block to global mem
    g_odata[i] = sdata[tid];
  }
}

template <class T>
void comp_fulldmshape(int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                      int32_t ninflu, int32_t diminflu, T *comm, int32_t N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int32_t smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  fulldmshape_krnl<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, ninflu,
                                                       diminflu, comm, N);

  carma_check_msg("fulldmshape_kernel<<<>>> execution failed\n");
}

template void comp_fulldmshape<float>(int32_t threads, int32_t blocks, float *d_idata,
                                      float *d_odata, int32_t ninflu, int32_t diminflu,
                                      float *comm, int32_t N);

template void comp_fulldmshape<double>(int32_t threads, int32_t blocks, double *d_idata,
                                       double *d_odata, int32_t ninflu,
                                       int32_t diminflu, double *comm, int32_t N);

template <class T>
__global__ void getIF_krnl(T *IF, float *dmshape, int32_t *indx_pup, int64_t nb_pts,
                           int64_t column, int64_t nb_col) {
  int32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  while (tid < nb_pts) {
    IF[column * nb_pts + tid] = dmshape[indx_pup[tid]];
    tid += blockDim.x * gridDim.x;
  }
}
template <class T>
__global__ void getIFfull_krnl(T *IF, float *dmshape, int64_t nb_pts, int64_t column,
                               int64_t nb_col) {
  int32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
  while (tid < nb_pts) {
    IF[column * nb_pts + tid] = dmshape[tid];
    tid += blockDim.x * gridDim.x;
  }
}
template <class T>
int32_t getIF(T *IF, float *dmshape, int32_t *indx_pup, int64_t nb_pts, int32_t column,
          int64_t nb_col, int32_t puponly, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, nb_pts, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  if (puponly)
    getIF_krnl<T>
        <<<grid, threads>>>(IF, dmshape, indx_pup, nb_pts, column, nb_col);
  else
    getIFfull_krnl<T><<<grid, threads>>>(IF, dmshape, nb_pts, column, nb_col);
  carma_check_msg("getIF_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int32_t getIF<float>(float *IF, float *dmshape, int32_t *indx_pup, int64_t nb_pts,
                          int32_t column, int64_t nb_col, int32_t puponly,
                          CarmaDevice *device);
template int32_t getIF<double>(double *IF, float *dmshape, int32_t *indx_pup,
                           int64_t nb_pts, int32_t column, int64_t nb_col, int32_t puponly,
                           CarmaDevice *device);

__global__ void do_statmat_krnl(float *statcov, float *xpos, float *ypos,
                                float norm, int64_t dim, int64_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t i, j;
  while (tid < N) {
    i = tid / dim;
    j = tid - i * dim;
    statcov[i * dim + j] =
        6.88 *
        powf(sqrtf((xpos[i] - xpos[j]) * (xpos[i] - xpos[j]) +
                   (ypos[i] - ypos[j]) * (ypos[i] - ypos[j])),
             5. / 3.) *
        norm;
    tid += blockDim.x * gridDim.x;
  }
}
int32_t dm_dostatmat(float *statcov, int64_t dim, float *xpos, float *ypos, float norm,
                 CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  int32_t N = (dim * dim);
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  do_statmat_krnl<<<grid, threads>>>(statcov, xpos, ypos, norm, dim, N);
  carma_check_msg("do_statcov_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fill_filtermat_krnl(float *filter, int32_t nactu, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    filter[tid] =
        tid % (nactu + 1) ? (float)-1. / nactu : (float)(1. - 1. / nactu);
    tid += blockDim.x * gridDim.x;
  }
}

int32_t fill_filtermat(float *filter, int32_t nactu, int32_t N, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_filtermat_krnl<<<grid, threads>>>(filter, nactu, N);
  carma_check_msg("fill_filtmat_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void convertToCom_krnl(uint16_t *volts, float *com, int32_t N,
                                  float volt_min, float volt_max, uint16_t val_max) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    com[tid] = float(volts[tid]) / float(val_max) * (volt_max - volt_min) + volt_min;
    tid += blockDim.x * gridDim.x;
  }
}

int32_t convertToCom(uint16_t *volts, float *com, int32_t N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device) {
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  convertToCom_krnl<<<grid, threads>>>(volts, com, N, volt_min, volt_max, val_max);
  carma_check_msg("convertToCom_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
