// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_dm.cu
//! \ingroup   libsutra
//! \class     SutraDm
//! \brief     this class provides the dm features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#include <sutra_dm.h>
#include "carma_utils.cuh"

#include <cuda.h>

#include <stdio.h>

__device__ inline struct doubleint getInterval(const int pos,
                                               const int *iStart_t) {
  const int start = iStart_t[pos];
  return {start, iStart_t[pos + 1] - start};
}

template <class T>
__device__ inline T get_data(const int pos, const T *cmdVector,
                            const T *influData, const int *iPos) {
  return cmdVector[iPos[pos]] * influData[pos];
}

template <class T>
__global__ void compShape2(T *outData, const T *cmdVector, const T *influData,
                           const int *iPos, const int *iStart_t, const int N) {
  const int id = threadIdx.x + blockIdx.x * blockDim.x;

  if (id < N) {
    const struct doubleint interval = getInterval(id, iStart_t);

    T sum = 0;

    for (int pos = interval.start; pos < interval.start + interval.nbInflu;
         ++pos)
      sum += get_data(pos, cmdVector, influData, iPos);

    outData[id] = sum;
  }
}

template <class T>
void comp_dmshape2(T *outData, const T *cmdVector, const T *influData,
                   const int *iStart_t, const int *iPos, const int roiLength,
                   const dim3 threads, const dim3 blocks, const int shared) {
  compShape2<T><<<blocks, threads>>>(outData, cmdVector, influData, iPos,
                                     iStart_t, roiLength);

  carma_check_msg("comp_dmshape2<<<>>> execution failed\n");
}

template void comp_dmshape2<float>(float *outData, const float *cmdVector,
                                   const float *influData, const int *iStart_t,
                                   const int *iPos, const int roiLength,
                                   const dim3 threads, const dim3 blocks,
                                   const int shared);

template void comp_dmshape2<double>(double *outData, const double *cmdVector,
                                    const double *influData,
                                    const int *iStart_t, const int *iPos,
                                    const int roiLength, const dim3 threads,
                                    const dim3 blocks, const int shared);

template <class T>
__global__ void dmshape_krnl(T *g_idata, T *g_odata, int *pos, int *istart,
                             int *npts, T *comm, unsigned int n, int N) {
  // T *sdata = SharedMemory<T>();
  // load shared mem
  // unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int local_istart = istart[i];
    int local_npts = npts[i];

    // sdata[tid] = 0;
    T value = 0;

    if (local_npts > 0) {
      for (int cc = 0; cc < local_npts; cc++) {
        int lpos = pos[local_istart + cc];
        int ninflu = lpos / n;
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
void comp_dmshape(int threads, int blocks, T *d_idata, T *d_odata, int *pos,
                  int *istart, int *npts, T *comm, unsigned int n, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  // int smemSize =
  //    (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

  dmshape_krnl<T><<<dimGrid, dimBlock /*, smemSize*/>>>(
      d_idata, d_odata, pos, istart, npts, comm, n, N);

  carma_check_msg("dmshape_kernel<<<>>> execution failed\n");
}

template void comp_dmshape<float>(int threads, int blocks, float *d_idata,
                                  float *d_odata, int *pos, int *istart,
                                  int *npts, float *comm, unsigned int n,
                                  int N);

template void comp_dmshape<double>(int threads, int blocks, double *d_idata,
                                   double *d_odata, int *pos, int *istart,
                                   int *npts, double *comm, unsigned int n,
                                   int N);

template <class T>
__global__ void oneactu_krnl_fast(T *g_idata, T *g_odata, int nactu, T ampli,
                                  int *xoff, int *yoff, int dim_im,
                                  int dim_influ, int N) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int iy = i / dim_influ;
    int ix = i - iy * dim_influ;
    int ixactu = ix + xoff[nactu];
    int iyactu = iy + yoff[nactu];

    // write result for this block to global mem
    if ((ixactu > -1) && (ixactu < dim_im) && (iyactu > -1) &&
        (iyactu < dim_im)) {
      int tid = ixactu + iyactu * dim_im;
      g_odata[tid] = ampli * g_idata[i + nactu * N];
    }
    i += blockDim.x * gridDim.x;
  }
}

template <class T>
void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu,
             T ampli, int *xoff, int *yoff, int dim_im, int dim_influ, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  oneactu_krnl_fast<T><<<dimGrid, dimBlock>>>(d_idata, d_odata, nactu, ampli,
                                              xoff, yoff, dim_im, dim_influ, N);

  carma_check_msg("oneactu_kernel<<<>>> execution failed\n");
}

template void oneactu<float>(int threads, int blocks, float *d_idata,
                             float *d_odata, int nactu, float ampli, int *xoff,
                             int *yoff, int dim_im, int dim_influ, int N);

template void oneactu<double>(int threads, int blocks, double *d_idata,
                              double *d_odata, int nactu, double ampli,
                              int *xoff, int *yoff, int dim_im, int dim_influ,
                              int N);

template <class T>
__global__ void oneactu_krnl_fast(T *g_idata, T *g_odata, int nactu, T ampli,
                                  int dim_im, int dim_influ, int N) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < N) {
    int iy = i / dim_influ;
    int ix = i - iy * dim_influ;

    // write result for this block to global mem
    if ((ix > -1) && (ix < dim_im) && (iy > -1) && (iy < dim_im)) {
      int tid = ix + iy * dim_im;
      g_odata[tid] = ampli * g_idata[i + nactu * dim_influ * dim_influ];
    }
    i += blockDim.x * gridDim.x;
  }
}

template <class T>
void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu,
             T ampli, int dim_im, int dim_influ, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  oneactu_krnl_fast<T><<<dimGrid, dimBlock>>>(d_idata, d_odata, nactu, ampli,
                                              dim_im, dim_influ, N);

  carma_check_msg("oneactu_kernel<<<>>> execution failed\n");
}

template void oneactu<float>(int threads, int blocks, float *d_idata,
                             float *d_odata, int nactu, float ampli, int dim_im,
                             int dim_influ, int N);

template void oneactu<double>(int threads, int blocks, double *d_idata,
                              double *d_odata, int nactu, double ampli,
                              int dim_im, int dim_influ, int N);

template <class T>
__global__ void fulldmshape_krnl(T *g_idata, T *g_odata, int ninflu,
                                 int diminflu, T *comm, int N) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  // TODO : replace if by while smartly
  if (i < N) {
    sdata[tid] = 0;

    for (int cc = 0; cc < ninflu; cc++) {
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
void comp_fulldmshape(int threads, int blocks, T *d_idata, T *d_odata,
                      int ninflu, int diminflu, T *comm, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  fulldmshape_krnl<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, ninflu,
                                                       diminflu, comm, N);

  carma_check_msg("fulldmshape_kernel<<<>>> execution failed\n");
}

template void comp_fulldmshape<float>(int threads, int blocks, float *d_idata,
                                      float *d_odata, int ninflu, int diminflu,
                                      float *comm, int N);

template void comp_fulldmshape<double>(int threads, int blocks, double *d_idata,
                                       double *d_odata, int ninflu,
                                       int diminflu, double *comm, int N);

template <class T>
__global__ void getIF_krnl(T *IF, float *dmshape, int *indx_pup, long nb_pts,
                           long column, long nb_col) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  while (tid < nb_pts) {
    IF[column * nb_pts + tid] = dmshape[indx_pup[tid]];
    tid += blockDim.x * gridDim.x;
  }
}
template <class T>
__global__ void getIFfull_krnl(T *IF, float *dmshape, long nb_pts, long column,
                               long nb_col) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  while (tid < nb_pts) {
    IF[column * nb_pts + tid] = dmshape[tid];
    tid += blockDim.x * gridDim.x;
  }
}
template <class T>
int getIF(T *IF, float *dmshape, int *indx_pup, long nb_pts, int column,
          long nb_col, int puponly, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
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
template int getIF<float>(float *IF, float *dmshape, int *indx_pup, long nb_pts,
                          int column, long nb_col, int puponly,
                          CarmaDevice *device);
template int getIF<double>(double *IF, float *dmshape, int *indx_pup,
                           long nb_pts, int column, long nb_col, int puponly,
                           CarmaDevice *device);

__global__ void do_statmat_krnl(float *statcov, float *xpos, float *ypos,
                                float norm, long dim, long N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int i, j;
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
int dm_dostatmat(float *statcov, long dim, float *xpos, float *ypos, float norm,
                 CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  int N = (dim * dim);
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  do_statmat_krnl<<<grid, threads>>>(statcov, xpos, ypos, norm, dim, N);
  carma_check_msg("do_statcov_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fill_filtermat_krnl(float *filter, int nactu, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    filter[tid] =
        tid % (nactu + 1) ? (float)-1. / nactu : (float)(1. - 1. / nactu);
    tid += blockDim.x * gridDim.x;
  }
}

int fill_filtermat(float *filter, int nactu, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_filtermat_krnl<<<grid, threads>>>(filter, nactu, N);
  carma_check_msg("fill_filtmat_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void convertToCom_krnl(uint16_t *volts, float *com, int N,
                                  float volt_min, float volt_max, uint16_t val_max) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    com[tid] = float(volts[tid]) / float(val_max) * (volt_max - volt_min) + volt_min;
    tid += blockDim.x * gridDim.x;
  }
}

int convertToCom(uint16_t *volts, float *com, int N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  convertToCom_krnl<<<grid, threads>>>(volts, com, N, volt_min, volt_max, val_max);
  carma_check_msg("convertToCom_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
