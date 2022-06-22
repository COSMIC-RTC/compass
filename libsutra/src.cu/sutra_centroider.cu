// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider.cu
//! \ingroup   libsutra
//! \class     SutraCentroider
//! \brief     this class provides the centroider features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#include <sutra_centroider.h>
#include "carma_utils.cuh"

/**********************************
  _  __                    _      *
 | |/ /___ _ __ _ __   ___| |___  *
 | ' // _ \ '__| '_ \ / _ \ / __| *
 | . \  __/ |  | | | |  __/ \__ \ *
 |_|\_\___|_|  |_| |_|\___|_|___/ *
                                  *
 **********************************/
template <class T>
__global__ void convert_krnl(T *odata, T *idata, float offset, float scale,
                             int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = (idata[tid] - offset) * scale;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
__global__ void fillvalidMask_krnl(T *d_validMask, int *validx, int *validy,
                                   int npix, int size, int nelem_thread) {
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
      // d_validMask[idim] = blockIdx.x;
      atomicAdd(d_validMask + idim, blockIdx.x + 1);
    }
  }
}

template <class T>
__global__ void calibration_validPix_sh_krnl(T *img_raw, float *img_cal,
                                             float *dark, float *flat,
                                             int *lutPix, int *validx,
                                             int *validy, int npix, int size,
                                             int nelem_thread) {
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
      img_cal[idim] = (float(img_raw[lutPix[idim]]) - dark[idim]) * flat[idim];
    }
  }
}

template <class Tin>
__global__ void calibration_validPix_pyr_krnl(Tin *img_raw, float *img_cal,
                                              float *dark, float *flat,
                                              int *lutPix, int *validx,
                                              int *validy, int nvalid,
                                              int img_sizex) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < nvalid) {
    int pos = validx[tid] + validy[tid] * img_sizex;
    img_cal[pos] = (float(img_raw[lutPix[pos]]) - dark[pos]) * flat[pos];
    tid += blockDim.x * gridDim.x;
  }
}

template <class Tin>
__global__ void calib_krnl(Tin *img_raw, float *img_cal, float *dark,
                           float *flat, int *lutPix, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    img_cal[tid] = (float(img_raw[lutPix[tid]]) - dark[tid]) * flat[tid];
    tid += blockDim.x * gridDim.x;
  }
}

/*
 *  _                           _
 * | |    __ _ _   _ _ __   ___| |__   ___ _ __ ___
 * | |   / _` | | | | '_ \ / __| '_ \ / _ \ '__/ __|
 * | |__| (_| | |_| | | | | (__| | | |  __/ |  \__ \
 * |_____\__,_|\__,_|_| |_|\___|_| |_|\___|_|  |___/
 *
 */
template <class T>
int convert_centro(T *d_odata, T *d_idata, float offset, float scale, int N,
                   CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  convert_krnl<<<grid, threads>>>(d_odata, d_idata, offset, scale, N);

  carma_check_msg("convert_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int convert_centro<float>(float *d_odata, float *d_idata, float offset,
                                   float scale, int N, CarmaDevice *device);
template int convert_centro<double>(double *d_odata, double *d_idata,
                                    float offset, float scale, int N,
                                    CarmaDevice *device);

int fill_validMask(int size, int npix, int blocks, int *d_validMask,
                   int *validx, int *validy, CarmaDevice *device) {
  int maxThreads = device->get_properties().maxThreadsPerBlock;
  int threads = npix * npix;
  unsigned int nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  threads /= nelem_thread;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  fillvalidMask_krnl<<<dimGrid, dimBlock>>>(d_validMask, validx, validy, npix,
                                            size, nelem_thread);

  carma_check_msg("fillvalidMask_krnl<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template <class Tin>
int calibration_validPix_sh(int npix, int size, int blocks, Tin *img_raw,
                            float *img_cal, float *dark, float *flat,
                            int *lutPix, int *validx, int *validy,
                            CarmaDevice *device, cudaStream_t stream) {
  int maxThreads = device->get_properties().maxThreadsPerBlock;
  int threads = npix * npix;
  unsigned int nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  threads /= nelem_thread;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  calibration_validPix_sh_krnl<<<dimGrid, dimBlock, 0, stream>>>(
      img_raw, img_cal, dark, flat, lutPix, validx, validy, npix, size,
      nelem_thread);

  carma_check_msg("calibration_validPix_sh_krnl<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int calibration_validPix_sh<float>(int npix, int size, int block,
                                            float *img_raw, float *img_cal,
                                            float *dark, float *flat,
                                            int *lutPix, int *validx,
                                            int *validy, CarmaDevice *device,
                                            cudaStream_t stream);
template int calibration_validPix_sh<uint16_t>(
    int npix, int size, int block, uint16_t *img_raw, float *img_cal,
    float *dark, float *flat, int *lutPix, int *validx, int *validy,
    CarmaDevice *device, cudaStream_t stream);

template <class Tin>
int calibration_validPix_pyr(Tin *img_raw, float *img_cal, float *dark,
                             float *flat, int *lutPix, int *validx, int *validy,
                             int nvalid, int img_sizex, CarmaDevice *device,
                             cudaStream_t stream) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  calibration_validPix_pyr_krnl<<<grid, threads, 0, stream>>>(
      img_raw, img_cal, dark, flat, lutPix, validx, validy, nvalid, img_sizex);

  carma_check_msg("calibration_validPix_pyr_krnl<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int calibration_validPix_pyr<float>(float *img_raw, float *img_cal,
                                             float *dark, float *flat,
                                             int *lutPix, int *validx,
                                             int *validy, int nvalid,
                                             int img_sizex, CarmaDevice *device,
                                             cudaStream_t stream);
template int calibration_validPix_pyr<uint16_t>(
    uint16_t *img_raw, float *img_cal, float *dark, float *flat, int *lutPix,
    int *validx, int *validy, int nvalid, int img_sizex, CarmaDevice *device,
    cudaStream_t stream);

template <class Tin>
int calibration(Tin *img_raw, float *img_cal, float *dark, float *flat,
                int *lutPix, int N, CarmaDevice *device, cudaStream_t stream) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  calib_krnl<<<grid, threads, 0, stream>>>(img_raw, img_cal, dark, flat, lutPix, N);

  carma_check_msg("calib_krnl<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int calibration<float>(float *img_raw, float *img_cal, float *dark,
                                float *flat, int *lutPix, int N,
                                CarmaDevice *device, cudaStream_t stream);
template int calibration<uint16_t>(uint16_t *img_raw, float *img_cal,
                                   float *dark, float *flat, int *lutPix, int N,
                                   CarmaDevice *device, cudaStream_t stream);
