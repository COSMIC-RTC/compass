// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_centroider.cu
//! \ingroup   libsutra
//! \class     SutraCentroider
//! \brief     this class provides the centroider features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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

/*
 _                           _
 | |    __ _ _   _ _ __   ___| |__   ___ _ __ ___
 | |   / _` | | | | '_ \ / __| '_ \ / _ \ '__/ __|
 | |__| (_| | |_| | | | | (__| | | |  __/ |  \__ \
|_____\__,_|\__,_|_| |_|\___|_| |_|\___|_|  |___/

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
__global__ void calib_krnl(Tin *img_raw, float *img_cal, float *dark,
                           float *flat, int *lutPix, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    img_cal[tid] = (float(img_raw[lutPix[tid]]) - dark[tid]) * flat[tid];
    tid += blockDim.x * gridDim.x;
  }
}

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
