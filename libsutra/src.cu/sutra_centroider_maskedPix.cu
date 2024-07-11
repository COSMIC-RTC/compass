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

//! \file      sutra_centroider_maskedPix.cu
//! \ingroup   libsutra
//! \class     SutraCentroiderMaskedPix
//! \brief     this class provides the centroider_maskedPix features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_centroider_maskedPix.hpp>
#include <carma_utils.cuh>

template <class T>
__global__ void get_maskedPix_krnl(T *g_odata, T *ref, float *g_idata,
                                   int32_t *subindx, int32_t *subindy,
                                   float *intensities, int32_t ns, int32_t nslopes) {
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nslopes) {
    int32_t i2 = subindx[i] + subindy[i] * ns;
    if(intensities[0] != 0)
      g_odata[i] = T(g_idata[i2] / intensities[0] * nslopes) - ref[i];
    else
      g_odata[i] = - ref[i];
  }
}

template <class T>
void get_masked_pix(T *d_odata, T *ref, float *d_idata, int32_t *subindx,
                  int32_t *subindy, float *intensities, int32_t ns, int32_t nslopes,
                  CarmaDevice *device, cudaStream_t stream) {
  // cout << "hello cu" << endl;

  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nslopes, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  get_maskedPix_krnl<T><<<grid, threads, 0, stream>>>(d_odata, ref, d_idata, subindx,
                                           subindy, intensities, ns, nslopes);

  carma_check_msg("get_maskedPix_kernel<<<>>> execution failed\n");
}

template void get_masked_pix<float>(float *d_odata, float *ref, float *d_idata,
                                  int32_t *subindx, int32_t *subindy,
                                  float *intensities, int32_t ns, int32_t nslopes,
                                  CarmaDevice *device, cudaStream_t stream);

#ifdef CAN_DO_HALF
template void get_masked_pix<half>(half *d_odata, half *ref, float *d_idata,
                                 int32_t *subindx, int32_t *subindy, float *intensities,
                                 int32_t ns, int32_t nslopes, CarmaDevice *device, cudaStream_t stream);
#endif

__global__ void fill_intensities_krnl(float *g_odata, float *g_idata,
                                      int32_t *subindx, int32_t *subindy, int32_t ns,
                                      int32_t nslopes) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < nslopes) {
    int32_t i = subindx[tid] + subindy[tid] * ns;
    g_odata[tid] = g_idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

void fill_intensities(float *d_odata, float *d_idata, int32_t *subindx,
                      int32_t *subindy, int32_t ns, int32_t nslopes, CarmaDevice *device, cudaStream_t stream) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nslopes, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_intensities_krnl<<<grid, threads, 0, stream>>>(d_odata, d_idata, subindx, subindy,
                                           ns, nslopes);

  carma_check_msg("fill_intensities_kernel<<<>>> execution failed\n");
}

template <typename T>
__global__ void pyr_fill_selected_pix_krnl(T *img, int32_t img_sizex, T *pix,
                                           int32_t *subindx, int32_t *subindy,
                                           int32_t nvalid) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < nvalid) {
    int32_t pos = subindx[tid] + subindy[tid] * img_sizex;
    img[pos] = pix[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template <typename T>
void pyr_fill_selected_pix(T *img, int32_t img_sizex, T *pix, int32_t *subindx,
                           int32_t *subindy, int32_t nvalid, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  pyr_fill_selected_pix_krnl<T>
      <<<grid, threads>>>(img, img_sizex, pix, subindx, subindy, nvalid);
  carma_check_msg("pyr_fill_selected_pix_krnl<<<>>> execution failed\n");
}

template void pyr_fill_selected_pix<float>(float *img, int32_t img_sizex,
                                           float *pix, int32_t *subindx,
                                           int32_t *subindy, int32_t nvalid,
                                           CarmaDevice *device);
template void pyr_fill_selected_pix<double>(double *img, int32_t img_sizex,
                                            double *pix, int32_t *subindx,
                                            int32_t *subindy, int32_t nvalid,
                                            CarmaDevice *device);
#ifdef CAN_DO_HALF
template void pyr_fill_selected_pix<half>(half *img, int32_t img_sizex, half *pix,
                                          int32_t *subindx, int32_t *subindy,
                                          int32_t nvalid, CarmaDevice *device);
#endif

template <typename T>
__global__ void pyr_fill_mask_krnl(T *img, int32_t img_sizex,
                                           int32_t *subindx, int32_t *subindy,
                                           int32_t nvalid) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < nvalid) {
    int32_t pos = subindx[tid] + subindy[tid] * img_sizex;
    img[pos] = T(1);
    tid += blockDim.x * gridDim.x;
  }
}

template <typename T>
void pyr_fill_mask(T *img, int32_t img_sizex, int32_t *subindx,
                           int32_t *subindy, int32_t nvalid, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  pyr_fill_mask_krnl<T>
      <<<grid, threads>>>(img, img_sizex, subindx, subindy, nvalid);
  carma_check_msg("pyr_fill_mask_krnl<<<>>> execution failed\n");
}

template void pyr_fill_mask<float>(float *img, int32_t img_sizex,
                                           int32_t *subindx,
                                           int32_t *subindy, int32_t nvalid,
                                           CarmaDevice *device);
template void pyr_fill_mask<double>(double *img, int32_t img_sizex,
                                            int32_t *subindx,
                                            int32_t *subindy, int32_t nvalid,
                                            CarmaDevice *device);
#ifdef CAN_DO_HALF
template void pyr_fill_mask<half>(half *img, int32_t img_sizex,
                                          int32_t *subindx, int32_t *subindy,
                                          int32_t nvalid, CarmaDevice *device);
#endif
