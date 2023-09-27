// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_maskedPix.cu
//! \ingroup   libsutra
//! \class     SutraCentroiderMaskedPix
//! \brief     this class provides the centroider_maskedPix features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_centroider_maskedPix.h>
#include <carma_utils.cuh>

template <class T>
__global__ void get_maskedPix_krnl(T *g_odata, T *ref, float *g_idata,
                                   int *subindx, int *subindy,
                                   float *intensities, int ns, int nslopes) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nslopes) {
    int i2 = subindx[i] + subindy[i] * ns;
    if(intensities[0] != 0)
      g_odata[i] = T(g_idata[i2] / intensities[0] * nslopes) - ref[i];
    else
      g_odata[i] = - ref[i];
  }
}

template <class T>
void get_masked_pix(T *d_odata, T *ref, float *d_idata, int *subindx,
                  int *subindy, float *intensities, int ns, int nslopes,
                  CarmaDevice *device, cudaStream_t stream) {
  // cout << "hello cu" << endl;

  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nslopes, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  get_maskedPix_krnl<T><<<grid, threads, 0, stream>>>(d_odata, ref, d_idata, subindx,
                                           subindy, intensities, ns, nslopes);

  carma_check_msg("get_maskedPix_kernel<<<>>> execution failed\n");
}

template void get_masked_pix<float>(float *d_odata, float *ref, float *d_idata,
                                  int *subindx, int *subindy,
                                  float *intensities, int ns, int nslopes,
                                  CarmaDevice *device, cudaStream_t stream);

#ifdef CAN_DO_HALF
template void get_masked_pix<half>(half *d_odata, half *ref, float *d_idata,
                                 int *subindx, int *subindy, float *intensities,
                                 int ns, int nslopes, CarmaDevice *device, cudaStream_t stream);
#endif

__global__ void fill_intensities_krnl(float *g_odata, float *g_idata,
                                      int *subindx, int *subindy, int ns,
                                      int nslopes) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < nslopes) {
    int i = subindx[tid] + subindy[tid] * ns;
    g_odata[tid] = g_idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

void fill_intensities(float *d_odata, float *d_idata, int *subindx,
                      int *subindy, int ns, int nslopes, CarmaDevice *device, cudaStream_t stream) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nslopes, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  fill_intensities_krnl<<<grid, threads, 0, stream>>>(d_odata, d_idata, subindx, subindy,
                                           ns, nslopes);

  carma_check_msg("fill_intensities_kernel<<<>>> execution failed\n");
}

template <typename T>
__global__ void pyr_fill_selected_pix_krnl(T *img, int img_sizex, T *pix,
                                           int *subindx, int *subindy,
                                           int nvalid) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < nvalid) {
    int pos = subindx[tid] + subindy[tid] * img_sizex;
    img[pos] = pix[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template <typename T>
void pyr_fill_selected_pix(T *img, int img_sizex, T *pix, int *subindx,
                           int *subindy, int nvalid, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  pyr_fill_selected_pix_krnl<T>
      <<<grid, threads>>>(img, img_sizex, pix, subindx, subindy, nvalid);
  carma_check_msg("pyr_fill_selected_pix_krnl<<<>>> execution failed\n");
}

template void pyr_fill_selected_pix<float>(float *img, int img_sizex,
                                           float *pix, int *subindx,
                                           int *subindy, int nvalid,
                                           CarmaDevice *device);
template void pyr_fill_selected_pix<double>(double *img, int img_sizex,
                                            double *pix, int *subindx,
                                            int *subindy, int nvalid,
                                            CarmaDevice *device);
#ifdef CAN_DO_HALF
template void pyr_fill_selected_pix<half>(half *img, int img_sizex, half *pix,
                                          int *subindx, int *subindy,
                                          int nvalid, CarmaDevice *device);
#endif

template <typename T>
__global__ void pyr_fill_mask_krnl(T *img, int img_sizex,
                                           int *subindx, int *subindy,
                                           int nvalid) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < nvalid) {
    int pos = subindx[tid] + subindy[tid] * img_sizex;
    img[pos] = T(1);
    tid += blockDim.x * gridDim.x;
  }
}

template <typename T>
void pyr_fill_mask(T *img, int img_sizex, int *subindx,
                           int *subindy, int nvalid, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, nvalid, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  pyr_fill_mask_krnl<T>
      <<<grid, threads>>>(img, img_sizex, subindx, subindy, nvalid);
  carma_check_msg("pyr_fill_mask_krnl<<<>>> execution failed\n");
}

template void pyr_fill_mask<float>(float *img, int img_sizex,
                                           int *subindx,
                                           int *subindy, int nvalid,
                                           CarmaDevice *device);
template void pyr_fill_mask<double>(double *img, int img_sizex,
                                            int *subindx,
                                            int *subindy, int nvalid,
                                            CarmaDevice *device);
#ifdef CAN_DO_HALF
template void pyr_fill_mask<half>(half *img, int img_sizex,
                                          int *subindx, int *subindy,
                                          int nvalid, CarmaDevice *device);
#endif
