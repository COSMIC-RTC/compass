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
                   carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  convert_krnl<<<grid, threads>>>(d_odata, d_idata, offset, scale, N);

  carmaCheckMsg("convert_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int convert_centro<float>(float *d_odata, float *d_idata, float offset,
                                   float scale, int N, carma_device *device);
template int convert_centro<double>(double *d_odata, double *d_idata,
                                    float offset, float scale, int N,
                                    carma_device *device);

int fill_validMask(int size, int npix, int blocks, int *d_validMask,
                   int *validx, int *validy, carma_device *device) {
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

  carmaCheckMsg("fillvalidMask_krnl<<<>>> execution failed\n");
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
                int *lutPix, int N, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  calib_krnl<<<grid, threads>>>(img_raw, img_cal, dark, flat, lutPix, N);

  carmaCheckMsg("calib_krnl<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int calibration<float>(float *img_raw, float *img_cal, float *dark,
                                float *flat, int *lutPix, int N,
                                carma_device *device);
template int calibration<uint16_t>(uint16_t *img_raw, float *img_cal,
                                   float *dark, float *flat, int *lutPix, int N,
                                   carma_device *device);
