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
__global__ void convert_krnl(T *odata, T *idata, T offset, T scale, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = (idata[tid] - offset) * scale;
    tid += blockDim.x * gridDim.x;
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
int convert_centro(T *d_odata, T *d_idata, T offset, T scale, int N,
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
                                    double offset, double scale, int N,
                                    carma_device *device);
