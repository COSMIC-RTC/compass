#include <sutra_centroider_maskedPix.h>
#include <carma_utils.cuh>

template <class T>
__global__ void get_maskedPix_krnl(T *g_odata, T *g_idata, int *subindx,
                                   int *subindy, T *intensities, int ns,
                                   int nslopes) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nslopes) {
    int i2 = subindx[i] + subindy[i] * ns;
    g_odata[i] = g_idata[i2] / intensities[0];
  }
}

template <class T>
void getMaskedPix(T *d_odata, T *d_idata, int *subindx, int *subindy,
                  T *intensities, int ns, int nslopes, carma_device *device) {
  // cout << "hello cu" << endl;

  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nslopes, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  get_maskedPix_krnl<T><<<grid, threads>>>(d_odata, d_idata, subindx, subindy,
                                           intensities, ns, nslopes);

  carmaCheckMsg("get_maskedPix_kernel<<<>>> execution failed\n");
}

template void getMaskedPix<float>(float *d_odata, float *d_idata, int *subindx,
                                  int *subindy, float *intensities, int ns,
                                  int nslopes, carma_device *device);
template void getMaskedPix<double>(double *d_odata, double *d_idata,
                                   int *subindx, int *subindy,
                                   double *intensities, int ns, int nslopes,
                                   carma_device *device);

template <class T>
__global__ void fill_intensities_krnl(T *g_odata, T *g_idata, int *subindx,
                                      int *subindy, int ns, int nslopes) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < nslopes) {
    int i = subindx[tid] + subindy[tid] * ns;
    g_odata[tid] = g_idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
void fill_intensities(T *d_odata, T *d_idata, int *subindx, int *subindy,
                      int ns, int nslopes, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nslopes, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  fill_intensities_krnl<T>
      <<<grid, threads>>>(d_odata, d_idata, subindx, subindy, ns, nslopes);

  carmaCheckMsg("fill_intensities_kernel<<<>>> execution failed\n");
}

template void fill_intensities<float>(float *intensities, float *cube,
                                      int *subindx, int *subindy, int ns,
                                      int nslopes, carma_device *device);
template void fill_intensities<double>(double *intensities, double *cube,
                                       int *subindx, int *subindy, int ns,
                                       int nslopes, carma_device *device);
