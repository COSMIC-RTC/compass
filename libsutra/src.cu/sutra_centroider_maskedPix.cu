#include <sutra_centroider_maskedPix.h>
#include <carma_utils.cuh>

template <class T>
__global__ void get_maskedPix_krnl(T *g_odata, T *ref, float *g_idata, int *subindx,
                                   int *subindy, float *intensities, int ns,
                                   int nslopes) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nslopes) {
    int i2 = subindx[i] + subindy[i] * ns;
    g_odata[i] = T(g_idata[i2] / intensities[0]) - ref[i];
  }
}

template <class T>
void getMaskedPix(T *d_odata, T *ref, float *d_idata, int *subindx, int *subindy,
                  float *intensities, int ns, int nslopes, carma_device *device) {
  // cout << "hello cu" << endl;

  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nslopes, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  get_maskedPix_krnl<T><<<grid, threads>>>(d_odata, ref, d_idata, subindx,
                                           subindy, intensities, ns, nslopes);

  carmaCheckMsg("get_maskedPix_kernel<<<>>> execution failed\n");
}

template void getMaskedPix<float>(float *d_odata, float *ref, float *d_idata,
                                  int *subindx, int *subindy,
                                  float *intensities, int ns, int nslopes,
                                  carma_device *device);

#ifdef CAN_DO_HALF
template void getMaskedPix<half>(half *d_odata, half *ref, float *d_idata,
  int *subindx, int *subindy,
  float *intensities, int ns, int nslopes,
  carma_device *device);
#endif

__global__ void fill_intensities_krnl(float *g_odata, float *g_idata, int *subindx,
                                      int *subindy, int ns, int nslopes) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < nslopes) {
    int i = subindx[tid] + subindy[tid] * ns;
    g_odata[tid] = g_idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

void fill_intensities(float *d_odata, float *d_idata, int *subindx, int *subindy,
                      int ns, int nslopes, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nslopes, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  fill_intensities_krnl
      <<<grid, threads>>>(d_odata, d_idata, subindx, subindy, ns, nslopes);

  carmaCheckMsg("fill_intensities_kernel<<<>>> execution failed\n");
}
