#include <sutra_centroider_maskedPix.h>
#include <carma_utils.cuh>

template <class T>
__global__ void get_maskedPix_krnl(T *g_odata, T *g_idata, int *subindx,
                                   int *subindy, T *subsum, int ns,
                                   int nslopes) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nslopes) {
    int i2 = subindx[i] + subindy[i] * ns;
    g_odata[i] = g_idata[i2] / subsum[0];
  }
}

template <class T>
void getMaskedPix(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum,
                  int ns, int nslopes, carma_device *device) {
  // cout << "hello cu" << endl;

  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nslopes, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  get_maskedPix_krnl<T><<<grid, threads>>>(d_odata, d_idata, subindx, subindy,
                                           subsum, ns, nslopes);

  carmaCheckMsg("get_maskedPix_kernel<<<>>> execution failed\n");
}

template void getMaskedPix<float>(float *d_odata, float *d_idata, int *subindx,
                                  int *subindy, float *subsum, int ns,
                                  int nslopes, carma_device *device);
template void getMaskedPix<double>(double *d_odata, double *d_idata,
                                   int *subindx, int *subindy, double *subsum,
                                   int ns, int nslopes, carma_device *device);

template <class T>
__global__ void fill_subsum_krnl(T *g_odata, T *g_idata, int *subindx,
                                 int *subindy, int ns, int nslopes) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < nslopes) {
    int i = subindx[tid] + subindy[tid] * ns;
    g_odata[tid] = g_idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
void fill_subsum(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns,
                 int nslopes, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nslopes, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  fill_subsum_krnl<T>
      <<<grid, threads>>>(d_odata, d_idata, subindx, subindy, ns, nslopes);

  carmaCheckMsg("fill_subsum_kernel<<<>>> execution failed\n");
}

template void fill_subsum<float>(float *subsum, float *cube, int *subindx,
                                 int *subindy, int ns, int nslopes,
                                 carma_device *device);
template void fill_subsum<double>(double *subsum, double *cube, int *subindx,
                                  int *subindy, int ns, int nslopes,
                                  carma_device *device);
