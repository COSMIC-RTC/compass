#include <yoga_template.h>

template<class T>
__global__ void multim_krnl(T *odata, T* idata, int N) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] *= idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template<class T>
int multim(T *d_odata, T *i_data, int N) {

  //cerr << "hello : " << __LINE__ << endl;
  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, 0);

  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount * 8;
  int nThreads = (N + nBlocks - 1) / nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads - 1) / nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  multim_krnl<<<grid, threads>>>(d_odata, i_data, N);

  cutilCheckMsg("multim_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int
multim<float>(float *d_odata, float *d_idata, int N);
template int
multim<double>(double *d_odata, double *d_idata, int N);

