#include <sutra_centroider_wcog.h>
#include <carma_utils.cuh>

template <class T>
__global__ void fillweights_krnl(T *d_out, T *weights, int Npix, int N) {
  int nim, idx;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / Npix;
    idx = tid - nim * Npix;
    d_out[tid] = weights[idx];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int fillweights(T *d_out, T *d_in, int npix, int N, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  fillweights_krnl<<<grid, threads>>>(d_out, d_in, npix * npix, N);
  carmaCheckMsg("<<<fillweights_krnl>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int fillweights<float>(float *d_out, float *d_in, int npix, int N,
                                carma_device *device);

template int fillweights<double>(double *d_out, double *d_in, int npix, int N,
                                 carma_device *device);

template <class T>
__global__ void centroidx(T *g_idata, T *g_odata, T *alpha, T *weights,
                          unsigned int n, unsigned int N, T scale, T offset,
                          unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  // unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  // unsigned int x = (tid % n) + 1;
  unsigned int x;
  int idim;
  sdata[tid] = 0.0f;
  for (int cc = 0; cc < nelem_thread; cc++) {
    x = ((tid * nelem_thread + cc) % n);
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += g_idata[idim] * x * weights[idim];
    else
      sdata[tid] += 0.0f;
  }
  __syncthreads();

  // sdata[tid] = (i < N) ? g_idata[i] * x * weights[i] : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);
  // if(tid == 0)
  //	printf("blockIdx %d sdata %f \n",blockIdx.x,sdata[tid]);
  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] =
        ((sdata[tid] * 1.0 / (alpha[blockIdx.x] + 1.e-6)) - offset) * scale;
}

template <class T>
__global__ void centroidy(T *g_idata, T *g_odata, T *alpha, T *weights,
                          unsigned int n, unsigned int N, T scale, T offset,
                          unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  // unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  // unsigned int y = (tid / n) + 1;
  unsigned int y;
  int idim;
  sdata[tid] = 0.0f;
  for (int cc = 0; cc < nelem_thread; cc++) {
    y = ((tid * nelem_thread + cc) / n);
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += g_idata[idim] * y * weights[idim];
    else
      sdata[tid] += 0.0f;
  }

  // sdata[tid] = (i < N) ? g_idata[i] * y * weights[i] : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] =
        ((sdata[tid] * 1.0 / (alpha[blockIdx.x] + 1.e-6)) - offset) * scale;
}

template <class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
                   T *d_odata, T *alpha, T *weights, T scale, T offset,
                   carma_device *device) {
  int maxThreads = device->get_properties().maxThreadsPerBlock;
  unsigned int nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  threads /= nelem_thread;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  centroidx<T><<<dimGrid, dimBlock, smemSize>>>(
      d_idata, d_odata, alpha, weights, n, size, scale, offset, nelem_thread);

  carmaCheckMsg("centroidx_kernel<<<>>> execution failed\n");

  centroidy<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, &(d_odata[blocks]),
                                                alpha, weights, n, size, scale,
                                                offset, nelem_thread);

  carmaCheckMsg("centroidy_kernel<<<>>> execution failed\n");
}

template void get_centroids<float>(int size, int threads, int blocks, int n,
                                   float *d_idata, float *d_odata, float *alpha,
                                   float *weights, float scale, float offset,
                                   carma_device *device);

template void get_centroids<double>(int size, int threads, int blocks, int n,
                                    double *d_idata, double *d_odata,
                                    double *alpha, double *weights,
                                    double scale, double offset,
                                    carma_device *device);
