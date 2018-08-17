#include <sutra_centroider_tcog.h>
#include <carma_utils.cuh>

template <class T>
__global__ void centroidx(T *g_idata, T *g_odata, T *alpha, T thresh,
                          unsigned int n, unsigned int N, T scale, T offset,
                          unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  // unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  // unsigned int x = (tid % n) + 1;
  unsigned int x;
  int idim;
  sdata[tid] = 0;
  for (int cc = 0; cc < nelem_thread; cc++) {
    x = ((tid * nelem_thread + cc) % n);
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += (g_idata[idim] > thresh) ? g_idata[idim] * x : 0;
    else
      sdata[tid] += 0;
  }

  // if (i < N)
  //   sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] * x : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] =
        ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6)) - offset) * scale;
}

template <class T>
__global__ void centroidy(T *g_idata, T *g_odata, T *alpha, T thresh,
                          unsigned int n, unsigned int N, T scale, T offset,
                          unsigned int nelem_thread) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  // unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  // unsigned int y = (tid / n) + 1;
  unsigned int y;
  int idim;
  sdata[tid] = 0;
  for (int cc = 0; cc < nelem_thread; cc++) {
    y = ((tid * nelem_thread + cc) / n);
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < N)
      sdata[tid] += (g_idata[idim] > thresh) ? g_idata[idim] * y : 0;
    else
      sdata[tid] += 0;
  }

  // if (i < N)
  //   sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] * y : 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] =
        ((sdata[0] * 1.0 / (alpha[blockIdx.x] + 1.e-6)) - offset) * scale;
}

template <class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
                   T *d_odata, T *alpha, T thresh, T scale, T offset,
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
      d_idata, d_odata, alpha, thresh, n, size, scale, offset, nelem_thread);

  carmaCheckMsg("centroidx_kernel<<<>>> execution failed\n");

  centroidy<T><<<dimGrid, dimBlock, smemSize>>>(d_idata, &(d_odata[blocks]),
                                                alpha, thresh, n, size, scale,
                                                offset, nelem_thread);

  carmaCheckMsg("centroidy_kernel<<<>>> execution failed\n");
}

template void get_centroids<float>(int size, int threads, int blocks, int n,
                                   float *d_idata, float *d_odata, float *alpha,
                                   float thresh, float scale, float offset,
                                   carma_device *device);

template void get_centroids<double>(int size, int threads, int blocks, int n,
                                    double *d_idata, double *d_odata,
                                    double *alpha, double thresh, double scale,
                                    double offset, carma_device *device);
