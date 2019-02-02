#include <sutra_centroider_bpcog.h>
#include <carma_utils.cuh>

template <typename T, int BLOCK_THREADS>
__launch_bounds__(BLOCK_THREADS) __global__
    void centroids(T *d_img, T *d_centroids, T *ref, int *validx, int *validy,
                   T *d_intensities, int nbpix, unsigned int npix,
                   unsigned int size, T scale, T offset,
                   unsigned int nelem_thread) {
  // Specialize BlockRadixSort for a 1D block of BLOCK_THREADS threads owning 1
  // item each
  typedef cub::BlockRadixSort<T, BLOCK_THREADS, 1> BlockRadixSortT;
  typedef cub::BlockReduce<T, BLOCK_THREADS> BlockReduce;

  // Allocate shared memory for BlockRadixSort
  __shared__ typename BlockRadixSortT::TempStorage temp_storageSort;
  __shared__ typename BlockReduce::TempStorage temp_storageSum;
  __shared__ T threshold;

  T idata = 0;
  T xdata = 0;
  T ydata = 0;

  unsigned int tid = threadIdx.x;
  unsigned int xvalid = validx[blockIdx.x];
  unsigned int yvalid = validy[blockIdx.x];
  unsigned int x = tid % npix;
  unsigned int y = tid / npix;
  int idim = (x + xvalid) + (y + yvalid) * size;

  T items[1];
  items[0] = d_img[idim];

  __syncthreads();
  BlockRadixSortT(temp_storageSort).SortDescending(items);

  if (tid == nbpix) threshold = items[0];

  __syncthreads();

  if (idim < size * size) {
    const T data_thresh = d_img[idim] - threshold;
    idata += (data_thresh > 0) ? data_thresh : 0;
    xdata += (data_thresh > 0) ? data_thresh * x : 0;
    ydata += (data_thresh > 0) ? data_thresh * y : 0;
    d_img[idim] = data_thresh;
  }

  __syncthreads();
  T intensity = BlockReduce(temp_storageSum).Sum(idata, blockDim.x);
  T slopex = BlockReduce(temp_storageSum).Sum(xdata, blockDim.x);
  T slopey = BlockReduce(temp_storageSum).Sum(ydata, blockDim.x);
  if (tid == 0) {
    d_centroids[blockIdx.x] =
        ((slopex * 1.0 / (intensity + 1.e-6)) - offset) * scale -
        ref[blockIdx.x];
    d_centroids[blockIdx.x + gridDim.x] =
        ((slopey * 1.0 / (intensity + 1.e-6)) - offset) * scale -
        ref[blockIdx.x + gridDim.x];
    d_intensities[blockIdx.x] = intensity;
  }
}

template <class T>
void get_centroids(int size, int threads, int blocks, int npix, T *d_img,
                   T *d_centroids, T *ref, int *validx, int *validy,
                   T *intensities, int nbpix, T scale, T offset,
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
  if (threads <= 16)
    centroids<T, 16><<<dimGrid, dimBlock>>>(d_img, d_centroids, ref, validx,
                                            validy, intensities, nbpix, npix,
                                            size, scale, offset, nelem_thread);
  else if (threads <= 32)
    centroids<T, 32><<<dimGrid, dimBlock>>>(d_img, d_centroids, ref, validx,
                                            validy, intensities, nbpix, npix,
                                            size, scale, offset, nelem_thread);

  else if (threads <= 64)
    centroids<T, 64><<<dimGrid, dimBlock>>>(d_img, d_centroids, ref, validx,
                                            validy, intensities, nbpix, npix,
                                            size, scale, offset, nelem_thread);
  else if (threads <= 128)
    centroids<T, 128><<<dimGrid, dimBlock>>>(d_img, d_centroids, ref, validx,
                                             validy, intensities, nbpix, npix,
                                             size, scale, offset, nelem_thread);
  else if (threads <= 256)
    centroids<T, 256><<<dimGrid, dimBlock>>>(d_img, d_centroids, ref, validx,
                                             validy, intensities, nbpix, npix,
                                             size, scale, offset, nelem_thread);
  else if (threads <= 512)
    centroids<T, 512><<<dimGrid, dimBlock>>>(d_img, d_centroids, ref, validx,
                                             validy, intensities, nbpix, npix,
                                             size, scale, offset, nelem_thread);
  else
    printf("SH way too big !!!\n");

  carmaCheckMsg("centroids_kernel<<<>>> execution failed\n");
}

template void get_centroids<float>(int size, int threads, int blocks, int npix,
                                   float *d_img, float *d_centroids, float *ref,
                                   int *validx, int *validy, float *intensities,
                                   int nbpix, float scale, float offset,
                                   carma_device *device);

template void get_centroids<double>(int size, int threads, int blocks, int npix,
                                    double *d_img, double *d_centroids,
                                    double *ref, int *validx, int *validy,
                                    double *intensities, int nbpix,
                                    double scale, double offset,
                                    carma_device *device);

template <class T>
__device__ inline void sortmax_krnl(T *sdata, unsigned int *values, int size,
                                    int n) {
  if (!((size & (size - 1)) == 0)) {
    unsigned int s;
    if ((size & 1) != 0)
      s = size / 2 + 1;  //(size&1)==size%2
    else
      s = size / 2;
    unsigned int s_old = size;
    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2 * s < s_old) && (s != 0)) s += 1;
    }
  } else {
    // do reduction in shared mem
    for (unsigned int s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        if (sdata[n] < sdata[n + s]) {
          mswap(values[n], values[n + s]);
          mswap(sdata[n], sdata[n + s]);
        }
      }
      __syncthreads();
    }
  }
}

template <class T>
__global__ void sortmax(T *g_idata, T *g_odata, unsigned int *values, int nmax,
                        int Npix, int size, int nelem_thread) {
  extern __shared__ uint svalues[];
  T *sdata = (T *)&svalues[Npix];
  /*
    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    svalues[tid] = tid;
    sdata[tid] = g_idata[i];
  */
  unsigned int tid = threadIdx.x;
  // unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  // unsigned int y = (tid / n) + 1;
  int idim;
  int sdim;

  for (int cc = 0; cc < nelem_thread; cc++) {
    idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    sdim = tid * nelem_thread + cc;
    if (idim < size) {
      sdata[sdim] = g_idata[idim];
      svalues[sdim] = sdim;
    }
  }

  __syncthreads();

  for (int cc = 0; cc < nmax; cc++) {
    for (int cpt = 0; cpt < nelem_thread; cpt++) {
      sdim = tid * nelem_thread + cpt;
      if (sdim >= cc)
        sortmax_krnl(&(sdata[cc]), &(svalues[cc]), Npix - cc, sdim - cc);
      __syncthreads();
    }
  }
  for (int cpt = 0; cpt < nelem_thread; cpt++) {
    sdim = tid * nelem_thread + cpt;
    if (sdim < nmax) {
      g_odata[nmax * blockIdx.x + sdim] = sdata[sdim];  // - sdata[nmax - 1];
      values[nmax * blockIdx.x + sdim] = svalues[sdim];
    }
  }
  /*
   __syncthreads();
   if ((blockIdx.x == 0) && (tid < nmax))
   printf("tid %d sdata %f \n",tid,g_odata[tid]);
   */
}

template <class T>
void subap_sortmax(int threads, int blocks, T *d_idata, T *d_odata,
                   unsigned int *values, int nmax, carma_device *device) {
  int maxThreads = device->get_properties().maxThreadsPerBlock;
  unsigned int nelem_thread = 1;
  while ((threads / nelem_thread > maxThreads) ||
         (threads % nelem_thread != 0)) {
    nelem_thread++;
  }

  dim3 dimBlock(threads / nelem_thread, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  size_t smemSize = threads * (sizeof(T) + sizeof(uint));
  sortmax<T><<<dimGrid, dimBlock, smemSize>>>(
      d_idata, d_odata, values, nmax, threads, threads * blocks, nelem_thread);

  carmaCheckMsg("sortmax_kernel<<<>>> execution failed\n");
}
template void subap_sortmax<float>(int threads, int blocks, float *d_idata,
                                   float *d_odata, unsigned int *values,
                                   int nmax, carma_device *device);
template void subap_sortmax<double>(int threads, int blocks, double *d_idata,
                                    double *d_odata, unsigned int *values,
                                    int nmax, carma_device *device);

template <class T>
__global__ void centroid_bpix(int nsub, int n, T *g_idata, unsigned int *values,
                              T *g_odata, T scale, T offset) {
  extern __shared__ uint svalues[];
  T *sdata = (T *)&svalues[blockDim.x];
  T intensities;
  // T minimum;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  svalues[tid] = values[i];
  sdata[tid] = g_idata[i] - g_idata[blockIdx.x * blockDim.x + blockDim.x - 1];

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);

  __syncthreads();
  // get the sum per subap
  if (tid == 0) intensities = (abs(sdata[tid]) > 1.e-6 ? sdata[tid] : 0.0f);

  __syncthreads();

  // Reload sdata
  sdata[tid] = g_idata[i] - g_idata[blockIdx.x * blockDim.x + blockDim.x - 1];

  __syncthreads();

  // compute the centroid on the first part of the array
  sdata[tid] *= ((svalues[tid] % n));
  // x centroid
  __syncthreads();
  reduce_krnl(sdata, blockDim.x, tid);
  //__syncthreads();
  if (tid == 0)
    g_odata[blockIdx.x] =
        (intensities != 0.0f ? ((sdata[tid] / intensities) - offset) * scale
                             : 0.0f);
  __syncthreads();
  sdata[tid] = g_idata[i] - g_idata[blockIdx.x * blockDim.x + blockDim.x - 1];

  __syncthreads();

  // compute the centroid on the first part of the array
  sdata[tid] *= (svalues[tid] / n);
  // y centroid
  __syncthreads();

  reduce_krnl(sdata, blockDim.x, tid);
  //__syncthreads();
  if (tid == 0)
    g_odata[blockIdx.x + nsub] =
        (intensities != 0.0f ? ((sdata[tid] / intensities) - offset) * scale
                             : 0.0f);
}

template <class T>
void subap_bpcentro(int threads, int blocks, int npix, T *d_idata,
                    unsigned int *values, T *d_odata, T scale, T offset) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = threads * (sizeof(T) + sizeof(uint));

  centroid_bpix<T><<<dimGrid, dimBlock, smemSize>>>(
      blocks, npix, d_idata, values, d_odata, scale, offset);

  carmaCheckMsg("centroid_bpix<<<>>> execution failed\n");
}
template void subap_bpcentro<float>(int threads, int blocks, int npix,
                                    float *d_idata, unsigned int *values,
                                    float *d_odata, float scale, float offset);
template void subap_bpcentro<double>(int threads, int blocks, int npix,
                                     double *d_idata, unsigned int *values,
                                     double *d_odata, double scale,
                                     double offset);
