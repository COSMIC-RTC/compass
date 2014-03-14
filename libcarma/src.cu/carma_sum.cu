#include <carma_obj.h>

bool
isPow2(unsigned int x) {
  return ((x & (x - 1)) == 0);
}

unsigned int
nextPow2(unsigned int x) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}

void
sumGetNumBlocksAndThreads(int n, int device, int &blocks, int &threads) {

  struct cudaDeviceProp deviceProperties;

  cudaGetDeviceProperties(&deviceProperties, device);

  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int maxBlocks = deviceProperties.multiProcessorCount;

  threads = (n < maxThreads * 2) ? nextPow2((n + 1) / 2) : maxThreads;
  blocks = (n + (threads * 2 - 1)) / (threads * 2);
  blocks = MIN(maxBlocks, blocks);

}

template<class T>
  struct SharedMemory {
    __device__
    inline
    operator T*() {
      extern __shared__ int __smem[];
      return (T*) __smem;
    }

    __device__
    inline
    operator const T*() const {
      extern __shared__ int __smem[];
      return (T*) __smem;
    }
  };

// specialize for double to avoid unaligned memory 
// access compile errors
template<>
  struct SharedMemory<double> {
    __device__
    inline
    operator double*() {
      extern __shared__ double __smem_d[];
      return (double*) __smem_d;
    }

    __device__
    inline
    operator const double*() const {
      extern __shared__ double __smem_d[];
      return (double*) __smem_d;
    }
  };
/*
 Parallel sum reduction using shared memory
 - takes log(n) steps for n input elements
 - uses n threads
 - only works for power-of-2 arrays
 */

/*
 This version adds multiple elements per thread sequentially.  This reduces the overall
 cost of the algorithm while keeping the work complexity O(n) and the step complexity 
 O(log n).
 (Brent's Theorem optimization)

 Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory. 
 In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.  
 If blockSize > 32, allocate blockSize*sizeof(T) bytes.
 */
// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T, unsigned int blockSize, bool nIsPow2>
  __global__ void
  reduce6(T *g_idata, T *g_odata, unsigned int n) {
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockSize * 2 + threadIdx.x;
    unsigned int gridSize = blockSize * 2 * gridDim.x;

    T mySum = 0;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n) {
      mySum += g_idata[i];
      // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
      if (nIsPow2 || i + blockSize < n)
        mySum += g_idata[i + blockSize];
      i += gridSize;
    }

    // each thread puts its local sum into shared memory 
    sdata[tid] = mySum;
    __syncthreads();

    // do reduction in shared mem
    if (blockSize >= 1024) {
      if (tid < 512) {
        sdata[tid] = mySum = mySum + sdata[tid + 512];
      }
      __syncthreads();
    }
    if (blockSize >= 512) {
      if (tid < 256) {
        sdata[tid] = mySum = mySum + sdata[tid + 256];
      }
      __syncthreads();
    }
    if (blockSize >= 256) {
      if (tid < 128) {
        sdata[tid] = mySum = mySum + sdata[tid + 128];
      }
      __syncthreads();
    }
    if (blockSize >= 128) {
      if (tid < 64) {
        sdata[tid] = mySum = mySum + sdata[tid + 64];
      }
      __syncthreads();
    }

#ifndef __DEVICE_EMULATION__
    if (tid < 32)
#endif
        {
      // now that we are using warp-synchronous programming (below)
      // we need to declare our shared memory volatile so that the compiler
      // doesn't reorder stores to it and induce incorrect behavior.
      volatile T* smem = sdata;
      if (blockSize >= 64) {
        smem[tid] = mySum = mySum + smem[tid + 32];
        __syncthreads();
      }
      if (blockSize >= 32) {
        smem[tid] = mySum = mySum + smem[tid + 16];
        __syncthreads();
      }
      if (blockSize >= 16) {
        smem[tid] = mySum = mySum + smem[tid + 8];
        __syncthreads();
      }
      if (blockSize >= 8) {
        smem[tid] = mySum = mySum + smem[tid + 4];
        __syncthreads();
      }
      if (blockSize >= 4) {
        smem[tid] = mySum = mySum + smem[tid + 2];
        __syncthreads();
      }
      if (blockSize >= 2) {
        smem[tid] = mySum = mySum + smem[tid + 1];
        __syncthreads();
      }
    }

    // write result for this block to global mem 
    if (tid == 0)
      g_odata[blockIdx.x] = sdata[0];
  }

template<class T>
  void
  reduce(int size, int threads, int blocks, T *d_idata, T *d_odata) {
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize =
        (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

    if (isPow2(size)) {
      switch (threads) {
      case 1024:
        reduce6<T, 1024, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 512:
        reduce6<T, 512, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 256:
        reduce6<T, 256, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 128:
        reduce6<T, 128, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 64:
        reduce6<T, 64, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      case 32:
        reduce6<T, 32, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      case 16:
        reduce6<T, 16, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      case 8:
        reduce6<T, 8, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      case 4:
        reduce6<T, 4, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      case 2:
        reduce6<T, 2, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      case 1:
        reduce6<T, 1, true> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      }
    } else {
      switch (threads) {
      case 1024:
        reduce6<T, 1024, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 512:
        reduce6<T, 512, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 256:
        reduce6<T, 256, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 128:
        reduce6<T, 128, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 64:
        reduce6<T, 64, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 32:
        reduce6<T, 32, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 16:
        reduce6<T, 16, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
            d_odata, size);
        break;
      case 8:
        reduce6<T, 8, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      case 4:
        reduce6<T, 4, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      case 2:
        reduce6<T, 2, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      case 1:
        reduce6<T, 1, false> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
            size);
        break;
      }
    }

  }

template void
reduce<int>(int size, int threads, int blocks, int *d_idata, int *d_odata);

template void
reduce<float>(int size, int threads, int blocks, float *d_idata,
    float *d_odata);

template void
reduce<double>(int size, int threads, int blocks, double *d_idata,
    double *d_odata);
