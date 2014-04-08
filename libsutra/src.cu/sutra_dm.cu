#include <sutra_dm.h>

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory {
  __device__
  inline operator T*() {
    extern __shared__ int __smem[];
    return (T*) __smem;
  }

  __device__
  inline operator const T*() const {
    extern __shared__ int __smem[];
    return (T*) __smem;
  }
};

// specialize for double to avoid unaligned memory 
// access compile errors
template<>
struct SharedMemory<double> {
  __device__
  inline operator double*() {
    extern __shared__ double __smem_d[];
    return (double*) __smem_d;
  }

  __device__
  inline operator const double*() const {
    extern __shared__ double __smem_d[];
    return (double*) __smem_d;
  }
};

template<class T>
__global__ void dmshape_krnl(T *g_idata, T *g_odata, int *pos, int *istart,
    int *npts, T *comm, unsigned int n, int N) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    int local_istart = istart[i];
    int local_npts = npts[i];

    sdata[tid] = 0;

    if (local_npts > 0) {
      for (int cc = 0; cc < local_npts; cc++) {
        int lpos = pos[local_istart + cc];
        int ninflu = lpos / n;
        sdata[tid] += comm[ninflu] * g_idata[lpos];
      }
    }
  }

  __syncthreads();

  if (i < N) {
    // write result for this block to global mem
    g_odata[i] = sdata[tid];
  }
}

template<class T>
void comp_dmshape(int threads, int blocks, T *d_idata, T *d_odata, int *pos,
    int *istart, int *npts, T *comm, unsigned int n, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  dmshape_krnl<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, pos,
      istart, npts, comm, n, N);

  cutilCheckMsg("dmshape_kernel<<<>>> execution failed\n");

}

template void
comp_dmshape<float>(int threads, int blocks, float *d_idata, float *d_odata,
    int *pos, int *istart, int *npts, float *comm, unsigned int n, int N);

template void
comp_dmshape<double>(int threads, int blocks, double *d_idata, double *d_odata,
    int *pos, int *istart, int *npts, double *comm, unsigned int n, int N);

template<class T>
__global__ void oneactu_krnl(T *g_idata, T *g_odata, int nactu, T ampli,
    int *xoff, int *yoff, int dim_im, int dim_influ, int N) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    int iy = i / dim_im;
    int ix = i - iy * dim_im;
    int ixactu = ix - xoff[nactu];
    int iyactu = iy - yoff[nactu];

    // write result for this block to global mem
    if ((ixactu > -1) && (ixactu < dim_influ) && (iyactu > -1)
        && (iyactu < dim_influ)) {
      int tid = ixactu + iyactu * dim_influ + nactu * dim_influ * dim_influ;
      g_odata[i] = ampli * g_idata[tid];
    }
  }
}

template<class T>
__global__ void oneactu_krnl_fast(T *g_idata, T *g_odata, int nactu, T ampli,
    int *xoff, int *yoff, int dim_im, int dim_influ, int N) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    int iy = i / dim_influ;
    int ix = i - iy * dim_influ;
    int ixactu = ix + xoff[nactu];
    int iyactu = iy + yoff[nactu];

    // write result for this block to global mem
    if ((ixactu > -1) && (ixactu < dim_im) && (iyactu > -1)
        && (iyactu < dim_im)) {
      int tid = ixactu + iyactu * dim_im;
      g_odata[tid] = ampli * g_idata[i + nactu * dim_influ * dim_influ];
    }
  }
}

template<class T>
void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu,
    T ampli, int *xoff, int *yoff, int dim_im, int dim_influ, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  //oneactu_krnl<T><<< dimGrid, dimBlock >>>(d_idata,d_odata, nactu, ampli, xoff, yoff, dim_im, dim_influ, N);
  oneactu_krnl_fast<T> <<<dimGrid, dimBlock>>>(d_idata, d_odata, nactu, ampli,
      xoff, yoff, dim_im, dim_influ, N);

  cutilCheckMsg("oneactu_kernel<<<>>> execution failed\n");

}

template void
oneactu<float>(int threads, int blocks, float *d_idata, float *d_odata,
    int nactu, float ampli, int *xoff, int *yoff, int dim_im, int dim_influ,
    int N);

template void
oneactu<double>(int threads, int blocks, double *d_idata, double *d_odata,
    int nactu, double ampli, int *xoff, int *yoff, int dim_im, int dim_influ,
    int N);

template<class T>
__global__ void oneactu_krnl_fast(T *g_idata, T *g_odata, int nactu, T ampli,
    int dim_im, int dim_influ, int N) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    int iy = i / dim_influ;
    int ix = i - iy * dim_influ;

    // write result for this block to global mem
    if ((ix > -1) && (ix < dim_im) && (iy > -1) && (iy < dim_im)) {
      int tid = ix + iy * dim_im;
      g_odata[tid] = ampli * g_idata[i + nactu * dim_influ * dim_influ];
    }
  }
}

template<class T>
void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu,
    T ampli, int dim_im, int dim_influ, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  //oneactu_krnl<T><<< dimGrid, dimBlock >>>(d_idata,d_odata, nactu, ampli, xoff, yoff, dim_im, dim_influ, N);
  oneactu_krnl_fast<T> <<<dimGrid, dimBlock>>>(d_idata, d_odata, nactu, ampli,
      dim_im, dim_influ, N);

  cutilCheckMsg("oneactu_kernel<<<>>> execution failed\n");

}

template void
oneactu<float>(int threads, int blocks, float *d_idata, float *d_odata,
    int nactu, float ampli, int dim_im, int dim_influ, int N);

template void
oneactu<double>(int threads, int blocks, double *d_idata, double *d_odata,
    int nactu, double ampli, int dim_im, int dim_influ, int N);

template<class T>
__global__ void fulldmshape_krnl(T *g_idata, T *g_odata, int ninflu,
    int diminflu, T *comm, int N) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    sdata[tid] = 0;

    for (int cc = 0; cc < ninflu; cc++) {
      sdata[tid] += comm[cc] * g_idata[i + cc * diminflu];
    }
  }
  __syncthreads();

  if (i < N) {
    // write result for this block to global mem
    g_odata[i] = sdata[tid];
  }
}

template<class T>
void comp_fulldmshape(int threads, int blocks, T *d_idata, T *d_odata,
    int ninflu, int diminflu, T *comm, int N) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  fulldmshape_krnl<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata,
      ninflu, diminflu, comm, N);

  cutilCheckMsg("fulldmshape_kernel<<<>>> execution failed\n");

}

template void
comp_fulldmshape<float>(int threads, int blocks, float *d_idata, float *d_odata,
    int ninflu, int diminflu, float *comm, int N);

template void
comp_fulldmshape<double>(int threads, int blocks, double *d_idata,
    double *d_odata, int ninflu, int diminflu, double *comm, int N);

