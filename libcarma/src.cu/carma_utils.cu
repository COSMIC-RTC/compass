#include "carma_utils.h"

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
__inline__ __device__ void reduce_krnl( T *sdata, int size, int n) {
  if (size & (size - 1)) { // test if size is not a power of 2
    unsigned int s;
    if ((size & 1) != 0)
      s = size / 2 + 1; //(size&1)==size%2
    else
      s = size / 2;
    unsigned int s_old = size;
    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        sdata[n] += sdata[n + s];
      }
      //__threadfence_block();
      //__threadfence();
      __syncthreads();
     s_old = s;
      s /= 2;
      if ((2 * s < s_old) && (s != 0))
        s += 1;
    }
  } else {
    // do reduction in shared mem
    for (unsigned int s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        sdata[n] += sdata[n + s];
      }
      //__threadfence_block();
      //__threadfence();
      __syncthreads();
    }
  }
}

template <class T_data>
__global__ void find_nnz_krnl(T_data *d_data, int *subsum, int N) {
	int *sdata = SharedMemory<int>();
	int tid = threadIdx.x + blockDim.x*blockIdx.x;
	int sid = threadIdx.x;

	//Load shared memory with 1 if d_data[tid]!= 0, with 0 else
	sdata[sid] = (tid<N ?(d_data[tid] != 0) : 0);

	__syncthreads();
	reduce_krnl(sdata,blockDim.x,sid);
	__syncthreads();

	if(threadIdx.x == 0)
		subsum[blockIdx.x] = sdata[0];

}
template <class T_data>
int
find_nnz(T_data *d_data, int N, int *d_nnz, int device) {
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	int smemSize = nthreads * sizeof(int);

	int *subsum;
	int subsum_c[nblocks];
	cutilSafeCall(
		      cudaMalloc((void** )&(subsum), sizeof(int) * nblocks));
	find_nnz_krnl<<<grid, threads, smemSize>>>(d_data,subsum,N);
	cutilCheckMsg("find_nnz_krnl<<<>>> execution failed\n");
	cutilSafeCall(
				cudaMemcpy(subsum_c,subsum,nblocks*sizeof(int),cudaMemcpyDeviceToHost));
	cutilSafeCall(
				cudaFree(subsum));
	int nnz = 0;
	for (int i=0 ; i<nblocks ; i++)
		nnz += subsum_c[i];
	cutilSafeCall(
			cudaMemcpy(d_nnz,&nnz,sizeof(int),cudaMemcpyHostToDevice));
	return nnz;
}

template
int find_nnz<float>(float *d_data, int N, int *d_nnz, int device);
template
int find_nnz<double>(double *d_data, int N, int *d_nnz, int device);
