#include "carma_utils.h"
#include "carma_utils.cuh"

template <class T_data>
__global__ void find_nnz_krnl(T_data *d_data, int *d_nnz, int N) {
	int *sdata = SharedMemory<int>();
	int tid = threadIdx.x + blockDim.x*blockIdx.x;
	int sid = threadIdx.x;
  if(blockIdx.x == 0)
    d_nnz[0] = 0;

	//Load shared memory with 1 if d_data[tid]!= 0, with 0 else
	sdata[sid] = (tid<N ?(d_data[tid] != 0) : 0);

	__syncthreads();
	reduce_krnl(sdata,blockDim.x,sid);
	__syncthreads();

	if(threadIdx.x == 0)
//		subsum[blockIdx.x] = sdata[0];
    atomicAdd(d_nnz, sdata[0]);

}
template <class T_data>
int
find_nnz(T_data *d_data, int N, int *d_nnz, int device) {
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	int smemSize = nthreads * sizeof(int);

	find_nnz_krnl<<<grid, threads, smemSize>>>(d_data,d_nnz,N);
	cutilCheckMsg("find_nnz_krnl<<<>>> execution failed\n");
	int nnz = 0;
	cutilSafeCall(
				cudaMemcpy(&nnz,d_nnz,sizeof(int),cudaMemcpyDeviceToHost));
	return nnz;
}

template
int find_nnz<float>(float *d_data, int N, int *d_nnz, int device);
template
int find_nnz<double>(double *d_data, int N, int *d_nnz, int device);
