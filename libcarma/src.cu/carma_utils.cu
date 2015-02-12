#include "carma_utils.h"
#include "carma_utils.cuh"
#include <thrust/sort.h>
#include <thrust/device_ptr.h>

template <class T_data>
__global__ void find_nnz_krnl(T_data *d_data, int *colind, int *d_nnz, int N) {
	int *sdata = SharedMemory<int>();
	int tid = threadIdx.x + blockDim.x*blockIdx.x;
	int sid = threadIdx.x;
  if(tid == 0)
    d_nnz[0] = 0;

	//Load shared memory with 1 if d_data[tid]!= 0, with 0 else
	if(tid<N) {
		sdata[sid] = (d_data[tid] != 0);
		colind[tid] = (sdata[sid]) ? tid : N+tid; // Init colind for further sort
	} else {
		sdata[sid] = 0;
	}
	__syncthreads();
	reduce_krnl(sdata,blockDim.x,sid);
	__syncthreads();

	if(threadIdx.x == 0)
//		subsum[blockIdx.x] = sdata[0];
		atomicAdd(d_nnz, sdata[0]);

}
template <class T_data>
int
find_nnz(T_data *d_data, int *colind, int N, int *d_nnz, int &h_nnz, carma_device *device) {
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	int smemSize = nthreads * sizeof(int);

	find_nnz_krnl<<<grid, threads, smemSize>>>(d_data,colind,d_nnz,N);
	cutilCheckMsg("find_nnz_krnl<<<>>> execution failed\n");

	// wrap raw pointer with a device_ptr
	thrust::device_ptr<int> dev_ptr(colind);

	thrust::sort(dev_ptr,dev_ptr+N);
	cutilSafeCall(
				cudaMemcpy(&h_nnz,d_nnz,sizeof(int),cudaMemcpyDeviceToHost));

	return EXIT_SUCCESS;
}

template
int find_nnz<float>(float *d_data, int *colind,int N, int *d_nnz, int &h_nnz, carma_device *device);
template
int find_nnz<double>(double *d_data, int *colind, int N, int *d_nnz, int &h_nnz, carma_device *device);



template <class T_data>
__global__ void fill_sparse_vect_krnl(T_data *dense_data, int *colind_sorted, T_data *values, int *colind, int *rowind, int nnz) {
	int tid = threadIdx.x + blockDim.x*blockIdx.x;
	if(tid == 0)
		rowind[0] = 0;
	if(tid == 1)
		rowind[1] = nnz;

	//Load shared memory with 1 if d_data[tid]!= 0, with 0 else
	if(tid<nnz) {
		values[tid] = dense_data[colind_sorted[tid]];
		colind[tid] = colind_sorted[tid];
	}
	__syncthreads();
}
template <class T_data>
int
fill_sparse_vect(T_data *dense_data, int *colind_sorted, T_data *values, int *colind, int *rowind, int nnz, carma_device *device) {
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, nnz, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);

	fill_sparse_vect_krnl<<<grid, threads>>>(dense_data, colind_sorted, values, colind, rowind, nnz);
	cutilCheckMsg("find_nnz_krnl<<<>>> execution failed\n");
	return EXIT_SUCCESS;
}
template
int
fill_sparse_vect<float>(float *dense_data, int *colind_sorted, float *values, int *colind, int *rowind, int nnz, carma_device *device);
template
int
fill_sparse_vect<double>(double *dense_data, int *colind_sorted, double *values, int *colind, int *rowind, int nnz, carma_device *device);

__global__ void
floattodouble_krnl(float *i_data, double *o_data, int N){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < N){
		o_data[tid] = (double)i_data[tid];
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void
doubletofloat_krnl(double *i_data, float *o_data, int N){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < N){
		o_data[tid] = (float)i_data[tid];
		tid += blockDim.x * gridDim.x;
	}
}

int
floattodouble(float *i_data, double *o_data, int N, carma_device *device){
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	floattodouble_krnl<<<grid , threads>>>(i_data,o_data,N);
	cutilCheckMsg("floattodouble_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}

int
doubletofloat(double *i_data, float *o_data, int N, carma_device *device){
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	doubletofloat_krnl<<<grid , threads>>>(i_data,o_data,N);
	cutilCheckMsg("floattodouble_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}
