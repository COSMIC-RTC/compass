

#include "kp_cu_cuda.h"
#include <iostream>

void kp_cu_cudaMemcpy(void* ptr_src,const void* ptr_dest, size_t nb_octets,enum cudaMemcpyKind direction)
{
	cudaError_t cudaStat;
	cudaStat = cudaMemcpy(ptr_src, ptr_dest, nb_octets, direction);
	if( cudaStat != cudaSuccess )
	{
		cerr << "error | kp_cu_cudaMemcpy | cudaMemcpy failed : "<< cudaGetErrorString(cudaStat)<<endl;
		throw KP_CUDA_MEMCPY_FAILED;
		//exit(EXIT_FAILURE);
	}
}
void kp_cu_cudaMemset(void* ptr_dev, int val, size_t nb_octets)
{
	cudaError_t cudaStat;
	cudaStat = cudaMemset(ptr_dev, val, nb_octets);
	if( cudaStat != cudaSuccess )
	{
		 cerr << "error | kp_cu_cudaMemset | cudaMemset failed : "<< cudaGetErrorString(cudaStat)<<endl;
		 throw KP_CUDA_MEMSET_FAILED;
		 //exit(EXIT_FAILURE);
	}
}
void kp_cu_cudaMalloc(void** ptr, size_t nb_octets)
{
	cudaError_t cudaStat ;
	cudaStat = cudaMalloc(ptr, nb_octets);
	if(cudaStat != cudaSuccess) 
	{
		 cerr << "error | kp_cu_cudaMalloc | cudaMalloc failed : "<< cudaGetErrorString(cudaStat) << endl;
		 throw KP_CUDA_MALLOC_FAILED;
		 //exit(EXIT_FAILURE);
	}
}
void kp_cu_cudaFree(void* ptr_dev)
{
	cudaError_t cudaStat ;
	cudaStat = cudaFree(ptr_dev);
	if(cudaStat != cudaSuccess) 
	{
		 cerr << "error | kp_cu_cudaFree | cudaFree failed : " << cudaGetErrorString(cudaStat) << endl;
		 throw KP_CUDA_FREE_FAILED;
		 //exit(EXIT_FAILURE);
	}
}
