#ifndef __KP_CU_CUDA_H__
#define __KP_CU_CUDA_H__

#include "cublas_v2.h"
#include "kp_exceptions.h"
#include <iostream>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>

using namespace std;

void kp_cu_cudaMemcpy(void* ptr_src,const void* ptr_dest, size_t nb_octets,enum cudaMemcpyKind direction);

void kp_cu_cudaMemset(void* ptr_dev, int val, size_t nb_octets);

void kp_cu_cudaMalloc(void** ptr, size_t nb_octets);

void kp_cu_cudaFree(void* ptr_dev);


#endif
