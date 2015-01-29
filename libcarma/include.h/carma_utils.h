/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 * 
 */

#ifndef _CARMA_UTILS_H_
#define _CARMA_UTILS_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cufft.h>

#include <driver_types.h>
#include <vector_types.h>

#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cufft.h>

#ifdef DEBUG
#define DEBUG_TRACE(fmt, args...) fprintf(stderr, "[%s@%d]: " fmt "\n", __FILE__, __LINE__, ## args)
#else
#define DEBUG_TRACE(fmt, args...) fprintf(stderr, "[%s@%d]: " fmt "\n", __FILE__, __LINE__, ## args)
#endif

#define SCAST(type, new_var, var) type new_var=static_cast<type>(var)

////////////////////////////////////////////////////////////////////////////
//! CUT bool type
////////////////////////////////////////////////////////////////////////////
enum CUTBoolean {
  CUTFalse = 0, CUTTrue = 1
};

// We define these calls here, so the user doesn't need to include __FILE__ and __LINE__
// The advantage is the developers gets to use the inline function so they can debug
//#ifdef DEBUG
#define CUSafeCall(err)              __CUSafeCall        (err, __FILE__, __LINE__)
#define cutilSafeCallNoSync(err)     __cudaSafeCallNoSync(err, __FILE__, __LINE__)
#define cutilSafeCall(err)           __cudaSafeCall      (err, __FILE__, __LINE__)
#define cutilSafeThreadSync()        __cudaSafeThreadSync(__FILE__, __LINE__)
#define cufftSafeCall(err)           __cufftSafeCall     (err, __FILE__, __LINE__)
#define cutilCheckError(err)         __cutilCheckError   (err, __FILE__, __LINE__)
#define cutilCheckMsg(msg)           __cutilCheckMsg     (msg, __FILE__, __LINE__)
#define cutilSafeMalloc(mallocCall)  __cutilSafeMalloc   ((mallocCall), __FILE__, __LINE__)
//#else
//#define CUSafeCall(err)              err
//#define cutilSafeCallNoSync(err)     err
//#define cutilSafeCall(err)           err
//#define cutilSafeThreadSync()        cudaThreadSynchronize()
//#define cufftSafeCall(err)           err
//#define cutilCheckError(err)         err
//#define cutilCheckMsg(msg)
//#define cutilSafeMalloc(mallocCall)  (mallocCall)
//#endif

#define MIN(a,b) ((a < b) ? a : b)
#define MAX(a,b) ((a > b) ? a : b)

inline unsigned int nextPow2(unsigned int x) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}

inline bool isPow2(unsigned int x) {
  return ((x & (x - 1)) == 0);
}

class carma_device;
void getNumBlocksAndThreads(carma_device *device, int n, int &blocks, int &threads);
template <class T_data>
int  find_nnz(T_data *d_data, int *tmp_colind, int N, int *d_nnz, int &h_nnz, carma_device *device);
template <class T_data>
int
fill_sparse_vect(T_data *dense_data, int *colind_sorted, T_data *values, int *colind, int *rowind, int nnz, carma_device *device);

void carma_start_profile();
void carma_stop_profile();

// NOTE: "(%s:%i) : " allows Eclipse to directly jump to the file at the right line
// when the user double clicks on the error line in the Output pane. Like any compile error.

inline void __cudaSafeCallNoSync(cudaError err, const char *file,
    const int line) {
  if (cudaSuccess != err) {
    fprintf(stderr, "(%s:%i) : cudaSafeCallNoSync() Runtime API error : %s.\n",
        file, line, cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

inline void __cudaSafeCall(cudaError err, const char *file, const int line) {
  if (cudaSuccess != err) {
    fprintf(stderr, "(%s:%i) : cudaSafeCall() Runtime API error : %s.\n", file,
        line, cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

inline void __CUSafeCall(CUresult err, const char *file, const int line) {
  if (CUDA_SUCCESS != err) {
    /* char *text; */
    /* if(cuGetErrorName(err, &text) == CUDA_ERROR_INVALID_VALUE )  */
      fprintf(stderr, "(%s:%i) : CUSafeCall() Driver API not recognized error : %d.\n", file,
	      line, err);
    /* else */
    /*   fprintf(stderr, "(%s:%i) : CUSafeCall() Driver API error : %s (%d).\n", file, */
    /* 	      line, text, err); */
    exit(EXIT_FAILURE);
  }
}

inline void __cudaSafeThreadSync(const char *file, const int line) {
  cudaError err = cudaThreadSynchronize();
  if (cudaSuccess != err) {
    fprintf(stderr,
        "(%s:%i) : cudaThreadSynchronize() Driver API error : %s.\n", file,
        line, cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

inline void __cufftSafeCall(cufftResult err, const char *file, const int line) {
  if (CUFFT_SUCCESS != err) {
    fprintf(stderr, "(%s:%i) : cufftSafeCall() CUFFT error.\n", file, line);
    exit(EXIT_FAILURE);
  }
}

inline void __cutilCheckError(CUTBoolean err, const char *file,
    const int line) {
  if (CUTTrue != err) {
    fprintf(stderr, "(%s:%i) : CUTIL CUDA error.\n", file, line);
    exit(EXIT_FAILURE);
  }
}

inline void __cutilCheckMsg(const char *errorMessage, const char *file,
    const int line) {
  cudaError_t err = cudaGetLastError();
  if (cudaSuccess != err) {
    fprintf(stderr, "(%s:%i) : cutilCheckMsg() CUTIL CUDA error : %s : %s.\n",
        file, line, errorMessage, cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
#ifdef DEBUG
  err = cudaThreadSynchronize();
  if( cudaSuccess != err) {
    fprintf(stderr, "(%s:%i) : cutilCheckMsg cudaThreadSynchronize error: %s : %s.\n",
        file, line, errorMessage, cudaGetErrorString( err) );
    exit(EXIT_FAILURE);
  }
#endif
}
inline void __cutilSafeMalloc(void *pointer, const char *file, const int line) {
  if (!(pointer)) {
    fprintf(stderr, "(%s:%i) : cutilSafeMalloc host malloc failure\n", file,
        line);
    exit(EXIT_FAILURE);
  }
}


#endif // _CARMA_UTILS_H_
