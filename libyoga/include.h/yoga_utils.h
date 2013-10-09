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
 
#ifndef _YOGA_UTILS_H_
#define _YOGA_UTILS_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cufft.h>

#include <driver_types.h>
#include <vector_types.h>

#include <cuda_runtime_api.h>

#include <cufft.h>

////////////////////////////////////////////////////////////////////////////
//! CUT bool type
////////////////////////////////////////////////////////////////////////////
enum CUTBoolean 
  {
    CUTFalse = 0,
    CUTTrue = 1
  };

// We define these calls here, so the user doesn't need to include __FILE__ and __LINE__
// The advantage is the developers gets to use the inline function so they can debug
#define cutilSafeCallNoSync(err)     __cudaSafeCallNoSync(err, __FILE__, __LINE__)
#define cutilSafeCall(err)           __cudaSafeCall      (err, __FILE__, __LINE__)
#define cutilSafeThreadSync()        __cudaSafeThreadSync(__FILE__, __LINE__)
#define cufftSafeCall(err)           __cufftSafeCall     (err, __FILE__, __LINE__)
#define cutilCheckError(err)         __cutilCheckError   (err, __FILE__, __LINE__)
#define cutilCheckMsg(msg)           __cutilCheckMsg     (msg, __FILE__, __LINE__)
#define cutilSafeMalloc(mallocCall)  __cutilSafeMalloc   ((mallocCall), __FILE__, __LINE__)

#define MIN(a,b) ((a < b) ? a : b)
#define MAX(a,b) ((a > b) ? a : b)

// Give a little more for Windows : the console window often disapears before we can read the message
#ifdef _WIN32
# if 1//ndef UNICODE
#  ifdef _DEBUG // Do this only in debug mode...
	inline void VSPrintf(FILE *file, LPCSTR fmt, ...)
	{
		size_t fmt2_sz	= 2048;
		char *fmt2		= (char*)malloc(fmt2_sz);
		va_list  vlist;
		va_start(vlist, fmt);
		while((_vsnprintf(fmt2, fmt2_sz, fmt, vlist)) < 0) // means there wasn't anough room
		{
			fmt2_sz *= 2;
			if(fmt2) free(fmt2);
			fmt2 = (char*)malloc(fmt2_sz);
		}
		OutputDebugStringA(fmt2);
		fprintf(file, fmt2);
		free(fmt2);
	}
#	define FPRINTF(a) VSPrintf a
#  else //debug
#	define FPRINTF(a) fprintf a
// For other than Win32
#  endif //debug
# else //unicode
// Unicode case... let's give-up for now and keep basic printf
#	define FPRINTF(a) fprintf a
# endif //unicode
#else //win32
#	define FPRINTF(a) fprintf a
#endif //win32

// NOTE: "%s(%i) : " allows Visual Studio to directly jump to the file at the right line
// when the user double clicks on the error line in the Output pane. Like any compile error.

inline void __cudaSafeCallNoSync( cudaError err, const char *file, const int line )
{
    if( cudaSuccess != err) {
        FPRINTF((stderr, "%s(%i) : cudaSafeCallNoSync() Runtime API error : %s.\n",
                file, line, cudaGetErrorString( err) ));
        exit(-1);
    }
}

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
    if( cudaSuccess != err) {
		FPRINTF((stderr, "%s(%i) : cudaSafeCall() Runtime API error : %s.\n",
                file, line, cudaGetErrorString( err) ));
        exit(-1);
    }
}

inline void __cudaSafeThreadSync( const char *file, const int line )
{
    cudaError err = cudaThreadSynchronize();
    if ( cudaSuccess != err) {
        FPRINTF((stderr, "%s(%i) : cudaThreadSynchronize() Driver API error : %s.\n",
                file, line, cudaGetErrorString( err) ));
        exit(-1);
    }
}

inline void __cufftSafeCall( cufftResult err, const char *file, const int line )
{
    if( CUFFT_SUCCESS != err) {
        FPRINTF((stderr, "%s(%i) : cufftSafeCall() CUFFT error.\n",
                file, line));
        exit(-1);
    }
}

inline void __cutilCheckError( CUTBoolean err, const char *file, const int line )
{
    if( CUTTrue != err) {
        FPRINTF((stderr, "%s(%i) : CUTIL CUDA error.\n",
                file, line));
        exit(-1);
    }
}

inline void __cutilCheckMsg( const char *errorMessage, const char *file, const int line )
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) {
        FPRINTF((stderr, "%s(%i) : cutilCheckMsg() CUTIL CUDA error : %s : %s.\n",
                file, line, errorMessage, cudaGetErrorString( err) ));
        exit(-1);
    }
#ifdef _DEBUG
    err = cudaThreadSynchronize();
    if( cudaSuccess != err) {
		FPRINTF((stderr, "%s(%i) : cutilCheckMsg cudaThreadSynchronize error: %s : %s.\n",
                file, line, errorMessage, cudaGetErrorString( err) ));
        exit(-1);
    }
#endif
}
inline void __cutilSafeMalloc( void *pointer, const char *file, const int line )
{
    if( !(pointer)) {
        FPRINTF((stderr, "%s(%i) : cutilSafeMalloc host malloc failure\n",
                file, line));
        exit(-1);
    }
}


#endif // _CUTIL_INLINE_FUNCTIONS_RUNTIME_H_
