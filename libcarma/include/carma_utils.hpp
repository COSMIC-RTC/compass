// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma_utils.hpp
//! \ingroup   libcarma
//! \brief     this file provides tools to CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef _CARMA_UTILS_H_
#define _CARMA_UTILS_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <driver_types.h>
#include <vector_types.h>
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <cuda.h>
#include <cuda_fp16.h>
#include <cuda_runtime_api.h>
#include <cufft.h>

#include <carma_indicators.hpp>

#ifdef USE_OCTOPUS
#include <Cacao.h>
#endif

#define CARMA_PI 3.1415926535897932384626433832

struct doubleint {
  int32_t start;
  int32_t nbInflu;
};

template <class T>
struct tuple_t {
  int32_t pos;
  T data;
};

namespace carma_utils {
template <typename T>
inline std::string to_string(const T &n) {
  std::ostringstream stm;
  stm << n;
  return stm.str();
}
template <typename T>
inline T from_string(const std::string &myString) {
  std::istringstream buffer(myString);
  T value;
  buffer >> value;
  return value;
}
void inline split(std::vector<std::string> &tokens, const std::string &text,
                  char sep) {
  std::string::size_type start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
}

}  // namespace carma_utils

#ifdef DEBUG
#define DEBUG_TRACE(fmt, args...) \
  fprintf(stderr, "[%s@%d]: " fmt "\n", __FILE__, __LINE__, ##args)
#else
#define DEBUG_TRACE(fmt, args...) \
  fprintf(stderr, "[%s@%d]: " fmt "\n", __FILE__, __LINE__, ##args)
#endif

#define CAST(type, new_var, var) type new_var = dynamic_cast<type>(var)
#define SCAST(type, new_var, var) type new_var = static_cast<type>(var)

// We define these calls here, so the user doesn't need to include __FILE__ and
// __LINE__ The advantage is the developers gets to use the inline function so
// they can debug
#ifdef DEBUG
#define carma_safe_call_no_sync(err) __carma_safe_call_no_sync(err, __FILE__, __LINE__)
#define carma_safe_call(err) __carma_safe_call((err), #err, __FILE__, __LINE__)
#define carma_safe_device_synchronize() \
  __carma_safe_device_synchronize(__FILE__, __LINE__)
#define carmafft_safe_call(err) __carmafft_safe_call(err, __FILE__, __LINE__)
#define carma_check_msg(msg) __carma_check_msg(msg, __FILE__, __LINE__)
#define carma_safe_malloc(mallocCall) \
  __carma_safe_malloc((mallocCall), __FILE__, __LINE__)
#else
#define carma_safe_call_no_sync(err) err
#define carma_safe_call(err) err
#define carma_safe_device_synchronize() cudaDeviceSynchronize()
#define carmafft_safe_call(err) err
#define cutil_check_error(err) err
#define carma_check_msg(msg)
#define cutil_safe_malloc(mallocCall) (mallocCall)
#endif

#ifndef MIN
#define MIN(a, b) ((a < b) ? a : b)
#endif
#ifndef MAX
#define MAX(a, b) ((a > b) ? a : b)
#endif

inline uint32_t next_pow2(uint32_t x) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}

inline bool is_pow2(uint32_t x) { return ((x & (x - 1)) == 0); }

class CarmaDevice;
void get_num_blocks_and_threads(CarmaDevice *device, int32_t n, int32_t &blocks,
                            int32_t &threads);
void sum_get_num_blocks_and_threads(int32_t n, CarmaDevice *device, int32_t &blocks,
                               int32_t &threads);
template <class T_data>
int32_t find_nnz(T_data *d_data, int32_t *tmp_colind, int32_t N, int32_t *d_nnz, int32_t &h_nnz,
             CarmaDevice *device);
template <class T_data>
int32_t fill_sparse_vect(T_data *dense_data, int32_t *colind_sorted, T_data *values,
                     int32_t *colind, int32_t *rowind, int32_t nnz, CarmaDevice *device);
int32_t float_to_double(float *idata, double *odata, int32_t N, CarmaDevice *device);
int32_t double_to_float(double *idata, float *odata, int32_t N, CarmaDevice *device);
int32_t print_mem_info();
template <typename T_data>
int32_t fill_array_with_value(T_data *d_data, T_data value, int32_t N,
                          CarmaDevice *device);

void carma_start_profile();
void carma_stop_profile();

// NOTE: "(%s:%i) : " allows Eclipse to directly jump to the file at the right
// line when the user double clicks on the error line in the Output pane. Like
// any compile error.

inline void __carma_safe_call_no_sync(cudaError err, const char *file,
                                  const int32_t line) {
  if (cudaSuccess != err) {
    fprintf(stderr, "(%s:%i) : carma_safe_call_no_sync() Runtime API error : %s.\n",
            file, line, cudaGetErrorString(err));
    // exit(EXIT_FAILURE);
    throw cudaGetErrorString(err);
  }
}

inline void __carma_safe_call(cudaError err, const char *code, const char *file,
                            const int32_t line) {
  if (cudaSuccess != err) {
    fprintf(stderr, "[%s:%i] %s\n carma_safe_call() Runtime API error : %s.\n",
            file, line, code, cudaGetErrorString(err));
    // exit(EXIT_FAILURE);
    throw cudaGetErrorString(err);
  }
}

inline void __carma_safe_device_synchronize(const char *file, const int32_t line) {
  cudaError err = cudaDeviceSynchronize();
  if (cudaSuccess != err) {
    fprintf(stderr,
            "(%s:%i) : cudaDeviceSynchronize() Driver API error : %s.\n", file,
            line, cudaGetErrorString(err));
    // exit(EXIT_FAILURE);
    throw cudaGetErrorString(err);
  }
}

inline void __carmafft_safe_call(cufftResult err, const char *file,
                               const int32_t line) {
  if (CUFFT_SUCCESS != err) {
    fprintf(stderr, "(%s:%i) : carmafft_safe_call() CUFFT error.\n", file, line);
    // exit(EXIT_FAILURE);
    throw "carmafft_safe_call() CUFFT error";
  }
}

inline void __carma_check_msg(const char *error_message, const char *file,
                            const int32_t line) {
  cudaError_t err = cudaGetLastError();
  if (cudaSuccess != err) {
    fprintf(stderr, "(%s:%i) : carma_check_msg() CUTIL CUDA error : %s : %s.\n",
            file, line, error_message, cudaGetErrorString(err));
    throw cudaGetErrorString(err);
  }
#ifdef DEBUG
  err = cudaDeviceSynchronize();
  if (cudaSuccess != err) {
    fprintf(stderr,
            "(%s:%i) : carma_check_msg cudaDeviceSynchronize error: %s : %s.\n",
            file, line, error_message, cudaGetErrorString(err));
    throw cudaGetErrorString(err);
  }
#endif
}
inline void __carma_safe_malloc(void *pointer, const char *file, const int32_t line) {
  if (!(pointer)) {
    fprintf(stderr, "(%s:%i) : cutil_safe_malloc host malloc failure\n", file,
            line);
    throw "cutil_safe_malloc() cutil_safe_malloc host malloc failure";
  }
}

#endif  // _CARMA_UTILS_H_
