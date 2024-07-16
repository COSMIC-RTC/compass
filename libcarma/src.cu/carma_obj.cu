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

//! \file      carma_obj.cu
//! \ingroup   libcarma
//! \brief     this file provides CarmaObj CUDA kernels
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <carma_obj.hpp>
#include <carma_utils.cuh>

/*
 short	2 bytes
 int32_t	4 bytes
 int64_t	4 bytes
 float	4 bytes
 double	8 bytes
 */
// for shmem : 1 float = 4 bytes => shmem can contain
// nb_elem = deviceProperties.sharedMemPerBlock/4 floats
// seems not to work on my laptop ... maybe my x server is already asking a lot
// ... worth a try on blast asap on my laptop, have to divide by 4 ...
//
//
template <class T>
__device__ T carma_sin(T data);
template <>
__device__ float carma_sin(float data) {
  return sinf(data);
}

template <>
__device__ double carma_sin(double data) {
  return sin(data);
}

template <class T>
__global__ void generic1d(T *odata, T *idata, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = carma_sin(2.0f * idata[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t launch_generic1d(T *d_odata, T *d_idata, int32_t N, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  generic1d<<<grid, threads>>>(d_odata, d_idata, N);

  return EXIT_SUCCESS;
}

template int32_t launch_generic1d<float>(float *d_odata, float *d_idata, int32_t N,
                                     CarmaDevice *device);

template int32_t launch_generic1d<double>(double *d_odata, double *d_idata, int32_t N,
                                      CarmaDevice *device);

template <class T>
__global__ void generic2d(T *odata, T *idata, int32_t N) {
  __shared__ T cache[BLOCK_SZ][BLOCK_SZ];

  int32_t x = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t y = threadIdx.y + blockIdx.y * blockDim.y;
  int32_t tid = x + y * blockDim.x * gridDim.x;

  cache[BLOCK_SZ - 1 - threadIdx.x][BLOCK_SZ - 1 - threadIdx.y] =
      carma_sin(2.0f * idata[tid]);

  __syncthreads();

  odata[tid] = cache[BLOCK_SZ - 1 - threadIdx.x][BLOCK_SZ - 1 - threadIdx.y];
}

template <class T>
int32_t launch_generic2d(T *d_odata, T *d_idata, int32_t N1, int32_t N2) {
  dim3 blocks(N1 / BLOCK_SZ, N2 / BLOCK_SZ), threads(BLOCK_SZ, BLOCK_SZ);
  int32_t N = N1 * N2;

  generic2d<<<blocks, threads>>>(d_odata, d_idata, N);

  return EXIT_SUCCESS;
}

template int32_t launch_generic2d<float>(float *d_odata, float *d_idata, int32_t N1,
                                     int32_t N2);
template int32_t launch_generic2d<double>(double *d_odata, double *d_idata, int32_t N1,
                                      int32_t N2);

template <class T_data>
__global__ void krnl_clip(T_data *data, T_data min, T_data max, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    data[tid] = carma_clip<T_data>(data[tid], min, max);
    tid += blockDim.x * gridDim.x;
  }
}

template <class T_data>
void clip_array(T_data *d_data, T_data min, T_data max, int32_t N,
                CarmaDevice *device, cudaStream_t stream) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  krnl_clip<<<grid, threads, 0, stream>>>(d_data, min, max, N);
  carma_check_msg("krnl_clip<<<>>> execution failed\n");
}

template void clip_array<int32_t>(int32_t *d_data, int32_t min, int32_t max, int32_t N,
                              CarmaDevice *device, cudaStream_t stream);
template void clip_array<uint32_t>(uint32_t *d_data, uint32_t min,
                                       uint32_t max, int32_t N,
                                       CarmaDevice *device, cudaStream_t stream);
template void clip_array<float>(float *d_data, float min, float max, int32_t N,
                                CarmaDevice *device, cudaStream_t stream);
template void clip_array<uint16_t>(uint16_t *d_data, uint16_t min, uint16_t max,
                                   int32_t N, CarmaDevice *device, cudaStream_t stream);

template void clip_array<double>(double *d_data, double min, double max, int32_t N,
                                 CarmaDevice *device, cudaStream_t stream);

template <>
void clip_array(cuFloatComplex *d_data, cuFloatComplex min, cuFloatComplex max,
                int32_t N, CarmaDevice *device, cudaStream_t stream) {
  throw "not implemented";
}
template <>
void clip_array(cuDoubleComplex *d_data, cuDoubleComplex min,
                cuDoubleComplex max, int32_t N, CarmaDevice *device, cudaStream_t stream) {
  throw "not implemented";
}
// template <>
// void clip_array(tuple_t<float> *d_data, tuple_t<float> min, tuple_t<float>
// max,
//                 int32_t N, CarmaDevice *device) {
//   throw "not implemented";
// }

template <class T>
__global__ void krnl_fillindex(T *odata, T *idata, int32_t *indx, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = idata[indx[tid]];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t fillindex(T *d_odata, T *d_idata, int32_t *indx, int32_t N, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  krnl_fillindex<<<grid, threads>>>(d_odata, d_idata, indx, N);
  carma_check_msg("krnl_fillindex<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t fillindex<float>(float *d_odata, float *d_idata, int32_t *indx, int32_t N,
                              CarmaDevice *device);

template int32_t fillindex<double>(double *d_odata, double *d_idata, int32_t *indx,
                               int32_t N, CarmaDevice *device);

template <class T>
__global__ void krnl_fillvalues(T *odata, T *val, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = val[0] / N;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t fillvalues(T *d_odata, T *val, int32_t N, CarmaDevice *device, cudaStream_t stream) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  krnl_fillvalues<<<grid, threads, 0, stream>>>(d_odata, val, N);
  carma_check_msg("krnl_fillvalues<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t fillvalues<float>(float *d_odata, float *val, int32_t N,
                               CarmaDevice *device, cudaStream_t stream);

template int32_t fillvalues<double>(double *d_odata, double *val, int32_t N,
                                CarmaDevice *device, cudaStream_t stream);

template int32_t fillvalues<uint32_t>(uint32_t *d_odata, uint32_t *val,
                                      int32_t N, CarmaDevice *device, cudaStream_t stream);

template <class T>
__global__ void getarray2d_krnl(T *odata, T *idata, int32_t tidx0, int32_t Ncol, int32_t NC,
                                int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t tidB;

  while (tid < N) {
    if (Ncol > 1)
      tidB = tidx0 + (tid / Ncol) * NC + (tid % Ncol);
    else
      tidB = tidx0 + tid * NC;
    odata[tid] = idata[tidB];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t getarray2d(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
               CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  getarray2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N);

  carma_check_msg("getarray2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t getarray2d<float>(float *d_odata, float *d_idata, int32_t x0, int32_t Ncol,
                               int32_t NC, int32_t N, CarmaDevice *device);

template int32_t getarray2d<double>(double *d_odata, double *d_idata, int32_t x0,
                                int32_t Ncol, int32_t NC, int32_t N, CarmaDevice *device);

template <class T>
__global__ void fillarray2d_krnl(T *odata, T *idata, int32_t tidx0, int32_t Ncol,
                                 int32_t NC, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t tidB;

  while (tid < N) {
    if (Ncol > 1)
      tidB = tidx0 + (tid / Ncol) * NC + (tid % Ncol);
    else
      tidB = tidx0 + tid * NC;
    odata[tidB] = idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t fillarray2d(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
                CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  fillarray2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N);

  carma_check_msg("fillarray2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t fillarray2d<float>(float *d_odata, float *d_idata, int32_t x0,
                                int32_t Ncol, int32_t NC, int32_t N, CarmaDevice *device);

template int32_t fillarray2d<double>(double *d_odata, double *d_idata, int32_t x0,
                                 int32_t Ncol, int32_t NC, int32_t N, CarmaDevice *device);

template <class T>
__global__ void plus_krnl(T *idata, T alpha, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    idata[tid] += alpha;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t carma_plus(T *d_odata, T alpha, int32_t N, CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  plus_krnl<<<grid, threads>>>(d_odata, alpha, N);

  carma_check_msg("plus_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t carma_plus<float>(float *d_odata, float alpha, int32_t N,
                               CarmaDevice *device);

template int32_t carma_plus<double>(double *d_odata, double alpha, int32_t N,
                                CarmaDevice *device);

template <class T>
__global__ void plusai_krnl(T *odata, T *idata, int32_t i, int32_t sgn, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    if (sgn == 1)
      odata[tid] += idata[i];
    else
      odata[tid] -= idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t carma_plusai(T *d_odata, T *i_data, int32_t i, int32_t sgn, int32_t N,
                 CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  plusai_krnl<<<grid, threads>>>(d_odata, i_data, i, sgn, N);

  carma_check_msg("plusai_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t carma_plusai<float>(float *d_odata, float *d_idata, int32_t i, int32_t sgn,
                                 int32_t N, CarmaDevice *device);

template int32_t carma_plusai<double>(double *d_odata, double *d_idata, int32_t i,
                                  int32_t sgn, int32_t N, CarmaDevice *device);

template <class T>
__global__ void fillarray2d2_krnl(T *odata, T *idata, int32_t tidx0, int32_t Ncol,
                                  int32_t NC, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t tidB;

  while (tid < N) {
    if (Ncol > 1)
      tidB = tidx0 + (tid / Ncol) * NC + (tid % Ncol);
    else
      tidB = tidx0 + tid * NC;
    odata[tidB] = idata[N - tid - 1];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t fillarray2d2(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
                 CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  fillarray2d2_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N);

  carma_check_msg("fillarray2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t fillarray2d2<float>(float *d_odata, float *d_idata, int32_t x0,
                                 int32_t Ncol, int32_t NC, int32_t N, CarmaDevice *device);

template int32_t fillarray2d2<double>(double *d_odata, double *d_idata, int32_t x0,
                                  int32_t Ncol, int32_t NC, int32_t N,
                                  CarmaDevice *device);

template <class T>
__global__ void kern_fill_sym_matrix(char src_uplo, T *data, int32_t Ncol, int32_t N) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;
  int32_t tidB;

  while (tid < N) {
    int32_t row = tid / Ncol;
    int32_t col = tid - row * Ncol;

    if (col > row) {
      tidB = row + col * Ncol;
      if (src_uplo == 'U')
        data[tid] = data[tidB];
      else if (src_uplo == 'L')
        data[tidB] = data[tid];
    }
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int32_t fill_sym_matrix(char src_uplo, T *d_data, int32_t Ncol, int32_t N,
                    CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  kern_fill_sym_matrix<<<grid, threads>>>(src_uplo, d_data, Ncol, N);

  carma_check_msg("kern_fill_sym_matrix<<<>>> execution failed");

  return EXIT_SUCCESS;
}

template int32_t fill_sym_matrix<float>(char uplo, float *d_data, int32_t Ncol, int32_t N,
                                    CarmaDevice *device);

template int32_t fill_sym_matrix<double>(char uplo, double *d_data, int32_t Ncol, int32_t N,
                                     CarmaDevice *device);

/**
 * @brief Kernel to extract a part of the image centred on center_pos
 *
 * @tparam T type of the image items
 * @param d_smallimg extracted small image of size extract_size*extract_size
 * @param d_fullimg full image of size fullimg_size*fullimg_size
 * @param fullimg_size size of the d_fullimg leading dimension
 * @param center_pos position of the center of d_smallimg in d_fullimg
 * @param extract_size size of the d_smallimg leading dimension
 * @param roll get pixels as if d_fullimg need to be roll
 */
template <class T>
__global__ void extract_krnl(T *d_smallimg, const T *d_fullimg,
                             int32_t fullimg_size, int32_t center_pos, int32_t extract_size,
                             bool roll) {
  int32_t tid = threadIdx.x + blockIdx.x * blockDim.x;

  int32_t max_y = center_pos / fullimg_size;
  int32_t max_x = center_pos % fullimg_size;

  int32_t j = tid / extract_size - extract_size / 2;
  int32_t i = tid % extract_size - extract_size / 2;

  int32_t y = max_y + j;
  int32_t x = max_x + i;

  bool inside = true;

  if (y < 0) {
    y += fullimg_size;
    inside = false;
  }
  if (y > fullimg_size - 1) {
    y -= fullimg_size;
    inside = false;
  }
  if (x < 0) {
    x += fullimg_size;
    inside = false;
  }
  if (x > fullimg_size - 1) {
    x -= fullimg_size;
    inside = false;
  }

  d_smallimg[tid] = (inside || roll) ? d_fullimg[x + y * fullimg_size] : 0;
}

/**
 * @brief Extract a part of the image centred on center_pos
 *
 * @tparam T type of the image items
 * @param d_smallimg extracted small image of size extract_size*extract_size
 * @param d_fullimg full image of size fullimg_size*fullimg_size
 * @param fullimg_size size of the d_fullimg leading dimension
 * @param center_pos position of the center of d_smallimg in d_fullimg
 * @param extract_size size of the d_smallimg leading dimension
 * @param roll get pixels as if d_fullimg need to be roll
 */
template <class T>
int32_t extract(T *d_smallimg, const T *d_fullimg, int32_t fullimg_size, int32_t center_pos,
            int32_t extract_size, bool roll) {
  extract_krnl<<<1, extract_size * extract_size>>>(
      d_smallimg, d_fullimg, fullimg_size, center_pos, extract_size, roll);

  carma_check_msg("extract_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int32_t extract(float *d_smallimg, const float *d_fullimg,
                     int32_t fullimg_size, int32_t center_pos, int32_t extract_size,
                     bool roll);
template int32_t extract(double *d_smallimg, const double *d_fullimg,
                     int32_t fullimg_size, int32_t center_pos, int32_t extract_size,
                     bool roll);
