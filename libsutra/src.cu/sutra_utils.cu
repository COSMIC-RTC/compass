// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_utils.cu
//! \ingroup   libsutra
//! \class     sutra_utils
//! \brief     this file provides utilities to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_utils.h>
#include "carma_utils.cuh"

int compute_nmaxhr(long nvalid) {
  // this is the big array => we use nmaxhr and treat it sequentially

  int mnmax = 500;
  int nmaxhr = mnmax;
  if (nvalid > 2 * mnmax) {
    int tmp0 = nvalid % mnmax;
    int tmp = 0;
    for (int cc = 1; cc < mnmax / 5; cc++) {
      tmp = nvalid % (mnmax + cc);
      if ((tmp > tmp0) || (tmp == 0)) {
        if (tmp == 0)
          tmp0 = 2 * mnmax;
        else
          tmp = tmp0;

        nmaxhr = mnmax + cc;
      }
    }
    return nmaxhr;
  }
  return nvalid;
}

__global__ void cfillrealp_krnl(cuFloatComplex *odata, float *idata, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid].x = idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int cfillrealp(cuFloatComplex *d_odata, float *d_idata, int N,
               CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  cfillrealp_krnl<<<grid, threads>>>(d_odata, d_idata, N);

  carma_check_msg("cfillrealp_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void cgetrealp_krnl(float *odata, cuFloatComplex *idata, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = idata[tid].x;
    tid += blockDim.x * gridDim.x;
  }
}

int cgetrealp(float *d_odata, cuFloatComplex *d_idata, int N,
              CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  cgetrealp_krnl<<<grid, threads>>>(d_odata, d_idata, N);

  carma_check_msg("cgetrealp_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void abs2_krnl(float *odata, cuFloatComplex *idata, int N) {
  cuFloatComplex cache;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    cache = idata[tid];
    odata[tid] = cache.x * cache.x + cache.y * cache.y;
    tid += blockDim.x * gridDim.x;
  }
}

int abs2(float *d_odata, cuFloatComplex *d_idata, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  abs2_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carma_check_msg("abs2_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void abs2_krnl(float *odata, cuFloatComplex *idata, int N,
                          float fact) {
  cuFloatComplex cache;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    cache = idata[tid];
    odata[tid] += (fact * (cache.x * cache.x + cache.y * cache.y));
    tid += blockDim.x * gridDim.x;
  }
}

int abs2(float *d_odata, cuFloatComplex *d_idata, int N, float fact,
         CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  abs2_krnl<<<grid, threads>>>(d_odata, d_idata, N, fact);
  carma_check_msg("abs2_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void abs2c_krnl(cuFloatComplex *odata, cuFloatComplex *idata,
                           int N) {
  cuFloatComplex cache;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    cache = idata[tid];
    odata[tid].x = cache.x * cache.x + cache.y * cache.y;
    odata[tid].y = 0.0;
    tid += blockDim.x * gridDim.x;
  }
}

int abs2c(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
          CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);
  // DEBUG_TRACE("N = %d, nb_threads = %d, nb_blocks = %d;",N , nb_threads, nb_blocks);
  abs2c_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carma_check_msg("abs2c_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void subapnorm_krnl(float *odata, float *idata, float *fact,
                               float *norm, float nphot, int n, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    if (norm[tid / n] != 0) {
      odata[tid] = idata[tid] * fact[tid / n] / norm[tid / n] * nphot;
    }
    tid += blockDim.x * gridDim.x;
  }
}

int subap_norm(float *d_odata, float *d_idata, float *fact, float *norm,
               float nphot, int n, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  subapnorm_krnl<<<grid, threads>>>(d_odata, d_idata, fact, norm, nphot, n, N);
  carma_check_msg("subapnorm_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void subapnormasync_krnl(float *odata, float *idata, float *fact,
                                    float *norm, float nphot, int n, int N,
                                    int istart) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  tid += istart;
  while (tid < N) {
    odata[tid] = idata[tid] * fact[tid / n] / norm[tid / n] * nphot;
    tid += blockDim.x * gridDim.x;
  }
}

int subap_norm_async(float *d_odata, float *d_idata, float *fact, float *norm,
                     float nphot, int n, int N, CarmaStreams *streams,
                     CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  int nstreams = streams->get_nb_streams();
  get_num_blocks_and_threads(device, N / nstreams, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  for (int i = 0; i < nstreams; i++) {
    subapnormasync_krnl<<<grid, threads, 0, streams->get_stream(i)>>>(
        d_odata, d_idata, fact, norm, nphot, n, N, i * nb_blocks * nb_threads);
    carma_check_msg("subapnormasync_kernel<<<>>> execution failed\n");
  }

  return EXIT_SUCCESS;
}

__global__ void krnl_fillindx(float *odata, float *idata, int *indx,
                              float alpha, float beta, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = (alpha * idata[indx[tid]]) + beta;
    tid += blockDim.x * gridDim.x;
  }
}

int fillindx(float *d_odata, float *d_idata, int *indx, float alpha, float beta,
             int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  krnl_fillindx<<<grid, threads>>>(d_odata, d_idata, indx, alpha, beta, N);

  carma_check_msg("fillindx_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int fillindx(float *d_odata, float *d_idata, int *indx, int N,
             CarmaDevice *device) {
  return fillindx(d_odata, d_idata, indx, 1.0f, 0.0f, N, device);
}
int fillindx(float *d_odata, float *d_idata, int *indx, float alpha, int N,
             CarmaDevice *device) {
  return fillindx(d_odata, d_idata, indx, alpha, 0.0f, N, device);
}
__global__ void fillarr2d_krnl(float *odata, float *idata, int tidx0, int Ncol,
                               int NC, int N, int dir) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int tidB;

  while (tid < N) {
    if (Ncol > 1)
      tidB = tidx0 + (tid / Ncol) * NC + (tid % Ncol);
    else
      tidB = tidx0 + tid * NC;
    if (dir > 0)
      odata[tidB] = idata[tid];
    else
      odata[tidB] = idata[N - 1 - tid];
    tid += blockDim.x * gridDim.x;
  }
}

int fillarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
              int dir, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  fillarr2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N, dir);

  carma_check_msg("fillarr2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
int fillarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
              CarmaDevice *device) {
  return fillarr2d(d_odata, d_idata, x0, Ncol, NC, N, 1, device);
}
__global__ void getarr2d_krnl(float *odata, float *idata, int tidx0, int Ncol,
                              int NC, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int tidB;

  while (tid < N) {
    if (Ncol > 1)
      tidB = tidx0 + (tid / Ncol) * NC + (tid % Ncol);
    else
      tidB = tidx0 + tid * NC;
    odata[tid] = idata[tidB];
    tid += blockDim.x * gridDim.x;
  }
}

int getarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
             CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  getarr2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N);

  carma_check_msg("getarr2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template <class T>
__global__ void addai_krnl(T *odata, T *idata, int i, int sgn, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    if (sgn == 1)
      odata[tid] += idata[i];
    else
      odata[tid] -= idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int addai(T *d_odata, T *i_data, int i, int sgn, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  addai_krnl<T><<<grid, threads>>>(d_odata, i_data, i, sgn, N);

  carma_check_msg("plusai_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int addai<float>(float *d_odata, float *i_data, int i, int sgn, int N,
                          CarmaDevice *device);
template int addai<double>(double *d_odata, double *i_data, int i, int sgn,
                           int N, CarmaDevice *device);

template <class T>
__global__ void roll_krnl(T *idata, int N, int M, int Nim) {
  T tmp;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < (N * M / 2)) {
    int x = tid % N;
    int y = tid / N;

    int xx = (x + N / 2) % N;
    int yy = (y + M / 2) % M;
    int tid2 = xx + yy * N;

    for (int ii = 0; ii < Nim; ii++) {
      tmp = idata[tid + ii * N * M];
      idata[tid + ii * N * M] = idata[tid2 + ii * N * M];
      idata[tid2 + ii * N * M] = tmp;
    }

    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int roll(T *idata, int N, int M, int nim, CarmaDevice *device) {
  long Ntot = N * M;
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, Ntot / 2, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  roll_krnl<T><<<grid, threads>>>(idata, N, M, nim);

  carma_check_msg("roll_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int roll<float>(float *idata, int N, int M, int nim,
                         CarmaDevice *device);

template int roll<double>(double *idata, int N, int M, int nim,
                          CarmaDevice *device);

template int roll<cuFloatComplex>(cuFloatComplex *idata, int N, int M, int nim,
                                  CarmaDevice *device);

template <class T>
__global__ void roll_krnl(T *idata, int N, int M) {
  T tmp;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < (N * M / 2)) {
    int x = tid % N;
    int y = tid / N;

    int xx = (x + N / 2) % N;
    int yy = (y + M / 2) % M;
    int tid2 = xx + yy * N;

    tmp = idata[tid];
    idata[tid] = idata[tid2];
    idata[tid2] = tmp;

    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int roll(T *idata, int N, int M, CarmaDevice *device) {
  long Ntot = N * M;
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, Ntot / 2, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  roll_krnl<T><<<grid, threads>>>(idata, N, M);

  carma_check_msg("roll_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int roll<float>(float *idata, int N, int M, CarmaDevice *device);

template int roll<double>(double *idata, int N, int M, CarmaDevice *device);

template int roll<cuFloatComplex>(cuFloatComplex *idata, int N, int M,
                                  CarmaDevice *device);

template <class T>
__global__ void roll_mult_krnl(T *odata, T *idata, int N, int M, T alpha) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < (N * M / 2)) {
    int x = tid % N;
    int y = tid / N;

    int xx = (x + N / 2) % N;
    int yy = (y + M / 2) % M;
    int tid2 = xx + yy * N;

    odata[tid] = alpha * idata[tid2];
    odata[tid2] = alpha * idata[tid];

    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int roll_mult(T *odata, T *idata, int N, int M, T alpha, CarmaDevice *device) {
  long Ntot = N * M;
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, Ntot / 2, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  roll_mult_krnl<T><<<grid, threads>>>(odata, idata, N, M, alpha);

  carma_check_msg("roll_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int roll_mult<float>(float *odata, float *idata, int N, int M,
                              float alpha, CarmaDevice *device);

template int roll_mult<double>(double *odata, double *idata, int N, int M,
                               double alpha, CarmaDevice *device);

template <class T>
__global__ void avg_krnl(T *data, T *p_sum, int N) {
  T *sdata = SharedMemory<T>();
  // Load shared memory
  int tid = threadIdx.x + blockDim.x * blockIdx.x;
  int sid = threadIdx.x;

  if (tid < N)
    sdata[sid] = data[tid];
  else
    sdata[sid] = 0;

  __syncthreads();

  reduce_krnl(sdata, blockDim.x, sid);

  __syncthreads();

  if (threadIdx.x == 0) p_sum[blockIdx.x] = sdata[0];
}

template <class T>
__global__ void remove_avg_krnl(T *data, int N, T avg) {
  int tid = threadIdx.x + blockDim.x * blockIdx.x;
  while (tid < N) {
    data[tid] -= avg;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int remove_avg(T *data, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);
  int smemSize = nb_threads * sizeof(T);

  T p_sum_c[nb_blocks];
  T *p_sum;
  carma_safe_call(cudaMalloc((void **)&(p_sum), sizeof(T) * nb_blocks));

  avg_krnl<<<grid, threads, smemSize>>>(data, p_sum, N);
  carma_check_msg("avg_krnl<<<>>> execution failed\n");
  carma_safe_call(
      cudaMemcpy(p_sum_c, p_sum, nb_blocks * sizeof(T), cudaMemcpyDeviceToHost));
  carma_safe_call(cudaFree(p_sum));

  T avg = 0;
  for (int i = 0; i < nb_blocks; i++) {
    avg += p_sum_c[i];
  }
  avg /= N;
  remove_avg_krnl<<<grid, threads>>>(data, N, avg);
  carma_check_msg("remove_avg_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int remove_avg<float>(float *data, int N, CarmaDevice *device);
template int remove_avg<double>(double *data, int N, CarmaDevice *device);

__global__ void conv_krnl(cuFloatComplex *odata, cuFloatComplex *idata, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  cuFloatComplex tmp;

  while (tid < N) {
    tmp.x = idata[tid].x * odata[tid].x - idata[tid].y * odata[tid].y;
    tmp.y = idata[tid].y * odata[tid].x + idata[tid].x * odata[tid].y;
    odata[tid] = tmp;
    tid += blockDim.x * gridDim.x;
  }
}

int convolve(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
             CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  conv_krnl<<<grid, threads>>>(d_odata, d_idata, N);

  carma_check_msg("conv_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void convmod_krnl(cuFloatComplex *odata, cuFloatComplex *idata,
                             int mod, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  cuFloatComplex tmp;

  while (tid < N) {
    tmp.x = (idata[tid].x * odata[tid].x - idata[tid].y * odata[tid].y) / mod;
    tmp.y = (idata[tid].y * odata[tid].x + idata[tid].x * odata[tid].y) / mod;
    odata[tid] = tmp;
    tid += blockDim.x * gridDim.x;
  }
}

int convolve_modulate(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int mod,
                      int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  convmod_krnl<<<grid, threads>>>(d_odata, d_idata, mod, N);

  carma_check_msg("conv_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template <class T>
__global__ void mult_krnl(T *i_data, T *scale, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    i_data[tid] = i_data[tid] * scale[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
__global__ void mult_krnl(T *i_data, T *scale, T gain, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    i_data[tid] = i_data[tid] * scale[tid] * gain;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
__global__ void mult_krnl(T *i_data, T gain, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    i_data[tid] = i_data[tid] * gain;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void mult_int_krnl(float *o_data, float *i_data, float gain, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    o_data[tid] = gain * i_data[tid] + o_data[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void mult_int_krnl(float *o_data, float *i_data, float *scale,
                              float gain, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    o_data[tid] = gain * (i_data[tid] * scale[tid]) + o_data[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void mult_int_krnl(float *o_data, float *i_data, float *scale,
                              float gain, int N, int istart) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  tid += istart;

  while (tid < N) {
    o_data[tid] = gain * (i_data[tid] * scale[tid]) + o_data[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void add_md_krnl(float *o_matrix, float *i_matrix, float *i_vector,
                            int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    o_matrix[tid * (N + 1)] = i_matrix[tid * (N + 1)] + i_vector[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int mult_vect(T *d_data, T *scale, int N, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  mult_krnl<<<grid, threads>>>(d_data, scale, N);

  carma_check_msg("mult_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}
template int mult_vect<float>(float *d_data, float *scale, int N,
                              CarmaDevice *device);
template int mult_vect<double>(double *d_data, double *scale, int N,
                               CarmaDevice *device);
#ifdef CAN_DO_HALF
template int mult_vect<half>(half *d_data, half *scale, int N,
                             CarmaDevice *device);
#endif

template <class T>
int mult_vect(T *d_data, T *scale, T gain, int N, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  mult_krnl<<<grid, threads>>>(d_data, scale, gain, N);

  carma_check_msg("mult_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int mult_vect<float>(float *d_data, float *scale, float gain, int N,
                              CarmaDevice *device);
template int mult_vect<double>(double *d_data, double *scale, double gain,
                               int N, CarmaDevice *device);
#ifdef CAN_DO_HALF
template int mult_vect<half>(half *d_data, half *scale, half gain, int N,
                             CarmaDevice *device);
#endif

template <class T>
int mult_vect(T *d_data, T gain, int N, CarmaDevice *device) {
  int nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  mult_krnl<<<grid, threads>>>(d_data, gain, N);

  carma_check_msg("mult_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int mult_vect<float>(float *d_data, float gain, int N,
                              CarmaDevice *device);
template int mult_vect<double>(double *d_data, double gain, int N,
                               CarmaDevice *device);
#ifdef CAN_DO_HALF
template int mult_vect<half>(half *d_data, half gain, int N,
                             CarmaDevice *device);
#endif

int mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
             CarmaDevice *device, CarmaStreams *streams) {
  int nb_threads = 0, nb_blocks = 0;

  int nstreams = streams->get_nb_streams();
  get_num_blocks_and_threads(device, N / nstreams, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  for (int i = 0; i < nstreams; i++) {
    mult_int_krnl<<<grid, threads, 0, streams->get_stream(i)>>>(
        o_data, i_data, scale, gain, N, i * nb_blocks * nb_threads);
    carma_check_msg("multint_kernel<<<>>> execution failed\n");
  }

  return EXIT_SUCCESS;
}

int mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
             CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  mult_int_krnl<<<grid, threads>>>(o_data, i_data, scale, gain, N);
  carma_check_msg("multint_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int mult_int(float *o_data, float *i_data, float gain, int N,
             CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;

  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);

  mult_int_krnl<<<grid, threads>>>(o_data, i_data, gain, N);
  carma_check_msg("multint_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int add_md(float *o_matrix, float *i_matrix, float *i_vector, int N,
           CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  add_md_krnl<<<grid, threads>>>(o_matrix, i_matrix, i_vector, N);
  carma_check_msg("add_md_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
