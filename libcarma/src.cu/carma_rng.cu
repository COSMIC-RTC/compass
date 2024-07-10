// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_rng.cu
//! \ingroup   libcarma
//! \brief     this file provides RNG CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24


#include <assert.h>
#include <carma_obj.hpp>
#include <curand_normal.h>
#include <curand_uniform.h>

#define CARMA_NYI_DEV                 \
  {                                   \
    printf("Method not implemented"); \
    assert(0);                        \
  }

// PRNG init kernel
__global__ void initPRNG(curandState *s, int32_t n, int32_t *seed, int32_t offset) {
  const int32_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) curand_init(seed[id], threadIdx.x, offset, &s[id]);
}

int32_t carma_prng_init(int32_t *seed, const int32_t nb_threads, const int32_t nb_blocks,
                    curandState *state) {
  dim3 grid(nb_blocks);
  dim3 threads(nb_threads);

  // Initialise RNG
  initPRNG<<<grid, threads>>>(state, nb_threads * nb_blocks, seed, nb_threads);
  carma_check_msg("initRNG<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template <class T>
__forceinline__ __device__ void carma_curand_uniform_gen(
    curandState *state, T (*fct)(curandState *), T *res, int32_t n, float beta) {
  const int32_t tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int32_t delta = blockDim.x * gridDim.x;
  for (int32_t idx = tidx; idx < n; idx += delta)
    res[idx] = beta * res[idx] + fct(&state[tidx]);
}

__forceinline__ __device__ void carma_curand_uniform_dev(curandState *s,
                                                         float *d, int32_t n,
                                                         float beta) {
  carma_curand_uniform_gen(s, curand_uniform, d, n, beta);
}

__forceinline__ __device__ void carma_curand_uniform_dev(curandState *s,
                                                         double *d, int32_t n,
                                                         float beta) {
  carma_curand_uniform_gen(s, curand_uniform_double, d, n, beta);
}

template <class T>
__global__ void carma_curand_uniform(curandState *s, T *d, int32_t n, float beta) {
  carma_curand_uniform_dev(s, d, n, beta);
}

template <>
__global__ void carma_curand_uniform(curandState *s, int32_t *d, int32_t n,
                                     float beta) CARMA_NYI_DEV;

template <>
__global__ void carma_curand_uniform(curandState *s, uint32_t *d, int32_t n,
                                     float beta) CARMA_NYI_DEV;

template <>
__global__ void carma_curand_uniform(curandState *s, cuFloatComplex *d, int32_t n,
                                     float beta) {
  carma_curand_uniform_gen(s, curand_uniform, (float *)d, n * 2, beta);
}

template <>
__global__ void carma_curand_uniform(curandState *s, cuDoubleComplex *d, int32_t n,
                                     float beta) {
  carma_curand_uniform_gen(s, curand_uniform_double, (double *)d, n * 2, beta);
}

template <class T>
__forceinline__ __device__ void carma_curand_normal_gen(curandState *state,
                                                        T (*fct)(curandState *),
                                                        T *res, int32_t n,
                                                        float alpha,
                                                        float beta) {
  const int32_t tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int32_t delta = blockDim.x * gridDim.x;
  for (int32_t idx = tidx; idx < n; idx += delta)
    res[idx] = beta * res[idx] + alpha * fct(&state[tidx]);
}

template <class T>
__forceinline__ __device__ void carma_curand_normal_dev(
    curandState *s, T *d, int32_t n, float alpha, float beta) CARMA_NYI_DEV

    template <>
    __forceinline__ __device__
    void carma_curand_normal_dev(curandState *s, float *d, int32_t n, float alpha,
                                 float beta) {
  carma_curand_normal_gen(s, curand_normal, d, n, alpha, beta);
}

template <>
__forceinline__ __device__ void carma_curand_normal_dev(curandState *s,
                                                        double *d, int32_t n,
                                                        float alpha,
                                                        float beta) {
  carma_curand_normal_gen(s, curand_normal_double, d, n, alpha, beta);
}

template <class T>
__global__ void carma_curand_normal(curandState *s, T *d, int32_t n, float alpha,
                                    float beta) {
  carma_curand_normal_dev(s, d, n, alpha, beta);
}

template <class T>
__global__ void carma_curand_poisson(curandState *state, T *res, int32_t n) {
  const int32_t tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int32_t delta = blockDim.x * gridDim.x;
  for (int32_t idx = tidx; idx < n; idx += delta)
    res[idx] = (T)curand_poisson(&state[tidx], (double)res[idx]);
}

template __global__ void carma_curand_poisson(curandState *s, int32_t *d, int32_t n);
template __global__ void carma_curand_poisson(curandState *s, uint32_t *d,
                                              int32_t n);
template __global__ void carma_curand_poisson(curandState *s, float *d, int32_t n);
template __global__ void carma_curand_poisson(curandState *s, double *d, int32_t n);
template <>
__global__ void carma_curand_poisson(curandState *s, cuFloatComplex *d,
                                     int32_t n) CARMA_NYI_DEV template <>
__global__ void carma_curand_poisson(curandState *s, cuDoubleComplex *d,
                                     int32_t n) CARMA_NYI_DEV

    /*
    template< class T_data, T_data (*ptr_sqrt)(T_data val),
        T_data (*ptr_log)(T_data val), T_data (*ptr_lgamma)(T_data val),
        T_data (*ptr_tan)(T_data val), T_data (*ptr_floor)(T_data val),
        T_data (*ptr_exp)(T_data val)>
    __global__ void carma_curand_montagn(curandState *state, T_data *res, int32_t n)
    { T_data xm; T_data tmp, sq, alxm, g, oldm = (-1.0); T_data em, t, y;

      const int32_t tidx = blockIdx.x * blockDim.x + threadIdx.x;
      const int32_t delta = blockDim.x * gridDim.x;
      for (int32_t idx = tidx; idx < n; idx += delta) {
        xm = res[idx];
        //xm = (T_data)results[idx];
        if (xm > 0.0f) {
          if (xm != oldm) {
            oldm = xm;
            sq = ptr_sqrt(2.0f * xm);
            alxm = ptr_log(xm);
            g = xm * alxm - ptr_lgamma(xm + 1.0f);
          }
          do {
            do {
              tmp = curand_uniform(&state[tidx]);
              y = ptr_tan(CARMA_PI * tmp);
              em = sq * y + xm;
            } while (em < 0.0f);
            em = ptr_floor(em);
            t = 0.9f * (1.0 + y * y) * ptr_exp(em * alxm - ptr_lgamma(em + 1.0f)
    - g); tmp = curand_uniform(&state[tidx]); } while (tmp > t); } else em =
    0.0f; res[idx] = xm;
      }
    }

    template<float, sqrtf, logf, lgammaf, tanf, floorf, expf>
    __global__ void carma_curand_montagn(curandState *state, float *res, int32_t n);
    template<double, sqrt, log, lgamma, tan, floor, exp>
    __global__ void carma_curand_montagn(curandState *state, double *res, int32_t
    n);
    */

    template <class T>
    __global__ void carma_curand_montagn_krn(curandState *state, T *res,
                                             int32_t n) CARMA_NYI_DEV

    template <>
    __global__
    void carma_curand_montagn_krn(curandState *state, float *res, int32_t n) {
  float xm;
  float tmp, sq, alxm, g, oldm = (-1.0);
  float em, t, y;

  const int32_t tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int32_t delta = blockDim.x * gridDim.x;
  for (int32_t idx = tidx; idx < n; idx += delta) {
    xm = res[idx];
    // xm = (float)results[idx];
    if (xm > 0.0f) {
      if (xm != oldm) {
        oldm = xm;
        sq = sqrtf(2.0f * xm);
        alxm = logf(xm);
        g = xm * alxm - lgammaf(xm + 1.0f);
      }
      do {
        do {
          tmp = curand_uniform(&state[tidx]);
          y = tanf(CARMA_PI * tmp);
          em = sq * y + xm;
        } while (em < 0.0f);
        em = floorf(em);
        t = 0.9f * (1.0 + y * y) * expf(em * alxm - lgammaf(em + 1.0f) - g);
        tmp = curand_uniform(&state[tidx]);
      } while (tmp > t);
    } else
      em = 0.0f;
    res[idx] = xm;
  }
}

template <>
__global__ void carma_curand_montagn_krn(curandState *state, double *res,
                                         int32_t n) {
  double tmp;
  const int32_t tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int32_t delta = blockDim.x * gridDim.x;
  for (int32_t idx = tidx; idx < n; idx += delta) {
    tmp = curand_uniform(&state[tidx]);
    res[idx] = tmp;
  }
}

template <class T>
int32_t carma_curand_montagn(curandState *state, T *d_odata, int32_t N,
                         CarmaDevice *device) {
  int32_t nb_blocks, nb_threads;
  get_num_blocks_and_threads(device, N, nb_blocks, nb_threads);

  dim3 grid(nb_blocks), threads(nb_threads);
  //  dim3 grid(128), threads(128);

  carma_curand_montagn_krn<<<grid, threads>>>(state, d_odata, N);

  return EXIT_SUCCESS;
}

template int32_t carma_curand_montagn<float>(curandState *state, float *d_odata,
                                         int32_t N, CarmaDevice *device);

template int32_t carma_curand_montagn<double>(curandState *state, double *d_odata,
                                          int32_t N, CarmaDevice *device);

template int32_t carma_curand_montagn<cuFloatComplex>(curandState *state,
                                                  cuFloatComplex *d_odata,
                                                  int32_t N, CarmaDevice *device);

template int32_t carma_curand_montagn<cuDoubleComplex>(curandState *state,
                                                   cuDoubleComplex *d_odata,
                                                   int32_t N, CarmaDevice *device);

template int32_t carma_curand_montagn<int32_t>(curandState *state, int32_t *d_odata, int32_t N,
                                       CarmaDevice *device);

template int32_t carma_curand_montagn<uint32_t>(curandState *state,
                                                uint32_t *d_odata, int32_t N,
                                                CarmaDevice *device);

template <>
int32_t carma_curand_montagn<uint16_t>(curandState *state, uint16_t *d_odata, int32_t N,
                                   CarmaDevice *device) {
  CARMA_NYI_DEV;
  return EXIT_FAILURE;
}

template <class T>
int32_t carma_prng_cu(T *results, const int32_t nb_threads, const int32_t nb_blocks,
                  curandState *state, char gtype, int32_t n, float alpha,
                  float beta) {
  // dim3 grid(1);
  dim3 threads(2 * nb_threads);

  if (gtype == 'U')
    carma_curand_uniform<<<nb_blocks, nb_threads>>>(state, results, n, beta);
  if (gtype == 'N')
    carma_curand_normal<<<nb_blocks, nb_threads>>>(state, results, n, alpha, beta);
  if (gtype == 'P') {
    carma_curand_poisson<<<nb_blocks, nb_threads>>>(state, results, n);
  }
  carma_check_msg("PRNG<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template int32_t carma_prng_cu(int32_t *results, const int32_t nb_threads, const int32_t nb_blocks,
                           curandState *state, char gtype, int32_t n, float alpha,
                           float beta);
template int32_t carma_prng_cu(uint32_t *results, const int32_t nb_threads,
                           const int32_t nb_blocks, curandState *state, char gtype,
                           int32_t n, float alpha, float beta);

template int32_t carma_prng_cu(float *results, const int32_t nb_threads,
                           const int32_t nb_blocks, curandState *state, char gtype,
                           int32_t n, float alpha, float beta);
template int32_t carma_prng_cu(double *results, const int32_t nb_threads,
                           const int32_t nb_blocks, curandState *state, char gtype,
                           int32_t n, float alpha, float beta);
template int32_t carma_prng_cu(cuFloatComplex *results, const int32_t nb_threads,
                           const int32_t nb_blocks, curandState *state, char gtype,
                           int32_t n, float alpha, float beta);
template int32_t carma_prng_cu(cuDoubleComplex *results, const int32_t nb_threads,
                           const int32_t nb_blocks, curandState *state, char gtype,
                           int32_t n, float alpha, float beta);

template <>
int32_t carma_prng_cu(uint16_t *results, const int32_t nb_threads, const int32_t nb_blocks,
                  curandState *state, char gtype, int32_t n, float alpha,
                  float beta) {
  CARMA_NYI_DEV;
  return EXIT_FAILURE;
}
