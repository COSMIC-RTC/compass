#include <carma_obj.h>

// PRNG init kernel
__global__ void initPRNG(curandState *s, int n, int *seed, int offset) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n)
    curand_init(seed[id], threadIdx.x, offset, &s[id]);
}

int carma_prng_init(int *seed, const int nThreads, const int nBlocks,
                    curandState *state) {

  dim3 grid(nBlocks);
  dim3 threads(nThreads);

  // Initialise RNG
  initPRNG<<<grid, threads>>>(state, nThreads * nBlocks, seed, nThreads);
  carmaCheckMsg("initRNG<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template<class T>
__global__ void
carma_curand_uniform(curandState *state, T *res, int n, float beta);

template<>
__global__ void carma_curand_uniform(curandState *s, float *d, int n,
                                     float beta) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta)
    d[idx] = beta * d[idx] + curand_uniform(&s[tidx]);
}

template<>
__global__ void carma_curand_uniform(curandState *s, double *d, int n,
                                     float beta) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta)
    d[idx] = beta * d[idx] + curand_uniform_double(&s[tidx]);
}

template<>
__global__ void carma_curand_uniform(curandState *s, cuFloatComplex *d, int n,
                                     float beta) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta) {
    d[idx].x = beta * d[idx].x + curand_uniform(&s[tidx]);
    d[idx].y = beta * d[idx].y + curand_uniform(&s[tidx]);
  }
}

template<>
__global__ void carma_curand_uniform(curandState *s, cuDoubleComplex *d, int n,
                                     float beta) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta) {
    d[idx].x = beta * d[idx].x + curand_uniform_double(&s[tidx]);
    d[idx].y = beta * d[idx].y + curand_uniform_double(&s[tidx]);
  }
}

template<class T>
__global__ void
carma_curand_normal(curandState *state, T *res, int n, float alpha, float beta);

template<>
__global__ void carma_curand_normal(curandState *s, float *d, int n,
                                    float alpha, float beta) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta)
    d[idx] = beta * d[idx] + alpha * curand_normal(&s[tidx]);
}

template<>
__global__ void carma_curand_normal(curandState *s, double *d, int n,
                                    float alpha, float beta) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta)
    d[idx] = beta * d[idx] + alpha * curand_normal_double(&s[tidx]);
}

template<>
__global__ void carma_curand_normal(curandState *s, cuFloatComplex *d, int n,
                                    float alpha, float beta) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta) {
    d[idx].x = beta * d[idx].x + alpha * curand_normal(&s[tidx]);
    d[idx].y = beta * d[idx].y + alpha * curand_normal(&s[tidx]);
  }
}

template<>
__global__ void carma_curand_normal(curandState *s, cuDoubleComplex *d, int n,
                                    float alpha, float beta) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta) {
    d[idx].x = beta * d[idx].x + alpha * curand_normal_double(&s[tidx]);
    d[idx].y = beta * d[idx].y + alpha * curand_normal_double(&s[tidx]);
  }
}

template<class T>
__global__ void
carma_curand_poisson(curandState *state, T *res, int n);


template<>
__global__ void carma_curand_poisson(curandState *s, float *d, int n) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta)
    d[idx] = (float)curand_poisson(&s[tidx],(double)d[idx]);
}

template<>
__global__ void carma_curand_poisson(curandState *s, double *d, int n) {
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta)
    d[idx] = (double)curand_poisson(&s[tidx],d[idx]);
}

/*

template< class T_data, T_data (*ptr_sqrt)(T_data val),
    T_data (*ptr_log)(T_data val), T_data (*ptr_lgamma)(T_data val),
    T_data (*ptr_tan)(T_data val), T_data (*ptr_floor)(T_data val),
    T_data (*ptr_exp)(T_data val)>
__global__ void carma_curand_montagn(curandState *state, T_data *res, int n) {
  T_data xm;
  T_data tmp, sq, alxm, g, oldm = (-1.0);
  T_data em, t, y;

  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta) {
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
          y = ptr_tan(3.1415926535897932384626433832f * tmp);
          em = sq * y + xm;
        } while (em < 0.0f);
        em = ptr_floor(em);
        t = 0.9f * (1.0 + y * y) * ptr_exp(em * alxm - ptr_lgamma(em + 1.0f) - g);
        tmp = curand_uniform(&state[tidx]);
      } while (tmp > t);
    } else
      em = 0.0f;
    res[idx] = xm;
  }
}

template<float, sqrtf, logf, lgammaf, tanf, floorf, expf>
__global__ void carma_curand_montagn(curandState *state, float *res, int n);
template<double, sqrt, log, lgamma, tan, floor, exp>
__global__ void carma_curand_montagn(curandState *state, double *res, int n);

*/
template<class T>
__global__ void
carma_curand_montagn_krn(curandState *state, T *res, int n);

template<>
__global__ void carma_curand_montagn_krn(curandState *state, float *res, int n) {
  float xm;
  float tmp, sq, alxm, g, oldm = (-1.0);
  float em, t, y;

  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta) {
    xm = res[idx];
    //xm = (float)results[idx];
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
          y = tanf(3.1415926535897932384626433832f * tmp);
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

template<>
__global__ void carma_curand_montagn_krn(curandState *state, double *res, int n) {
  double tmp;
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta) {
    tmp = curand_uniform(&state[tidx]);
    res[idx] = tmp;
  }
}

template<class T>
int carma_curand_montagn(curandState *state, T *d_odata, int N, carma_device *device) {

  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  carma_curand_montagn_krn<<<grid, threads>>>(state, d_odata, N);

  return EXIT_SUCCESS;
}

template int
carma_curand_montagn<float>(curandState *state, float *d_odata, int N, carma_device *device);

template int
carma_curand_montagn<double>(curandState *state, double *d_odata, int N, carma_device *device);


template<class T>
int
carma_prng_cu(T *results, const int nThreads, curandState *state, char gtype,
              int n, float alpha, float beta);

template<>
int carma_prng_cu(float *results, const int nThreads, const int nBlocks,
                  curandState *state, char gtype, int n, float alpha, float beta) {

  dim3 grid(1);
  dim3 threads(2 * nThreads);

  if (gtype == 'U')
    carma_curand_uniform<float> <<<nBlocks, nThreads>>>(state, results, n,
        beta);
  if (gtype == 'N')
    carma_curand_normal<float> <<<nBlocks, nThreads>>>(state, results, n,
        alpha, beta);
  if (gtype == 'P') {
    carma_curand_poisson<float> <<<nBlocks, nThreads>>>(state, results, n);
  }
  carmaCheckMsg("PRNG<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}
template<>
int carma_prng_cu(double *results, const int nThreads, const int nBlocks,
                  curandState *state, char gtype, int n, float alpha, float beta) {

  dim3 grid(1);
  dim3 threads(2 * nThreads);

  if (gtype == 'U')
    carma_curand_uniform<double> <<<nBlocks, nThreads>>>(state, results, n,
        beta);
  if (gtype == 'N')
    carma_curand_normal<double> <<<nBlocks, nThreads>>>(state, results, n,
        alpha, beta);
  if (gtype == 'P')
    carma_curand_poisson<double> <<<nBlocks, nThreads>>>(state, results, n);

  carmaCheckMsg("PRNG<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template<>
int carma_prng_cu(cuFloatComplex *results, const int nThreads,
                  const int nBlocks, curandState *state, char gtype, int n, float alpha,
                  float beta) {

  dim3 grid(1);
  dim3 threads(2 * nThreads);

  if (gtype == 'U')
    carma_curand_uniform<cuFloatComplex> <<<nBlocks, nThreads>>>(state,
        results, n, beta);
  if (gtype == 'N')
    carma_curand_normal<cuFloatComplex> <<<nBlocks, nThreads>>>(state,
        results, n, alpha, beta);
  //   if (gtype == 'P')
  //  carma_curand_poisson<float><<<1, nThreads>>>(state, (float *)results, 2*n);

  carmaCheckMsg("PRNG<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}
template<>
int carma_prng_cu(cuDoubleComplex *results, const int nThreads,
                  const int nBlocks, curandState *state, char gtype, int n, float alpha,
                  float beta) {

  dim3 grid(1);
  dim3 threads(2 * nThreads);

  if (gtype == 'U')
    carma_curand_uniform<cuDoubleComplex> <<<nBlocks, nThreads>>>(state,
        results, n, beta);
  if (gtype == 'N')
    carma_curand_normal<cuDoubleComplex> <<<nBlocks, nThreads>>>(state,
        results, n, alpha, beta);
  //   if (gtype == 'P')
  //  carma_curand_poisson<double><<<1, nThreads>>>(state, (double *)results, 2*n);

  carmaCheckMsg("PRNG<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}
