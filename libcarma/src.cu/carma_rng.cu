#include <carma_obj.h>

// PRNG init kernel
__global__ void initPRNG(curandState *s, int n, int *seed, int offset) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n)
    curand_init(seed[id], threadIdx.x, offset, &s[id]);
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
__global__ void carma_curand_poisson(curandState *state, float *res, int n) {
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
__global__ void carma_curand_poisson(curandState *state, double *res, int n) {
  double xm;

  double tmp, sq, alxm, g, oldm = (-1.0);
  double em, t, y;
  const int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  const int delta = blockDim.x * gridDim.x;
  for (int idx = tidx; idx < n; idx += delta) {
    //xm = (float)results[idx];
    if (xm > 0.0) {
      if (xm != oldm) {
        oldm = xm;
        sq = sqrt(2.0 * xm);
        alxm = log(xm);
        g = xm * alxm - lgamma(xm + 1.0);
      }
      do {
        do {
          tmp = curand_uniform(&state[tidx]);
          y = tan(3.1415926535897932384626433832f * tmp);
          em = sq * y + xm;
        } while (em < 0.0f);
        em = floor(em);
        t = 0.9 * (1.0 + y * y) * exp(em * alxm - lgamma(em + 1.0) - g);
        tmp = curand_uniform(&state[tidx]);
      } while (tmp > t);
    } else
      em = 0.0;
    res[idx] = xm;
  }
}

int carma_prng_init(int *seed, const int nThreads, const int nBlocks,
    curandState *state) {

  dim3 grid(nBlocks);
  dim3 threads(nThreads);

  // Initialise RNG  
  initPRNG<<<grid, threads>>>(state, nThreads * nBlocks, seed, nThreads);
  cutilCheckMsg("initRNG<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

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
  if (gtype == 'P')
    carma_curand_poisson<float> <<<nBlocks, nThreads>>>(state, results, n);

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

  return EXIT_SUCCESS;
}
