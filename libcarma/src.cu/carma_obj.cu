#include <carma_obj.h>
#include <carma_utils.cuh>

/*
 short	2 bytes
 int	4 bytes
 long	4 bytes
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
__global__ void generic1d(T *odata, T *idata, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = carma_sin(2.0f * idata[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int launch_generic1d(T *d_odata, T *d_idata, int N, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  generic1d<<<grid, threads>>>(d_odata, d_idata, N);

  return EXIT_SUCCESS;
}

template int launch_generic1d<float>(float *d_odata, float *d_idata, int N,
                                     carma_device *device);

template int launch_generic1d<double>(double *d_odata, double *d_idata, int N,
                                      carma_device *device);

template <class T>
__global__ void generic2d(T *odata, T *idata, int N) {
  __shared__ T cache[BLOCK_SZ][BLOCK_SZ];

  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  int tid = x + y * blockDim.x * gridDim.x;

  cache[BLOCK_SZ - 1 - threadIdx.x][BLOCK_SZ - 1 - threadIdx.y] =
      carma_sin(2.0f * idata[tid]);

  __syncthreads();

  odata[tid] = cache[BLOCK_SZ - 1 - threadIdx.x][BLOCK_SZ - 1 - threadIdx.y];
}

template <class T>
int launch_generic2d(T *d_odata, T *d_idata, int N1, int N2) {
  dim3 blocks(N1 / BLOCK_SZ, N2 / BLOCK_SZ), threads(BLOCK_SZ, BLOCK_SZ);
  int N = N1 * N2;

  generic2d<<<blocks, threads>>>(d_odata, d_idata, N);

  return EXIT_SUCCESS;
}

template int launch_generic2d<float>(float *d_odata, float *d_idata, int N1,
                                     int N2);
template int launch_generic2d<double>(double *d_odata, double *d_idata, int N1,
                                      int N2);

template <class T_data>
__global__ void krnl_clip(T_data *data, T_data min, T_data max, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    data[tid] = carma_clip<T_data>(data[tid], min, max);
    tid += blockDim.x * gridDim.x;
  }
}

template <class T_data>
void clip_array(T_data *d_data, T_data min, T_data max, int N,
                carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  krnl_clip<<<grid, threads>>>(d_data, min, max, N);
  carmaCheckMsg("krnl_clip<<<>>> execution failed\n");
}

template void
clip_array<int>(int *d_data, int min, int max, int N, carma_device *device);
template void
clip_array<unsigned int>(unsigned int *d_data, unsigned int min, unsigned int max, int N, carma_device *device);
template void
clip_array<float>(float *d_data, float min, float max, int N, carma_device *device);
template void
clip_array<double>(double *d_data, double min, double max, int N, carma_device *device);
template<> void
clip_array(cuFloatComplex *d_data, cuFloatComplex min, cuFloatComplex max, int N, carma_device *device) {
  throw "not implemented";
}
template<> void
clip_array(cuDoubleComplex *d_data, cuDoubleComplex min, cuDoubleComplex max, int N, carma_device *device) {
  throw "not implemented";
}

template <class T>
__global__ void krnl_fillindex(T *odata, T *idata, int *indx, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = idata[indx[tid]];
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int fillindex(T *d_odata, T *d_idata, int *indx, int N, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  krnl_fillindex<<<grid, threads>>>(d_odata, d_idata, indx, N);
  carmaCheckMsg("krnl_fillindex<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int fillindex<float>(float *d_odata, float *d_idata, int *indx, int N,
                              carma_device *device);

template int fillindex<double>(double *d_odata, double *d_idata, int *indx,
                               int N, carma_device *device);

template <class T>
__global__ void krnl_fillvalues(T *odata, T val, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = val / N;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int fillvalues(T *d_odata, T val, int N, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  krnl_fillvalues<<<grid, threads>>>(d_odata, val, N);
  carmaCheckMsg("krnl_fillvalues<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int fillvalues<float>(float *d_odata, float val, int N,
                               carma_device *device);

template int fillvalues<double>(double *d_odata, double val, int N,
                                carma_device *device);

template int fillvalues<unsigned int>(unsigned int *d_odata, unsigned int val,
                                      int N, carma_device *device);

template <class T>
__global__ void getarray2d_krnl(T *odata, T *idata, int tidx0, int Ncol, int NC,
                                int N) {
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

template <class T>
int getarray2d(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
               carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  getarray2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N);

  carmaCheckMsg("getarray2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int getarray2d<float>(float *d_odata, float *d_idata, int x0, int Ncol,
                               int NC, int N, carma_device *device);

template int getarray2d<double>(double *d_odata, double *d_idata, int x0,
                                int Ncol, int NC, int N, carma_device *device);

template <class T>
__global__ void fillarray2d_krnl(T *odata, T *idata, int tidx0, int Ncol,
                                 int NC, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int tidB;

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
int fillarray2d(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
                carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  fillarray2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N);

  carmaCheckMsg("fillarray2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int fillarray2d<float>(float *d_odata, float *d_idata, int x0,
                                int Ncol, int NC, int N, carma_device *device);

template int fillarray2d<double>(double *d_odata, double *d_idata, int x0,
                                 int Ncol, int NC, int N, carma_device *device);

template <class T>
__global__ void plus_krnl(T *idata, T alpha, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    idata[tid] += alpha;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T>
int carma_plus(T *d_odata, T alpha, int N, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  plus_krnl<<<grid, threads>>>(d_odata, alpha, N);

  carmaCheckMsg("plus_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int carma_plus<float>(float *d_odata, float alpha, int N,
                               carma_device *device);

template int carma_plus<double>(double *d_odata, double alpha, int N,
                                carma_device *device);

template <class T>
__global__ void plusai_krnl(T *odata, T *idata, int i, int sgn, int N) {
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
int carma_plusai(T *d_odata, T *i_data, int i, int sgn, int N,
                 carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  plusai_krnl<<<grid, threads>>>(d_odata, i_data, i, sgn, N);

  carmaCheckMsg("plusai_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int carma_plusai<float>(float *d_odata, float *d_idata, int i, int sgn,
                                 int N, carma_device *device);

template int carma_plusai<double>(double *d_odata, double *d_idata, int i,
                                  int sgn, int N, carma_device *device);

template <class T>
__global__ void fillarray2d2_krnl(T *odata, T *idata, int tidx0, int Ncol,
                                  int NC, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int tidB;

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
int fillarray2d2(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
                 carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  fillarray2d2_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N);

  carmaCheckMsg("fillarray2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int fillarray2d2<float>(float *d_odata, float *d_idata, int x0,
                                 int Ncol, int NC, int N, carma_device *device);

template int fillarray2d2<double>(double *d_odata, double *d_idata, int x0,
                                  int Ncol, int NC, int N,
                                  carma_device *device);

template <class T>
__global__ void kern_fill_sym_matrix(char src_uplo, T *data, int Ncol, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int tidB;

  while (tid < N) {
    int row = tid / Ncol;
    int col = tid - row * Ncol;

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
int fill_sym_matrix(char src_uplo, T *d_data, int Ncol, int N,
                    carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, N, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  kern_fill_sym_matrix<<<grid, threads>>>(src_uplo, d_data, Ncol, N);

  carmaCheckMsg("kern_fill_sym_matrix<<<>>> execution failed");

  return EXIT_SUCCESS;
}

template int fill_sym_matrix<float>(char uplo, float *d_data, int Ncol, int N,
                                    carma_device *device);

template int fill_sym_matrix<double>(char uplo, double *d_data, int Ncol, int N,
                                     carma_device *device);
