#include <sutra_ao_utils.h>
#include "carma_utils.cuh"

__global__ void cfillrealp_krnl(cuFloatComplex *odata, float *idata, int N) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid].x = idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int cfillrealp(cuFloatComplex *d_odata, float *d_idata, int N, carma_device *device) {

  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  cfillrealp_krnl<<<grid, threads>>>(d_odata, d_idata, N);

  cutilCheckMsg("cfillrealp_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void cgetrealp_krnl(float *odata, cuFloatComplex *idata, int N) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = idata[tid].x;
    tid += blockDim.x * gridDim.x;
  }
}

int cgetrealp(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device) {

  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  cgetrealp_krnl<<<grid, threads>>>(d_odata, d_idata, N);

  cutilCheckMsg("cgetrealp_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void abs2_krnl(float *odata, cuFloatComplex *idata, int N) {
  cuFloatComplex cache;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < N) {
    cache = idata[tid];
    odata[tid] = cache.x * cache.x + cache.y * cache.y;
  }
}

int abs2(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  abs2_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  cutilCheckMsg("abs2_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void abs2c_krnl(cuFloatComplex *odata, cuFloatComplex *idata,
    int N) {
  cuFloatComplex cache;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < N) {
    cache = idata[tid];
    odata[tid].x = cache.x * cache.x + cache.y * cache.y;
    odata[tid].y = 0.0;
  }
}

int abs2c(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);
  //DEBUG_TRACE("N = %d, nthreads = %d, nblocks = %d;",N , nthreads, nblocks);
  abs2c_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  cutilCheckMsg("abs2c_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void subapnorm_krnl(float *odata, float *idata, float *fact,
    float *norm, float nphot, int n, int N) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = idata[tid] * fact[tid / n] / norm[tid / n] * nphot;
    tid += blockDim.x * gridDim.x;
  }
}

int subap_norm(float *d_odata, float *d_idata, float *fact, float *norm,
    float nphot, int n, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  subapnorm_krnl<<<grid, threads>>>(d_odata, d_idata, fact, norm, nphot, n, N);
  cutilCheckMsg("subapnorm_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void subapnormasync_krnl(float *odata, float *idata, float *fact,
    float *norm, float nphot, int n, int N, int istart) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  tid += istart;
  if (tid < N) {
    odata[tid] = idata[tid] * fact[tid / n] / norm[tid / n] * nphot;
  }
}

int subap_norm_async(float *d_odata, float *d_idata, float *fact, float *norm,
    float nphot, int n, int N, carma_streams *streams, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  int nstreams = streams->get_nbStreams();
  getNumBlocksAndThreads(device, N / nstreams, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  for (int i = 0; i < nstreams; i++) {
    subapnormasync_krnl<<<grid, threads, 0, streams->get_stream(i)>>>(d_odata,
        d_idata, fact, norm, nphot, n, N, i * nblocks * nthreads);
    cutilCheckMsg("subapnormasync_kernel<<<>>> execution failed\n");
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
    int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  krnl_fillindx<<<grid, threads>>>(d_odata, d_idata, indx, alpha, beta, N);

  cutilCheckMsg("fillindx_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int fillindx(float *d_odata, float *d_idata, int *indx, int N, carma_device *device) {
  return fillindx(d_odata, d_idata, indx, 1.0f, 0.0f, N, device);

}
int fillindx(float *d_odata, float *d_idata, int *indx, float alpha, int N,
    carma_device *device) {
  return fillindx(d_odata, d_idata, indx, alpha, 0.0f, N, device);

}
__global__ void fillarr2d_krnl(float *odata, float *idata, int tidx0, int Ncol,
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

int fillarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
    carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  fillarr2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N);

  cutilCheckMsg("fillarr2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
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
    carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  getarr2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol, NC, N);

  cutilCheckMsg("getarr2d_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template<class T>
__global__ void addai_krnl(T *odata, T* idata, int i, int sgn, int N) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    if (sgn == 1)
      odata[tid] += idata[i];
    else
      odata[tid] -= idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

template<class T>
int addai(T *d_odata, T *i_data, int i, int sgn, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  addai_krnl<T><<<grid, threads>>>(d_odata, i_data, i, sgn, N);

  cutilCheckMsg("plusai_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template
int addai<float>(float *d_odata, float *i_data, int i, int sgn, int N, carma_device *device);
template
int addai<double>(double *d_odata, double *i_data, int i, int sgn, int N, carma_device *device);

template<class T>
__global__ void roll_krnl(T *idata, int N, int M, int Ntot) {
  T tmp;

  int tidt = threadIdx.x + blockIdx.x * blockDim.x;
  int nim = tidt / Ntot;

  int tid = tidt - nim * Ntot;

  while (tid < Ntot) {

    int x = tid % N;
    int y = tid / N;

    int xx = (x + N / 2) % N;
    int yy = (y + M / 2) % M;
    int tid2 = xx + yy * N;

    tmp = idata[tid + nim * (N * M)];
    idata[tid + nim * (N * M)] = idata[tid2 + nim * (N * M)];
    idata[tid2 + nim * (N * M)] = tmp;

    tid += blockDim.x * gridDim.x;
  }
}

template<class T>
int roll(T *idata, int N, int M, int nim, carma_device *device) {

  long Ntot = N * M * nim;
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, Ntot / 2, nBlocks, nThreads);

  dim3 grid(nBlocks), threads(nThreads);

  roll_krnl<T><<<grid, threads>>>(idata, N, M, Ntot / 2);

  cutilCheckMsg("roll_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;

}

template int
roll<float>(float *idata, int N, int M, int nim, carma_device *device);

template int
roll<double>(double *idata, int N, int M, int nim, carma_device *device);

template int
roll<cuFloatComplex>(cuFloatComplex *idata, int N, int M, int nim, carma_device *device);

template<class T>
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

template<class T>
int roll(T *idata, int N, int M, carma_device *device) {

  long Ntot = N * M;
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, Ntot / 2, nBlocks, nThreads);

  dim3 grid(nBlocks), threads(nThreads);

  roll_krnl<T><<<grid, threads>>>(idata, N, M);

  cutilCheckMsg("roll_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;

}

template int
roll<float>(float *idata, int N, int M, carma_device *device);

template int
roll<double>(double *idata, int N, int M, carma_device *device);

template int
roll<cuFloatComplex>(cuFloatComplex *idata, int N, int M, carma_device *device);

template<class T>
__global__ void avg_krnl(T *data, T *p_sum){
	T *sdata = SharedMemory<T>();
	//Load shared memory
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	int sid = threadIdx.x;

	sdata[sid] = data[tid];

	__syncthreads();

	reduce_krnl(sdata,blockDim.x,sid);

	__syncthreads();

	if (threadIdx.x == 0)
		p_sum[blockIdx.x] = sdata[0];
}

template<class T>
__global__ void remove_avg_krnl(T *data, int N, T avg){
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	while(tid < N){
		data[tid] -= avg;
		tid += blockDim.x * gridDim.x;
	}
}

template<class T>
int remove_avg(T *data, int N, carma_device *device){

	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	int smemSize = nthreads * sizeof(T);

	T p_sum_c[nblocks];
	T *p_sum;
	cutilSafeCall(
	      cudaMalloc((void** )&(p_sum), sizeof(T) * nblocks));

	avg_krnl<<< grid, threads, smemSize>>>(data,p_sum);
	cutilCheckMsg("avg_krnl<<<>>> execution failed\n");
	cutilSafeCall(
			cudaMemcpy(p_sum_c,p_sum,nblocks*sizeof(T),cudaMemcpyDeviceToHost));
	cutilSafeCall(
			cudaFree(p_sum));

	T avg = 0;
	for(int i=0 ; i<nblocks ; i++){
		avg += p_sum_c[i];
	}
	avg /= N;

	remove_avg_krnl<<< grid, threads >>>(data,N,avg);
	cutilCheckMsg("remove_avg_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}

template int
remove_avg<float>(float *data, int N, carma_device *device);
template int
remove_avg<double>(double *data, int N, carma_device *device);
