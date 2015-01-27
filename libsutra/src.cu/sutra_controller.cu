#include <sutra_controller.h>
#include <math_constants.h>
/*
 _  __                    _     
 | |/ /___ _ __ _ __   ___| |___ 
 | ' // _ \ '__| '_ \ / _ \ / __|
 | . \  __/ |  | | | |  __/ \__ \
|_|\_\___|_|  |_| |_|\___|_|___/
 
 */

__global__ void shift_krnl(float *data, int offset, int N) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    data[tid] = data[tid + offset * N];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void mult_krnl(float *i_data, float *scale, int N) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    i_data[tid] = i_data[tid] * scale[tid];
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void mult_krnl(float *i_data, float gain, int N) {

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

  if (tid < N) {
    o_data[tid] = gain * (i_data[tid] * scale[tid]) + o_data[tid];
  }
}

__global__ void add_md_krnl(float *o_matrix, float *i_matrix, float *i_vector,
    int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < N) {
    o_matrix[tid * (N + 1)] = i_matrix[tid * (N + 1)] + i_vector[tid];
  }
}

__global__ void
TT_filt_krnl(float *mat, int n, int N){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < N){
		tid % (n+1) ? mat[tid] *= -1.0f : mat[tid] = (float)(1.-mat[tid]);
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void
fill_filtmat_krnl(float *filter, int nactu, int N){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < N){
		filter[tid] =tid % (nactu+1) ? (float)-1./nactu : (float)(1.-1./nactu);
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void
fill_cmat_krnl(float *cmat, float *wtt, float *Mtt, long nact, long nslope, long N){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int i,j;
	while (tid < N) {
		i = tid/nact;
		j = tid - i*nact;
		//if(j < nact-2) cmat[tid] = 1.0f;
		cmat[tid] = (j < nact-2) ? wtt[j + i*(nact-2)] : Mtt[j - (nact - 2) + i*2];
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void
do_statcov_krnl(float *statcov, float *xpos, float *ypos, float norm, long dim, long N){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int i,j;
	while (tid < N) {
		i = tid/dim;
		j = tid - i*dim;
		statcov[i*dim + j] = (float)(6.88 * pow(sqrt(pow((double)(xpos[i]-xpos[j]),2) + pow((double)(ypos[i]-ypos[j]),2)),5./3.) * norm);
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void
floattodouble_krnl(float *i_data, double *o_data, int N){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < N){
		o_data[tid] = (double)i_data[tid];
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void
doubletofloat_krnl(double *i_data, float *o_data, int N){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < N){
		o_data[tid] = (float)i_data[tid];
		tid += blockDim.x * gridDim.x;
	}
}
__global__ void
pupphase_krnl(float *o_data, float *i_data, int *indx_pup, int N){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < N){
		o_data[tid] = i_data[indx_pup[tid]];
		tid += blockDim.x * gridDim.x;
	}
}

__device__ cuComplex
exp_complex(cuComplex z){
	cuComplex res;
	float t = expf(z.x);
	sincosf(z.y,&res.y,&res.x);
	res.x *= t;
	res.y *= t;

	return res;
}

__global__ void
compute_Hcor_krnl(float *o_data, int nrow, int ncol, float Fs, float Te, float gmin, float gmax, int delay){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int j = tid/nrow;
	int i = tid - j*nrow;

	float step = (Fs/2.0f - (Fs/(2*ncol)))/(ncol-1);
	float f = Fs/(2*ncol) + j*step;
	float G = gmin + i*(gmax-gmin)/(nrow-1);
	cuFloatComplex pTe = make_cuFloatComplex(0.0f,2*CUDART_PI_F*f*Te);
	cuFloatComplex moins_pTe = make_cuFloatComplex(0.0f,-2*CUDART_PI_F*f*Te);
	cuFloatComplex pTe2 = cuCmulf(pTe,pTe);
	cuFloatComplex UnMoinsepTe = make_cuFloatComplex(1.0f-exp_complex(moins_pTe).x,-exp_complex(moins_pTe).y);
	cuFloatComplex pdelay = make_cuFloatComplex(0.0f,-2*CUDART_PI_F*f*Te*delay);

	cuFloatComplex res = cuCdivf(UnMoinsepTe,pTe2);
	cuFloatComplex Hbo = cuCmulf(res,exp_complex(pdelay));
	Hbo.x = 1 + G*Hbo.x;
	Hbo.y *= G;
	float mod = cuCabsf(Hbo);

	o_data[tid] = 1/(mod * mod);

}

__global__ void
absnormfft_krnl(cuFloatComplex *idata, float *odata, int N, float norm){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	cuFloatComplex cache;
	while (tid < N){
		cache = idata[tid + 1]; // Reject first element (0 frequency)
		odata[tid] = (cache.x * cache.x + cache.y * cache.y) * norm;
		tid += blockDim.x * gridDim.x;
	}

}

/*
  _                           _
 | |    __ _ _   _ _ __   ___| |__   ___ _ __ ___ 
 | |   / _` | | | | '_ \ / __| '_ \ / _ \ '__/ __|
 | |__| (_| | |_| | | | | (__| | | |  __/ |  \__ \
 |_____\__,_|\__,_|_| |_|\___|_| |_|\___|_|  |___/
 
 */

int shift_buf(float *d_data, int offset, int N, carma_device *device) {

  int maxThreads = device->get_properties().maxThreadsPerBlock;
  int nBlocks = device->get_properties().multiProcessorCount * 8;
  int nThreads = (N + nBlocks - 1) / nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads - 1) / nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);

  shift_krnl<<<grid, threads>>>(d_data, offset, N);

  cutilCheckMsg("shift_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int mult_vect(float *d_data, float *scale, int N, carma_device *device) {
  int maxThreads = device->get_properties().maxThreadsPerBlock;
  int nBlocks = device->get_properties().multiProcessorCount * 8;
  int nThreads = (N + nBlocks - 1) / nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads - 1) / nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);

  mult_krnl<<<grid, threads>>>(d_data, scale, N);

  cutilCheckMsg("mult_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int mult_vect(float *d_data, float gain, int N, carma_device *device) {

  int maxThreads = device->get_properties().maxThreadsPerBlock;
  int nBlocks = device->get_properties().multiProcessorCount * 8;
  int nThreads = (N + nBlocks - 1) / nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads - 1) / nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);

  mult_krnl<<<grid, threads>>>(d_data, gain, N);

  cutilCheckMsg("mult_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
    carma_device *device, carma_streams *streams) {

  int nthreads = 0, nblocks = 0;

  int nstreams = streams->get_nbStreams();
  getNumBlocksAndThreads(device, N / nstreams, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  for (int i = 0; i < nstreams; i++) {
    mult_int_krnl<<<grid, threads, 0, streams->get_stream(i)>>>(o_data, i_data,
        scale, gain, N, i * nblocks * nthreads);
    cutilCheckMsg("multint_kernel<<<>>> execution failed\n");
  }

  return EXIT_SUCCESS;
}

int mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
    carma_device *device) {

  int nthreads = 0, nblocks = 0;

  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  mult_int_krnl<<<grid, threads>>>(o_data, i_data, scale, gain, N);
  cutilCheckMsg("multint_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int mult_int(float *o_data, float *i_data, float gain, int N,carma_device *device) {

  int nthreads = 0, nblocks = 0;

  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  mult_int_krnl<<<grid, threads>>>(o_data, i_data, gain, N);
  cutilCheckMsg("multint_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int add_md(float *o_matrix, float *i_matrix, float *i_vector, int N,
    carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  add_md_krnl<<<grid, threads>>>(o_matrix, i_matrix, i_vector, N);
  cutilCheckMsg("add_md_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int
fill_filtmat(float *filter, int nactu, int N, carma_device *device){
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);

	fill_filtmat_krnl<<<grid, threads>>>(filter, nactu, N);
	cutilCheckMsg("fill_filtmat_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}
int
TT_filt(float *mat, int n, carma_device *device){
	int nthreads = 0, nblocks = 0;
	int N = n*n;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);

	TT_filt_krnl<<<grid, threads>>>(mat,n,N);
	cutilCheckMsg("TT_filt_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}

int
fill_cmat(float *cmat, float *wtt, float *Mtt, int nactu, int nslopes, carma_device *device){
	int nthreads = 0, nblocks = 0;
	int N = nactu * nslopes;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);

	fill_cmat_krnl<<<grid, threads>>>(cmat,wtt,Mtt,nactu,nslopes,N);
	cutilCheckMsg("fill_cmat_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}

int
do_statmat(float *statcov, long dim, float *xpos, float *ypos, float norm, carma_device *device){
	int nthreads = 0, nblocks = 0;
	int N = (dim * dim);
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	do_statcov_krnl<<<grid , threads>>>(statcov,xpos,ypos,norm,dim,N);
	cutilCheckMsg("do_statcov_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}

int
floattodouble(float *i_data, double *o_data, int N, carma_device *device){
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	floattodouble_krnl<<<grid , threads>>>(i_data,o_data,N);
	cutilCheckMsg("floattodouble_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}

int
doubletofloat(double *i_data, float *o_data, int N, carma_device *device){
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	doubletofloat_krnl<<<grid , threads>>>(i_data,o_data,N);
	cutilCheckMsg("floattodouble_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}

int
get_pupphase(float *o_data, float *i_data, int *indx_pup, int Nphi, carma_device *device){
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, Nphi, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);
	pupphase_krnl<<<grid , threads>>>(o_data,i_data,indx_pup,Nphi);
	cutilCheckMsg("pupphase_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}

int
compute_Hcor_gpu(float *o_data, int nrow, int ncol, float Fs, float gmin, float gmax, int delay, carma_device *device){
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, nrow*ncol, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);

	compute_Hcor_krnl<<<grid , threads>>>(o_data,nrow,ncol,Fs,1.0f/Fs,gmin,gmax,delay);
	cutilCheckMsg("compute_Hcor_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}

int
absnormfft(cuFloatComplex *idata, float *odata, int N, float norm, carma_device *device){
	int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, N, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);

	absnormfft_krnl<<<grid , threads>>>(idata,odata,N,norm);
	cutilCheckMsg("absnormfft_krnl<<<>>> execution failed\n");
	return EXIT_SUCCESS;
}
