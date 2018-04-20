#include <sutra_gamora.h>

__global__ void fillamplikrnl(cuFloatComplex *amplipup, float *phase, int *wherephase,
                              float scale, int Npts, int nx, int Nx, int puponly) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int nim;
  int nline;
  int ncol;

  while (tid < Npts) {
    nim = wherephase[tid];
    nline = nim / nx;
    ncol = nim - nline*nx;
    nim = ncol + nline*Nx;
    if(puponly == 1) {
      amplipup[nim].x = 1.0f;
      amplipup[nim].y = 0.0f;
    } else if(puponly == 0) {
      amplipup[nim].x = cosf(-scale * phase[tid]);
      amplipup[nim].y = sinf(-scale * phase[tid]);
    } else if(puponly == 2) {
      amplipup[nim].x = phase[tid];
      amplipup[nim].y = 0.0f;

    }
    tid += blockDim.x * gridDim.x;
  }
}

int fill_amplipup(cuFloatComplex *amplipup, float *phase, int *wherephase,
                  float scale, int Npts, int nx, int Nx, int puponly, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, Npts, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  fillamplikrnl<<<grid, threads>>>(amplipup, phase, wherephase, scale, Npts, nx, Nx, puponly);
  carmaCheckMsg("fillamplikrnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void cumulpsf_krnl(float *odata, cuFloatComplex *idata, int N) {
  cuFloatComplex cache;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    cache = idata[tid];
    odata[tid] += (cache.x * cache.x + cache.y * cache.y);
    tid += blockDim.x * gridDim.x;
  }
}

int cumulpsf(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  cumulpsf_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carmaCheckMsg("cumulpsf_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void abs2complex_krnl(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N) {
  cuFloatComplex idata;
  cuFloatComplex odata;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < N) {
    idata = d_idata[tid];
    odata.x = idata.x * idata.x + idata.y * idata.y;
    odata.y = 0.0f;
    d_odata[tid] = odata;
    tid += blockDim.x * gridDim.x;
  }
}
int abs2complex(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  abs2complex_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carmaCheckMsg("abs2complex_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void modulus2_krnl(float *d_odata, cuFloatComplex *d_idata, int N) {
  cuFloatComplex idata;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < N) {
    idata = d_idata[tid];
    d_odata[tid] = idata.x * idata.x + idata.y * idata.y;
    tid += blockDim.x * gridDim.x;
  }
}
int modulus2(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  modulus2_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carmaCheckMsg("modulus2_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void real_krnl(float *d_odata, cuFloatComplex *d_idata, int N) {
  cuFloatComplex cache;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < N) {
    cache = d_idata[tid];
    d_odata[tid] = cache.x;
    tid += blockDim.x * gridDim.x;
  }
}

int real(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  real_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carmaCheckMsg("real_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillmask_krnl(float *d_odata, float *d_idata, int N, int norm) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < N) {
    d_odata[tid] = d_idata[tid] < (1e-5*norm) ? 0.0f : 1.0f;
    tid += blockDim.x * gridDim.x;
  }
}

int fill_mask(float *d_odata, float *d_idata, int N, int norm, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  fillmask_krnl<<<grid, threads>>>(d_odata, d_idata, N, norm);
  carmaCheckMsg("fillmask_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void pow2_krnl(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N) {
  cuFloatComplex idata;
  cuFloatComplex odata;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < N) {
    idata = d_idata[tid];
    odata.x = idata.x * idata.x;
    odata.y = 0.0f;
    d_odata[tid] = odata;
    tid += blockDim.x * gridDim.x;
  }
}
int pow2(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  pow2_krnl<<<grid, threads>>>(d_odata, d_idata, N);
  carmaCheckMsg("pow2_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillterm1_krnl(float *d_odata, cuFloatComplex *d_idata,  cuFloatComplex *d_pupfft, int N) {
  cuFloatComplex idata;
  cuFloatComplex pupfft;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < N) {
    idata = d_idata[tid];
    pupfft = d_pupfft[tid];
    d_odata[tid] = idata.x * pupfft.x + idata.y * pupfft.y;
    tid += blockDim.x * gridDim.x;
  }
}

int fill_term1(float *d_odata, cuFloatComplex *d_idata, cuFloatComplex *d_pupfft, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  fillterm1_krnl<<<grid, threads>>>(d_odata, d_idata, d_pupfft, N);
  carmaCheckMsg("fillterm1_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void add2Dphi_krnl(cuFloatComplex *d_odata, float *d_term1,  float *d_term2, float e, int N) {
  cuFloatComplex cache;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < N) {
    cache = d_odata[tid];
    cache.x += 2*((d_term1[tid] - d_term2[tid]) * e);
    cache.y = 0.0f;
    d_odata[tid] = cache;
    tid += blockDim.x * gridDim.x;
  }
}

int add2Dphi(cuFloatComplex *d_odata, float *d_term1, float *d_term2, float e, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  add2Dphi_krnl<<<grid, threads>>>(d_odata, d_term1, d_term2, e, N);
  carmaCheckMsg("add2Dphi_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void computeOTFvii_krnl(float *d_otfVii, cuFloatComplex * d_Dphi, float *d_otftel,
                                   float *d_mask, float scale, int N) {

  cuFloatComplex Dphi;
  float tmp;
  float den;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while(tid < N) {
    Dphi = d_Dphi[tid];
    den = d_otftel[tid] > 1e-9 ? 1.0f/d_otftel[tid] : 0.0f;
    tmp = Dphi.x * den * d_mask[tid] * scale * scale;
    d_otfVii[tid] = expf(-0.5 * tmp) * d_mask[tid];

    tid += blockDim.x * gridDim.x;
  }
}
int computeOTFvii(float *d_otfVii, cuFloatComplex * d_Dphi, float *d_otftel, float *d_mask,
                  float scale, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  computeOTFvii_krnl<<<grid, threads>>>(d_otfVii, d_Dphi, d_otftel, d_mask,scale, N);
  carmaCheckMsg("computeOTFvii_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void ifftscale_krnl(cuFloatComplex *d_odata, float scale, int N) {
  cuFloatComplex cache;
  int tid = threadIdx.x + blockIdx.x*blockDim.x;

  while(tid < N) {
    cache = d_odata[tid];
    cache.x = cache.x * scale;
    cache.y = cache.y * scale;
    d_odata[tid] = cache;
    tid += blockDim.x * gridDim.x;
  }
}
int ifftscale(cuFloatComplex *d_odata, float scale, int N, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  ifftscale_krnl<<<grid, threads>>>(d_odata,scale, N);
  carmaCheckMsg("iffscale_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
