#include <sutra_wfs.h>
#include <sutra_ao_utils.h>
#include "carma_utils.cuh"

__global__ void camplipup_krnl(cuFloatComplex *amplipup, float *phase,
    float *offset, float *mask, float scale, int *istart, int *jstart,
    int *ivalid, int *jvalid, int nphase, int nphase2, int npup, int Nfft,
    int N, int offset_phase) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x + offset_phase;

  while (tid < N) {
    int nim = tid / nphase2;
    int idim = tid - nim * nphase2;

    int idimx = idim % nphase; // nphase : size of the phase support in subaps
    int idimy = idim / nphase;

    int idphase = idimx + idimy * npup + istart[ivalid[nim]]
        + jstart[jvalid[nim]] * npup;
    // npup : size of the input phase screen

    int idx = idimx + idimy * Nfft + nim * Nfft * Nfft;

    amplipup[idx].x = (cosf(phase[idphase] * scale - offset[idim]))
        * mask[idphase];
    amplipup[idx].y = (sinf(phase[idphase] * scale - offset[idim]))
        * mask[idphase];
    tid += blockDim.x * gridDim.x;
  }
}

int fillcamplipup(cuFloatComplex *amplipup, float *phase, float *offset,
    float *mask, float scale, int *istart, int *jstart, int *ivalid,
    int *jvalid, int nphase, int npup, int Nfft, int Ntot, carma_device *device, int offset_phase=0)
// here amplipup is a cube of data of size nfft x nfft x nsubap
// phase is an array of size pupdiam x pupdiam
// offset is an array of size pdiam x pdiam
// mask is an array of size pupdiam x pupdiam
// number of thread required : pdiam x pdiam x nsubap
    {
  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, Ntot, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);
  /*
   int nthreads = 0,nblocks = 0;
   getNumBlocksAndThreads(device, Ntot, nblocks, nthreads);

   dim3 grid(nblocks), threads(nthreads);
   */

  int nphase2 = nphase * nphase;

  camplipup_krnl<<<grid, threads>>>(amplipup, phase, offset, mask, scale,
      istart, jstart, ivalid, jvalid, nphase, nphase2, npup, Nfft, Ntot, 0);//offset_phase);
  carmaCheckMsg("camplipup_krnl<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}


__global__ void bimg_krnl(float *bimage, float *bcube, int npix, int npix2,
    int nsub, int *ivalid, int *jvalid, float alpha, int N) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
   Npix = npix x npix
   */
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    int nim = tid / npix2;
    int tidim = tid - nim * npix2;
    int xim = tidim % npix;
    int yim = tidim / npix;

    int idbin = xim + yim * nsub + ivalid[nim] * npix
        + jvalid[nim] * npix * nsub;
    bimage[idbin] = alpha * bimage[idbin] + bcube[tid];

    tid += blockDim.x * gridDim.x;
  }
}

int fillbinimg(float *bimage, float *bcube, int npix, int nsub, int Nsub,
    int *ivalid, int *jvalid, bool add, carma_device *device) {
  int Npix = npix * npix;
  int N = Npix * nsub;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  float alpha;
  if (add)
    alpha = 1.0f;
  else
    alpha = 0.0f;

  bimg_krnl<<<grid, threads>>>(bimage, bcube, npix, Npix, Nsub, ivalid, jvalid,
      alpha, N);

  carmaCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void bimg_krnl_async(float *bimage, float *bcube, int npix,
    int npix2, int nsub, int *ivalid, int *jvalid, float alpha, int N,
    int idstart) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
   Npix = npix x npix
   */
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  tid += idstart;

  while (tid < N) {
    int nim = tid / npix2;
    int tidim = tid - nim * npix2;
    int xim = tidim % npix;
    int yim = tidim / npix;

    int idbin = xim + yim * nsub + ivalid[nim] * npix
        + jvalid[nim] * npix * nsub;
    bimage[idbin] = alpha * bimage[idbin] + bcube[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int fillbinimg_async(carma_host_obj<float> *image_telemetry, float *bimage,
    float *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid,
    int nim, bool add, carma_device *device) {
  float *hdata = image_telemetry->getData();
  int nstreams = image_telemetry->get_nbStreams();

  int Npix = npix * npix;
  int N = Npix * nsub;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  // here nstreams should be : final image size / npix
  dim3 threads(nthreads);
  dim3 grid(N / (nstreams * threads.x));
  float alpha;
  if (add)
    alpha = 1.0f;
  else
    alpha = 0.0f;

  // asynchronously launch nstreams kernels, each operating on its own portion of data
  for (int i = 0; i < nstreams; i++) {
    bimg_krnl_async<<<grid, threads, 0, image_telemetry->get_cudaStream_t(i)>>>(
        bimage, bcube, npix, Npix, Nsub, ivalid, jvalid, alpha, N,
        i * N / nstreams);
    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x will only
    //   commence executing when all previous CUDA calls in stream x have completed
    cudaMemcpyAsync(&(hdata[i * nim / nstreams]), &(bimage[i * nim / nstreams]),
        sizeof(float) * nim / nstreams, cudaMemcpyDeviceToHost,
        image_telemetry->get_cudaStream_t(i));
  }
  //cudaStreamSynchronize(image_telemetry->get_cudaStream_t(nstreams-1));
  carmaCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int fillbinimg_async(carma_streams *streams, carma_obj<float> *bimage,
    carma_obj<float> *bcube, int npix, int nsub, int Nsub, int *ivalid,
    int *jvalid, bool add, carma_device *device) {
  float *g_image = bimage->getData();
  float *g_cube = bcube->getData();
  int nstreams = streams->get_nbStreams();

  int Npix = npix * npix;
  int N = Npix * nsub;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  // here nstreams should be : final image size / npix
  dim3 threads(nthreads);
  dim3 grid(N / (nstreams * threads.x));

  float alpha;
  if (add)
    alpha = 1.0f;
  else
    alpha = 0.0f;

  // asynchronously launch nstreams kernels, each operating on its own portion of data
  for (int i = 0; i < nstreams; i++) {

    bimg_krnl_async<<<grid, threads, 0, streams->get_stream(i)>>>(g_image,
        g_cube, npix, Npix, Nsub, ivalid, jvalid, alpha, N, i * N / nstreams);
    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x will only
    //   commence executing when all previous CUDA calls in stream x have completed
  }

  carmaCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillbincube_krnl(float *bcube, cuFloatComplex *hrimage,
    int *indxpix, int Nfft, int Npix, int Nrebin, int N) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
   Npix = npix x npix
   */
  int npix, nsubap, nrebin;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  //int nim = tid / Npix;
  //int tidim = tid - nim * Npix;
  //int xim = tidim % 20;
  //int yim = tidim / 20;
  while (tid < N) {
    nsubap = tid / Npix;
    npix = tid % Npix;
   // if (xim>=6 && xim<14 && yim>=6 && yim<14){
    for (int i = 0; i < Nrebin; i++) {
      nrebin = indxpix[i + npix * Nrebin];
      bcube[tid] += hrimage[nrebin + Nfft * nsubap].x;
    }
   // }
   // else bcube[tid] = 0.;
    tid += blockDim.x * gridDim.x;
  }
}

int fillbincube(float *bcube, cuFloatComplex *hrimage, int *indxpix, int Nfft,
    int Npix, int Nrebin, int Nsub, carma_device *device) {
  int N = Npix * Nsub;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  fillbincube_krnl<<<grid, threads>>>(bcube, hrimage, indxpix, Nfft, Npix,
      Nrebin, N);
  carmaCheckMsg("fillbincube_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillbincube_krnl_async(float *bcube, cuFloatComplex *hrimage,
    int *indxpix, int Nfft, int Npix, int Nrebin, int N, int idstart) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
   Npix = npix x npix
   */
  int npix, nsubap, nrebin;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  tid += idstart;

  while (tid < N) {
    nsubap = tid / Npix;
    npix = tid % Npix;
    for (int i = 0; i < Nrebin; i++) {
      nrebin = indxpix[i + npix * Nrebin];
      bcube[tid] += hrimage[nrebin + Nfft * nsubap].x;
    }
    tid += blockDim.x * gridDim.x;
  }
}

int fillbincube_async(carma_streams *streams, float *bcube,
    cuFloatComplex *hrimage, int *indxpix, int Nfft, int Npix, int Nrebin,
    int Nsub, carma_device *device) {
  int N = Npix * Nsub;
  int nthreads = 0, nblocks = 0;
  int nstreams = streams->get_nbStreams();
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  // asynchronously launch nstreams kernels, each operating on its own portion of data
  for (int i = 0; i < nstreams; i++) {
    fillbincube_krnl_async<<<grid, threads, 0, streams->get_stream(i)>>>(bcube,
        hrimage, indxpix, Nfft, Npix, Nrebin, N, i * N / nstreams);
    carmaCheckMsg("fillbincubeasync_kernel<<<>>> execution failed\n");
  }
  return EXIT_SUCCESS;
}

__global__ void indexfill_krnl(cuFloatComplex *odata, cuFloatComplex *idata,
    int *indx, int ntot, int Ntot, int N) {
  int nim, idim;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / ntot;
    idim = tid - (nim * ntot);
    odata[indx[idim] + (nim * Ntot)].x = idata[tid].x;
    tid += blockDim.x * gridDim.x;
  }
}

int indexfill(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int *indx,
    int nx, int Nx, int N, carma_device *device) {

  int ntot = nx * nx;
  int Ntot = Nx * Nx;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  indexfill_krnl<<<grid, threads>>>(d_odata, d_idata, indx, ntot, Ntot, N);

  carmaCheckMsg("indexfill_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void conv_krnl(cuFloatComplex *odata, cuFloatComplex *idata, int N,
    int n) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  cuFloatComplex tmp;
  int nim;
  int tidim;
  //int nim = tid / n;
  //int tidim = tid - nim * n;

  while (tid < N) {
    nim = tid / n;
    tidim = tid - nim * n;
    tmp.x = idata[tidim].x * odata[tid].x - idata[tidim].y * odata[tid].y;
    tmp.y = idata[tidim].y * odata[tid].x + idata[tidim].x * odata[tid].y;
    odata[tid] = tmp;
    tid += blockDim.x * gridDim.x;
  }
}

int convolve_cube(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
    int n, carma_device *device) {
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  conv_krnl<<<grid, threads>>>(d_odata, d_idata, N, n);

  carmaCheckMsg("conv_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

template<class T>
__device__ void red_krnl(T *sdata, int size, int n) {
  if (!((size & (size - 1)) == 0)) {
    unsigned int s;
    if (size % 2 != 0)
      s = size / 2 + 1;
    else
      s = size / 2;
    unsigned int s_old = size;
    while (s > 0) {
      if ((n < s) && (n + s < s_old)) {
        sdata[n] += sdata[n + s];
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2 * s < s_old) && (s != 0))
        s += 1;
    }
  } else {
    // do reduction in shared mem
    for (unsigned int s = size / 2; s > 0; s >>= 1) {
      if (n < s) {
        sdata[n] += sdata[n + s];
      }
      __syncthreads();
    }
  }
}

template<class T>
__global__ void reduce2(T *g_idata, T *g_odata, unsigned int n, unsigned int nelem_thread =1) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = 0;
  for (int cc=0;cc < nelem_thread; cc++) {
    int idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < n)
      sdata[tid] +=g_idata[idim];
    else
      sdata[tid] += 0;
  }

  //sdata[tid] = (i < n) ? g_idata[i] : 0;

  __syncthreads();
  red_krnl(sdata, blockDim.x, tid);
  __syncthreads();

//  if (tid == 0)
//    printf("blockIdx.x %d \n", blockIdx.x);
//  __syncthreads();
//
  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = sdata[tid];
}

template<class T>
__global__ void reduce2(T *g_idata, T *g_odata, T thresh, unsigned int n, unsigned int nelem_thread =1) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = 0;
  for (int cc=0;cc < nelem_thread; cc++) {
    int idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < n) { 
      if (g_idata[idim] > thresh)  sdata[tid] += g_idata[idim];
      else sdata[tid] += 0;
    } else
      sdata[tid] += 0;
  }
  /*
  if (i < n)
    sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] : 0;
  */
  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);
  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = sdata[0];
}

template<class T>
__global__ void reduce2_new(T *g_idata, T *g_odata, T thresh, unsigned int n, unsigned int nelem_thread =1) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = 0;
  for (int cc=0;cc < nelem_thread; cc++) {
    int idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < n) {
      if (g_idata[idim] <= thresh)
    	  g_idata[idim] = 0;
	  sdata[tid] += g_idata[idim];
    } else
      sdata[tid] += 0;
  }
  /*
  if (i < n)
    sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] : 0;
  */
  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);
  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = sdata[0];
}

template<class T>
__global__ void reduce2(T *g_idata, T *g_odata, T *weights, unsigned int n, unsigned int nelem_thread=1) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  //unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  sdata[tid] = 0;
  for (int cc=0;cc < nelem_thread; cc++) {
    int idim = tid * nelem_thread + cc + (blockDim.x * nelem_thread) * blockIdx.x;
    if (idim < n) 
      sdata[tid] +=g_idata[idim] * weights[idim];
    else
      sdata[tid] += 0;
  }

//   sdata[tid] = (i < n) ? g_idata[i] * weights[i] : 0;

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = sdata[tid];
}

template<class T>
void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata, carma_device *device) {

	int maxThreads = device->get_properties().maxThreadsPerBlock;
	unsigned int nelem_thread = 1;
	while((threads/nelem_thread > maxThreads) || (threads % nelem_thread != 0)){
		nelem_thread++;
	}

  threads /= nelem_thread;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = threads * sizeof(T);
  reduce2<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, (unsigned int)size,nelem_thread);

  carmaCheckMsg("reduce2_kernel<<<>>> execution failed\n");
}
template void
subap_reduce<float>(int size, int threads, int blocks, float *d_idata,
    float *d_odata, carma_device *device);
template void
subap_reduce<double>(int size, int threads, int blocks, double *d_idata,
    double *d_odata, carma_device *device);

template<class T>
__global__ void reduce2_async(T *g_idata, T *g_odata, unsigned int n,
    int stream_offset) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  i += stream_offset * blockDim.x;
  //tid+=idstream*nblocks/nstreams;

  sdata[tid] = (i < n) ? g_idata[i] : 0;

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);
  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x + stream_offset] = sdata[0];
}

template<class T>
void subap_reduce_async(int threads, int blocks, carma_streams *streams,
    T *d_idata, T *d_odata) {
  int nstreams = streams->get_nbStreams();
  dim3 dimBlock(threads);
  dim3 dimGrid(blocks / nstreams);

  int nbelem = threads * blocks;
  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = threads * sizeof(T);
  for (int i = 0; i < nstreams; i++) {
    reduce2_async<T> <<<dimGrid, dimBlock, smemSize, (*streams)[i]>>>(d_idata,
        d_odata, nbelem, i * blocks / nstreams);
  }

  carmaCheckMsg("reduce_kernel<<<>>> execution failed\n");
}
template void
subap_reduce_async<float>(int threads, int blocks, carma_streams *streams,
    float *d_idata, float *d_odata);
template void
subap_reduce_async<double>(int threads, int blocks, carma_streams *streams,
    double *d_idata, double *d_odata);

template<class T>
void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata,
    T thresh, carma_device *device) {

  int maxThreads = device->get_properties().maxThreadsPerBlock;
  unsigned int nelem_thread = 1;
  while((threads/nelem_thread > maxThreads) || (threads % nelem_thread != 0)){
		nelem_thread++;
  }

  dim3 dimBlock(threads / nelem_thread, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =  threads * sizeof(T);
  reduce2<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, thresh,
      size,nelem_thread);
  carmaCheckMsg("reduce_kernel<<<>>> execution failed\n");
}
template void
subap_reduce<float>(int size, int threads, int blocks, float *d_idata,
    float *d_odata, float thresh, carma_device *device);

template void
subap_reduce<double>(int size, int threads, int blocks, double *d_idata,
    double *d_odata, double thresh, carma_device *device);

template<class T>
void subap_reduce_new(int size, int threads, int blocks, T *d_idata, T *d_odata,
    T thresh, carma_device *device) {

  int maxThreads = device->get_properties().maxThreadsPerBlock;
  unsigned int nelem_thread = 1;
  while((threads/nelem_thread > maxThreads) || (threads % nelem_thread != 0)){
		nelem_thread++;
  }

  dim3 dimBlock(threads / nelem_thread, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =  threads * sizeof(T);
  reduce2_new<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, thresh,
      size,nelem_thread);
  carmaCheckMsg("reduce_kernel<<<>>> execution failed\n");
}
template void
subap_reduce_new<float>(int size, int threads, int blocks, float *d_idata,
    float *d_odata, float thresh, carma_device *device);

template void
subap_reduce_new<double>(int size, int threads, int blocks, double *d_idata,
    double *d_odata, double thresh, carma_device *device);

template<class T>
void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata,
    T *weights, carma_device *device) {

	int maxThreads = device->get_properties().maxThreadsPerBlock;
	unsigned int nelem_thread = 1;
	while((threads/nelem_thread > maxThreads) || (threads % nelem_thread != 0)){
		nelem_thread++;
	}

  dim3 dimBlock(threads / nelem_thread, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = threads * sizeof(T);
  reduce2<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, weights,
      size,nelem_thread);
  carmaCheckMsg("reduce_kernel<<<>>> execution failed\n");
}
template void
subap_reduce<float>(int size, int threads, int blocks, float *d_idata,
    float *d_odata, float *weights, carma_device *device);

template void
subap_reduce<double>(int size, int threads, int blocks, double *d_idata,
    double *d_odata, double *weights, carma_device *device);

template<class T>
__global__ void reduce_phasex(T *g_idata, T *g_odata, int *indx, unsigned int n,
    T alpha) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * n * n;

  sdata[tid] = g_idata[indx[i + tid * n + n - 1]] - g_idata[indx[i + tid * n]];

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0) {
    g_odata[blockIdx.x] = sdata[0] / n * alpha;
  }
}

template<class T>
__global__ void reduce_phasey(T *g_idata, T *g_odata, int *indx, unsigned int n,
    T alpha)
// FIXME
// full parallelization would require to use 4 x threads
// instead of doing 4 operations per threads
    {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * n * n + threadIdx.x;

  sdata[tid] = g_idata[indx[i + (n - 1) * n]] - g_idata[indx[i]];

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0) {
    g_odata[blockIdx.x] = sdata[0] / n * alpha;
  }
}

template<class T>
__global__ void derive_phasex(T *g_idata, T *g_odata, int *indx, T *mask,
    T alpha, unsigned int n, unsigned int N, float *fluxPerSub) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    if (tid % n == 0) { //start of a column
      sdata[tid] = (g_idata[indx[i + 1]] - g_idata[indx[i]]) * mask[indx[i]];
    } else {
      if ((tid + 1) % n == 0) { //end of a column
        sdata[tid] = (g_idata[indx[i]] - g_idata[indx[i - 1]]) * mask[indx[i]];
      } else
        sdata[tid] = (g_idata[indx[i + 1]] - g_idata[indx[i - 1]])
            * mask[indx[i]];
    }
  }

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = sdata[0] / n * alpha / fluxPerSub[blockIdx.x] / 2;
}

template<class T>
__global__ void derive_phasey(T *g_idata, T *g_odata, int *indx, T *mask,
    T alpha, unsigned int n, unsigned int N, float *fluxPerSub) {
  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < N) {
    if (tid < n) { //start of a column
      sdata[tid] = (g_idata[indx[i + n]] - g_idata[indx[i]]) * mask[indx[i]];
    } else {
      if (tid >= n * (n - 1)) { //end of a column
        sdata[tid] = (g_idata[indx[i]] - g_idata[indx[i - n]]) * mask[indx[i]];
      } else
        sdata[tid] = (g_idata[indx[i + n]] - g_idata[indx[i - n]])
            * mask[indx[i]];
    }
  }

  __syncthreads();

  red_krnl(sdata, blockDim.x, tid);

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = sdata[0] / n / 2 * alpha / fluxPerSub[blockIdx.x];
}

template<class T>
void phase_reduce(int threads, int blocks, T *d_idata, T *d_odata, int *indx,
    T alpha) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  reduce_phasex<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, indx,
      threads, alpha);

  carmaCheckMsg("reduce_phasex_kernel<<<>>> execution failed\n");

  reduce_phasey<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
      &(d_odata[blocks]), indx, threads, alpha);

  carmaCheckMsg("reduce_phasey_kernel<<<>>> execution failed\n");
}

template void
phase_reduce<float>(int threads, int blocks, float *d_idata, float *d_odata,
    int *indx, float alpha);

template void
phase_reduce<double>(int threads, int blocks, double *d_idata, double *d_odata,
    int *indx, double alpha);

template<class T>
void phase_derive(int size, int threads, int blocks, int n, T *d_idata,
    T *d_odata, int *indx, T *mask, T alpha, float *fluxPerSub) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize =
      (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
  derive_phasex<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, indx,
      mask, alpha, n, size, fluxPerSub);

  carmaCheckMsg("phase_derivex_kernel<<<>>> execution failed\n");

  derive_phasey<T> <<<dimGrid, dimBlock, smemSize>>>(d_idata,
      &(d_odata[blocks]), indx, mask, alpha, n, size, fluxPerSub);

  carmaCheckMsg("phase_derivey_kernel<<<>>> execution failed\n");
}

template void
phase_derive<float>(int size, int threads, int blocks, int n, float *d_idata,
    float *d_odata, int *indx, float *mask, float alpha, float *fluxPerSub);

template void
phase_derive<double>(int size, int threads, int blocks, int n, double *d_idata,
    double *d_odata, int *indx, double *mask, double alpha, float *fluxPerSub);

template<class Tout, class Tin>
__global__ void pyrgetpup_krnl(Tout *g_odata, Tin *g_idata, Tout *offsets,
    Tin *pup, float lambda, unsigned int n) {
  // roll( pup * exp(i*phase) ) * offsets

  Tout *sdata = SharedMemory<Tout>();
  const float PI = 3.1415926535897932384626433;
  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  int x, y, xx, yy, i2;

  if (i < n * n / 2) {
    x = i % n;
    y = i / n;

    xx = (x + n / 2) % n;
    yy = (y + n / 2) % n;
    i2 = xx + yy * n;
    
    float mic2rad =(2*PI/lambda);
        
    sdata[2 * tid].x = cosf(g_idata[i]*mic2rad) * pup[i];
    sdata[2 * tid].y = sinf(g_idata[i]*mic2rad) * pup[i];

    sdata[2 * tid + 1].x = cosf(g_idata[i2]*mic2rad) * pup[i2];
    sdata[2 * tid + 1].y = sinf(g_idata[i2]*mic2rad) * pup[i2];
  }
  __syncthreads();

  if (i < n * n / 2) {

    g_odata[i].x = sdata[2 * tid + 1].x * offsets[i].x
        - sdata[2 * tid + 1].y * offsets[i].y;
    g_odata[i].y = sdata[2 * tid + 1].x * offsets[i].y
        + sdata[2 * tid + 1].y * offsets[i].x;

    g_odata[i2].x = sdata[2 * tid].x * offsets[i2].x
        - sdata[2 * tid].y * offsets[i2].y;
    g_odata[i2].y = sdata[2 * tid].x * offsets[i2].y
        + sdata[2 * tid].y * offsets[i2].x;
  }
}

template<class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tout *d_offsets, Tin *d_pup,
    int np, float lambda, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, np * np / 2, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  int smemSize = 2 * nThreads * sizeof(Tout);
  pyrgetpup_krnl<Tout, Tin> <<<grid, threads, smemSize>>>(d_odata, d_idata,
      d_offsets, d_pup, lambda, np);

  carmaCheckMsg("pyrgetpup_kernel<<<>>> execution failed\n");
}
template void
pyr_getpup<cuFloatComplex, float>(cuFloatComplex *d_odata, float *d_idata,
    cuFloatComplex *d_offsets, float *d_pup, int np, float lambda, carma_device *device);
template void
pyr_getpup<cuDoubleComplex, double>(cuDoubleComplex *d_odata, double *d_idata,
    cuDoubleComplex *d_offsets, double *d_pup, int np, float lambda, carma_device *device);

/*
////////////////////////////////////////////////////////////////////////:
New version of pyrgetpup for hr version of pyramid model
////////////////////////////////////////////////////////////////////////
 */


template<class Tout, class Tin>
__global__ void pyrgetpup_krnl(Tout *g_odata, Tin *g_idata,
			       Tin *pup, float lambda, float cx, float cy, unsigned int n, unsigned int N) {
  const float PI = 3.1415926535897932384626433;
  // load shared mem
  //const unsigned int tid = threadIdx.x;
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < n * n) {
    const float x = i % n;
    const float y = i / n;

    const float xx = (x - n / 2.0f);
    const float yy = (y - n / 2.0f);
        
    const float mic2rad =(2*PI/lambda);

    const float phi_modu =  cx * xx + cy * yy;

    const float offs =  (N - n)/2.0f;

    const int i2 = x + offs + (y + offs)* N;
 
    g_odata[i2].x = cosf(g_idata[i]*mic2rad + phi_modu) * pup[i];
    g_odata[i2].y = sinf(g_idata[i]*mic2rad + phi_modu) * pup[i];

  }
}

template<class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tin *d_pup,
    int np, int N, float lambda, float cx, float cy, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, np * np, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  pyrgetpup_krnl<Tout, Tin> <<<grid, threads>>>(d_odata, d_idata,
						d_pup, lambda, cx, cy, np, N);

  carmaCheckMsg("pyrgetpup_kernel<<<>>> execution failed\n");
}
template void
pyr_getpup<cuFloatComplex, float>(cuFloatComplex *d_odata, float *d_idata,
    float *d_pup, int np, int N, float lambda, float cx, float cy, carma_device *device);
template void
pyr_getpup<cuDoubleComplex, double>(cuDoubleComplex *d_odata, double *d_idata,
    double *d_pup, int np, int N, float lambda, float cx, float cy, carma_device *device);



template<class T>
__global__ void rollmod_krnl(T *g_odata, T *g_idata, T *g_mask, float cx, float cy,
    unsigned int n, unsigned int ns, unsigned int nim) {
  // roll( pup * exp(i*phase) ) * offsets

  // load shared mem
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  int xx, yy, i2;

  if (i < ns * ns * nim) {
    int n_im = i / (ns * ns);
    int tidim = i - n_im * ns * ns;

    xx = tidim % ns;
    yy = tidim / ns;

    int xx2 = 0;
    int yy2 = 0;

    //cx *= 0;
    //cy *= 0;

    // here we re-order the quadrant to simulate a roll :
    // upper right is 1, upper left is 2, bottom right is 3 and bottom left is 4
    if (n_im == 0) {
      // first image : upper right
      xx2 = (xx + (n - ns) - cx);
      yy2 = (yy + (n - ns) - cy);
    }
    if (n_im == 1) {
      // second image : upper left
      xx2 = (xx - cx);
      yy2 = (yy + (n - ns) - cy);
    }
    if (n_im == 2) {
      // third image : lower right
      xx2 = (xx + (n - ns) - cx);
      yy2 = (yy - cy);
    }
    if (n_im == 3) {
      // fourth image : lower left
      xx2 = (xx - cx);
      yy2 = (yy - cy);
    }

    if (xx2 < 0)
      xx2 = xx2 + n;
    else
      xx2 = xx2 % n;

    if (yy2 < 0)
      yy2 = yy2 + n;
    else
      yy2 = yy2 % n;

    i2 = xx2 + yy2 * n;

    if (i2 < n * n) {
      T tmp1, tmp2;
      tmp1 = g_idata[i2];
      tmp2 = g_mask[tidim];
      g_odata[i].x = tmp1.x * tmp2.x - tmp1.y * tmp2.y;
      g_odata[i].y = tmp1.x * tmp2.y + tmp1.y * tmp2.x;
    }
  }
}

template<class T>
void pyr_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int np,
    int ns, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, ns * ns * 4, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  rollmod_krnl<T> <<<grid, threads>>>(d_odata, d_idata, d_mask, cx, cy, np,
      ns, 4);

  carmaCheckMsg("rollmod_kernel<<<>>> execution failed\n");
}
template void
pyr_rollmod<cuFloatComplex>(cuFloatComplex *d_odata, cuFloatComplex*d_idata,
    cuFloatComplex*d_mask, float cx, float cy, int np, int ns, carma_device *device);
template void
pyr_rollmod<cuDoubleComplex>(cuDoubleComplex *d_odata, cuDoubleComplex *d_idata,
    cuDoubleComplex *d_mask, float cx, float cy, int np, int ns, carma_device *device);

//////////////////////////////////////////////////////////////
// ADDING PYR_ROLLMOD MODIFIED FOR ROOF-PRISM: ROOF_ROOLMOD //
//////////////////////////////////////////////////////////////

template<class T>
__global__ void roof_rollmod_krnl(T *g_odata, T *g_idata, T *g_mask, float cx, float cy,
    unsigned int n, unsigned int ns, unsigned int nim) {
  // roll( pup * exp(i*phase) ) * offsets

  // load shared mem
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  int xx, yy, i2;

  if (i < ns * ns * nim) {
    int n_im = i / (ns * ns);
    int tidim = i - n_im * ns * ns;

    xx = tidim % ns;
    yy = tidim / ns;

    int xx2 = 0;
    int yy2 = 0;

    //cx *= 0;
    //cy *= 0;

    // here we re-order the half-planes to simulate a roll :
    // top is 1 and bottom is 2 ; left is 3 and right is 4.
    if (n_im == 0) {
      // first image : top
      xx2 = (xx + (n - ns/2));
      yy2 = (yy + (n - ns) - cy);
    }
    if (n_im == 1) {
      // second image : bottom
      xx2 = (xx + (n - ns/2));
      yy2 = (yy - cy);
    }
    if (n_im == 2) {
      // third image : left
      xx2 = (xx + (n - ns) - cx);
      yy2 = (yy + (n - ns/2)); 
    }
    if (n_im == 3) {
      // fourth image : right
      xx2 = (xx - cx);
      yy2 = (yy + (n - ns/2));
    }

    if (xx2 < 0)
      xx2 = xx2 + n;
    else
      xx2 = xx2 % n;

    if (yy2 < 0)
      yy2 = yy2 + n;
    else
      yy2 = yy2 % n;

    i2 = xx2 + yy2 * n;

    if (i2 < n * n) {
      T tmp1, tmp2;
      tmp1 = g_idata[i2];
      tmp2 = g_mask[tidim];
      g_odata[i].x = tmp1.x * tmp2.x - tmp1.y * tmp2.y;
      g_odata[i].y = tmp1.x * tmp2.y + tmp1.y * tmp2.x;
    }
  }
}

template<class T>
void roof_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int np,
    int ns, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, ns * ns * 4, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  roof_rollmod_krnl<T> <<<grid, threads>>>(d_odata, d_idata, d_mask, cx, cy, np,
      ns, 4);

  carmaCheckMsg("roof_rollmod_kernel<<<>>> execution failed\n");
}
template void
roof_rollmod<cuFloatComplex>(cuFloatComplex *d_odata, cuFloatComplex*d_idata,
    cuFloatComplex*d_mask, float cx, float cy, int np, int ns, carma_device *device);
template void
roof_rollmod<cuDoubleComplex>(cuDoubleComplex *d_odata, cuDoubleComplex *d_idata,
    cuDoubleComplex *d_mask, float cx, float cy, int np, int ns, carma_device *device);

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template<class T>
__global__ void fillbinpyr_krnl(T *g_odata, T *g_idata, unsigned int nrebin,
    unsigned int n, unsigned int ns, unsigned int nim) {
  // bin2d(hrimg)

  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  int x, y;

  if (i < ns * ns * nim) {
    int n_im = i / (ns * ns);
    int tidim = i - n_im * (ns * ns);
    x = tidim % ns;
    y = tidim / ns;
    int xim = x * nrebin;
    int yim = y * nrebin;

    sdata[tid] = 0;

    for (int cpty = 0; cpty < nrebin; cpty++) {
      for (int cptx = 0; cptx < nrebin; cptx++) {
        int tidim2 = (xim + cptx) + (yim + cpty) * n + n_im * n * n;
        sdata[tid] += g_idata[tidim2];
      }
    }
//    sdata[tid] /= (nrebin * nrebin);
  }
  __syncthreads();

  if (i < ns * ns * nim) {
    g_odata[i] = sdata[tid];
  }
}

template<class T>
void pyr_fillbin(T *d_odata, T *d_idata, int nrebin, int np, int ns, int nim,
    carma_device *device) {
  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, ns * ns * nim, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  int smemSize = nThreads * sizeof(T);
  fillbinpyr_krnl<T> <<<grid, threads, smemSize>>>(d_odata, d_idata, nrebin,
      np, ns, nim);

  carmaCheckMsg("pyrgetpup_kernel<<<>>> execution failed\n");
}
template void
pyr_fillbin<float>(float *d_odata, float *d_idata, int nrebin, int np, int ns,
    int nim, carma_device *device);
template void
pyr_fillbin<double>(double *d_odata, double *d_idata, int nrebin, int np,
    int ns, int nim, carma_device *device);

template<class T>
__global__ void pyr_bimg_krnl(T *bimage, const T *bcube,
                              const int nxsub, const float alpha, const int N) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
   Npix = npix x npix
   */
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  const int nxsub2 = nxsub*nxsub;
  const int nximg = (2*nxsub+3);

  while (tid < N) {
    const int nim = tid / nxsub2;
    const int tidim = tid - nim * nxsub2;
    const int yim = tidim / nxsub;
    const int xim = tidim - yim * nxsub;

    int offset;
    switch(nim){
      case 0:
        offset = 1;
        break;
      case 1:
        offset = 2 + nxsub;
        break;
      case 2:
        offset = 1 + (nximg*(1+nxsub));
        break;
      default:
        offset = 1 + (nximg+1)*(1+nxsub);
        break;
    }
    const int idbin = offset + xim + (yim+1) * nximg;
    bimage[idbin] = alpha * bimage[idbin] + bcube[tid];
    tid += blockDim.x * gridDim.x;
  }
}

template<class T>
int pyr_fillbinimg(T *bimage, const T *bcube, const int nxsub,
                   const bool add, carma_device *device) {
  const int N = 4*nxsub*nxsub;
  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  T alpha;
  if (add)
    alpha = 1.0;
  else
    alpha = 0.0;

  pyr_bimg_krnl<<<grid, threads>>>(bimage, bcube, nxsub, alpha, N);

  carmaCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template
int pyr_fillbinimg<float>(float *bimage, const float *bcube, const int nxsub,
                          const bool add, carma_device *device);
template
int pyr_fillbinimg<double>(double *bimage, const double *bcube, const int nxsub,
                           const bool add, carma_device *device);

template<class T>
__global__ void pyr_bimg_krnl(T *oimage, const T *image, const int n, 
                              const int rebin, const float alpha, const int N) {
  /*
   indx is an array nrebin^2 * npix^2
   it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
   Npix = npix x npix
   */
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < n * n) {
    const int x = tid % n;
    const int y = tid / n;
    const int xx = x * rebin;
    const int yy = y * rebin;
    oimage[tid] *= alpha;
    for (int ii = 0; ii <  rebin; ii++) {
      for (int jj = 0; jj <  rebin; jj++) {
	const int xim = xx + ii;
	const int yim = yy + jj;
	const int iim = xim + yim * N;
	oimage[tid] = oimage[tid] + image[iim];
      }
    }
    oimage[tid] /= (rebin * rebin);
    tid += blockDim.x * gridDim.x;
  }
}

template<class T>
int pyr_fillbinimg(T *oimage, const T *image, const int n, const int N,
                   const int rebin, const bool add, carma_device *device) {

  int nthreads = 0, nblocks = 0;
  getNumBlocksAndThreads(device, n*n, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  T alpha;
  if (add)
    alpha = 1.0;
  else
    alpha = 0.0;

  pyr_bimg_krnl<<<grid, threads>>>(oimage, image, n, rebin, alpha, N);

  carmaCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template
int pyr_fillbinimg<float>(float *oimage, const float *image, const int n, const int N,
			  const int rebin, const bool add, carma_device *device);
template
int pyr_fillbinimg<double>(double *oimage, const double *image, const int n, const int N,
			  const int rebin, const bool add, carma_device *device);

//////////////////////////////////////////////////////////////
// ADDING PYR_FILLBIN MODIFIED FOR ROOF-PRISM: ROOF_FILLBIN //
//////////////////////////////////////////////////////////////

template<class T>
__global__ void fillbinroof_krnl(T *g_odata, T *g_idata, unsigned int nrebin,
    unsigned int n, unsigned int ns, unsigned int nim) {
  // bin2d(hrimg)

  T *sdata = SharedMemory<T>();

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  int x, y;

  if (i < ns * ns * nim) {
    int n_im = i / (ns * ns);
    int tidim = i - n_im * ns * ns;

    /// Reprendre ici...

    x = tidim % ns;
    y = tidim / ns;
    int xim = x * nrebin;
    int yim = y * nrebin;

    sdata[tid] = 0;

    for (int cpty = 0; cpty < nrebin; cpty++) {
      for (int cptx = 0; cptx < nrebin; cptx++) {
        int tidim2 = (xim + cptx) + (yim + cpty) * n + n_im * n * n;
        sdata[tid] += g_idata[tidim2];
      }
    }
    sdata[tid] /= (nrebin * nrebin);
  }
  __syncthreads();

  if (i < ns * ns * nim) {
    g_odata[i] = sdata[tid];
  }
}

template<class T>
void roof_fillbin(T *d_odata, T *d_idata, int nrebin, int np, int ns, int nim,
    carma_device *device) {
  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, ns * ns * nim, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  int smemSize = nThreads * sizeof(T);
  fillbinroof_krnl<T> <<<grid, threads, smemSize>>>(d_odata, d_idata, nrebin,
      np, ns, nim);

  carmaCheckMsg("pyrgetpup_kernel<<<>>> execution failed\n");
}
template void
roof_fillbin<float>(float *d_odata, float *d_idata, int nrebin, int np, int ns,
    int nim, carma_device *device);
template void
roof_fillbin<double>(double *d_odata, double *d_idata, int nrebin, int np,
    int ns, int nim, carma_device *device);

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template<class Tout, class Tin>
__global__ void abspyr_krnl(Tout *g_odata, Tin *g_idata, unsigned int ns,
    unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  int x, y;

  x = i % ns;
  y = i / ns;
  x = (x + ns / 2) % ns;
  y = (y + ns / 2) % ns;

  unsigned int i2 = x + y * ns;

  if (i < ns * ns) {
    for (int cpt = 0; cpt < nim; cpt++) {
      g_odata[i + cpt * ns * ns] = sqrtf(
          g_idata[i2 + cpt * ns * ns].x * g_idata[i2 + cpt * ns * ns].x
              + g_idata[i2 + cpt * ns * ns].y * g_idata[i2 + cpt * ns * ns].y);
      g_odata[i2 + cpt * ns * ns] = sqrtf(
          g_idata[i + cpt * ns * ns].x * g_idata[i + cpt * ns * ns].x
              + g_idata[i + cpt * ns * ns].y * g_idata[i + cpt * ns * ns].y);
    }
  }
}

template<class Tout, class Tin>
void pyr_abs(Tout *d_odata, Tin *d_idata, int ns, int nim, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, ns * ns / 2, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  abspyr_krnl<Tout, Tin> <<<grid, threads>>>(d_odata, d_idata, ns, nim);

  carmaCheckMsg("abspyr_kernel<<<>>> execution failed\n");
}
template void
pyr_abs<float, cuFloatComplex>(float *d_odata, cuFloatComplex *d_idata, int ns,
    int nim, carma_device *device);
template void
pyr_abs<double, cuDoubleComplex>(double *d_odata, cuDoubleComplex *d_idata,
    int ns, int nim, carma_device *device);

template<class Tin, class Tout>
__global__ void abs2pyr_krnl(Tout *g_odata, Tin *g_idata, Tout fact,
    unsigned int ns, unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  int x, y;

  x = i % ns;
  y = i / ns;
  x = (x + ns / 2) % ns;
  y = (y + ns / 2) % ns;

  unsigned int i2 = x + y * ns;
  //if (i2 < ns/2) i2 = i2 + ns/2 + ns/2*ns;

  Tin tmp1, tmp2;
  if (i < ns * ns / 2) {
    for (int cpt = 0; cpt < nim; cpt++) {
      tmp1 = g_idata[i + cpt * ns * ns];
      tmp2 = g_idata[i2 + cpt * ns * ns];

      g_odata[i + cpt * ns * ns] += (tmp2.x * tmp2.x + tmp2.y * tmp2.y) * fact;
      g_odata[i2 + cpt * ns * ns] += (tmp1.x * tmp1.x + tmp1.y * tmp1.y) * fact;

    }
  }
}

template<class Tin, class Tout>
void pyr_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int ns, int nim,
    carma_device *device) {
  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, ns * ns / 2, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  abs2pyr_krnl<Tin, Tout> <<<grid, threads>>>(d_odata, d_idata, fact, ns,
      nim);

  carmaCheckMsg("abs2pyr_kernel<<<>>> execution failed\n");
}

template void
pyr_abs2<cuFloatComplex, float>(float *d_odata, cuFloatComplex *d_idata,
    float fact, int ns, int nim, carma_device *device);
template void
pyr_abs2<cuDoubleComplex, double>(double *d_odata, cuDoubleComplex *d_idata,
    double fact, int ns, int nim, carma_device *device);

////////////////////////////////////////////////////////
// ADDING PYR_ABS2 MODIFIED FOR ROOF-PRISM: ROOF_ABS2 //
////////////////////////////////////////////////////////

template<class Tin, class Tout>
__global__ void abs2roof_krnl(Tout *g_odata, Tin *g_idata, Tout fact,
    unsigned int ns, unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  int x, y;

  // The following four lines shifts the index to induce a roll
  // along both axes.

  x = i % ns;
  y = i / ns;
  x = (x + ns / 2) % ns;
  y = (y + ns / 2) % ns;

  unsigned int i2 = x + y * ns;
  //if (i2 < ns/2) i2 = i2 + ns/2 + ns/2*ns;

  Tin tmp1, tmp2;
  if (i < ns * ns / 2) {
    for (int cpt = 0; cpt < nim; cpt++) {
      tmp1 = g_idata[i + cpt * ns * ns];
      tmp2 = g_idata[i2 + cpt * ns * ns];

      g_odata[i + cpt * ns * ns] += (tmp2.x * tmp2.x + tmp2.y * tmp2.y) * fact;
      g_odata[i2 + cpt * ns * ns] += (tmp1.x * tmp1.x + tmp1.y * tmp1.y) * fact;
      
    }
  }
}

template<class Tin, class Tout>
void roof_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int ns, int nim,
    carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, ns * ns / 2, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  abs2roof_krnl<Tin, Tout> <<<grid, threads>>>(d_odata, d_idata, fact, ns,
      nim);

  carmaCheckMsg("abs2pyr_kernel<<<>>> execution failed\n");
}

template void
roof_abs2<cuFloatComplex, float>(float *d_odata, cuFloatComplex *d_idata,
    float fact, int ns, int nim, carma_device *device);
template void
roof_abs2<cuDoubleComplex, double>(double *d_odata, cuDoubleComplex *d_idata,
    double fact, int ns, int nim, carma_device *device);

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template<class Tout, class Tin>
__global__ void submask_krnl(Tout *g_odata, Tin *g_mask, unsigned int n) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < n * n) {
    g_odata[i].x = g_odata[i].x * g_mask[i];
    g_odata[i].y = g_odata[i].y * g_mask[i];
  }
}

template<class Tout, class Tin>
void pyr_submask(Tout *d_odata, Tin *d_mask, int n, carma_device *device) {
  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, n * n, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  submask_krnl<Tout, Tin> <<<grid, threads>>>(d_odata, d_mask, n);

  carmaCheckMsg("submask_kernel<<<>>> execution failed\n");
}
template void
pyr_submask<cuFloatComplex, float>(cuFloatComplex *d_odata, float*d_mask, int n,
    carma_device *device);
template void
pyr_submask<cuDoubleComplex, double>(cuDoubleComplex *d_odata, double *d_idata,
    int n, carma_device *device);



template<class T>
__global__ void submaskpyr_krnl(T *g_odata, T *g_mask, unsigned int n) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    T tmp;
    tmp.x = g_odata[i].x * g_mask[i].x - g_odata[i].y* g_mask[i].y;
    tmp.y = g_odata[i].y * g_mask[i].x + g_odata[i].x *g_mask[i].y;

    g_odata[i].x = tmp.x;
    g_odata[i].y = tmp.y;

    i  += blockDim.x * gridDim.x;
  }
}

template<class T>
void pyr_submaskpyr(T *d_odata, T *d_mask, int n, carma_device *device) {
  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, n * n, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  submaskpyr_krnl<T> <<<grid, threads>>>(d_odata, d_mask, n);

  carmaCheckMsg("submask_kernel<<<>>> execution failed\n");
}
template void
pyr_submaskpyr<cuFloatComplex>(cuFloatComplex *d_odata, cuFloatComplex *d_mask, int n,
    carma_device *device);
template void
pyr_submaskpyr<cuDoubleComplex>(cuDoubleComplex *d_odata,cuDoubleComplex  *d_idata,
    int n, carma_device *device);


template<class Tout, class Tin>
__global__ void submask3d_krnl(Tout *g_odata, Tin *g_mask, unsigned int n,
    unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < n * n) {
    for (int cpt = 0; cpt < nim; cpt++) {
      g_odata[i + cpt * n * n].x = g_odata[i + cpt * n * n].x * g_mask[i];
      g_odata[i + cpt * n * n].y = g_odata[i + cpt * n * n].y * g_mask[i];
    }
  }
}

template<class Tout, class Tin>
void pyr_submask3d(Tout *d_odata, Tin *d_mask, int n, int nim, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, n*n, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  submask3d_krnl<Tout, Tin> <<<grid, threads>>>(d_odata, d_mask, n, nim);

  carmaCheckMsg("submask3d_kernel<<<>>> execution failed\n");
}
template void
pyr_submask3d<cuFloatComplex, float>(cuFloatComplex *d_odata, float*d_mask,
    int n, int nim, carma_device *device);
template void
pyr_submask3d<cuDoubleComplex, double>(cuDoubleComplex *d_odata,
    double *d_idata, int n, int nim, carma_device *device);

template<class T>
__global__ void subsum_krnl(T *g_odata, T *g_idata, int *subindx, int *subindy,
    unsigned int ns, unsigned int nvalid, unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int i2 = subindx[i] + subindy[i] * ns;
    g_odata[i] = 0;
    for (int cpt = 0; cpt < nim; cpt++) {
      g_odata[i] += g_idata[i2 + cpt * ns * ns];
    }
  }
}

template<class T>
void pyr_subsum(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns,
    int nvalid, int nim, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, nvalid, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  subsum_krnl<T> <<<grid, threads>>>(d_odata, d_idata, subindx, subindy, ns,
      nvalid, nim);

  carmaCheckMsg("subsum_kernel<<<>>> execution failed\n");
}

template void
pyr_subsum<float>(float *d_odata, float *d_idata, int *subindx, int *subindy,
    int ns, int nvalid, int nim, carma_device *device);
template void
pyr_subsum<double>(double *d_odata, double *d_idata, int *subindx, int *subindy,
    int ns, int nvalid, int nim, carma_device *device);



template<class T>
__global__ void subsum_krnl(T *g_odata, T *g_idata, int *subindx, int *subindy,
    unsigned int ns, unsigned int nvalid) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int iq1 = subindx[i] + subindy[i] * ns;
    int iq2 = subindx[i+nvalid] + subindy[i+nvalid] * ns;
    int iq3 = subindx[i+2*nvalid] + subindy[i+2*nvalid] * ns;
    int iq4 = subindx[i+3*nvalid] + subindy[i+3*nvalid] * ns;
    g_odata[i] = g_idata[iq1]+g_idata[iq2]+g_idata[iq3]+g_idata[iq4];
  }
}

template<class T>
void pyr_subsum(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns,
    int nvalid, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, nvalid, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  subsum_krnl<T> <<<grid, threads>>>(d_odata, d_idata, subindx, subindy, ns,nvalid);

  carmaCheckMsg("subsum_kernel<<<>>> execution failed\n");
}

template void
pyr_subsum<float>(float *d_odata, float *d_idata, int *subindx, int *subindy,
    int ns, int nvalid, carma_device *device);
template void
pyr_subsum<double>(double *d_odata, double *d_idata, int *subindx, int *subindy,
    int ns, int nvalid, carma_device *device);

////////////////////////////////////////////////////////////
// ADDING PYR_SUBSUM MODIFIED FOR HR pyramid              //
////////////////////////////////////////////////////////////

template<class T>
__global__ void subsum2_krnl(T *g_odata, T *g_idata, int *subindx, int *subindy,
    unsigned int ns, unsigned int nvalid) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int iq1 = subindx[i] + subindy[i] * ns;
    int iq2 = subindx[i+nvalid] + subindy[i+nvalid] * ns;
    int iq3 = subindx[i+2*nvalid] + subindy[i+2*nvalid] * ns;
    int iq4 = subindx[i+3*nvalid] + subindy[i+3*nvalid] * ns;
    g_odata[i] = g_idata[iq1] + g_idata[iq2] + g_idata[iq3] + g_idata[iq4];
  }
}

template<class T>
void pyr_subsum2(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns,
    int nvalid, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, nvalid, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  subsum2_krnl<T> <<<grid, threads>>>(d_odata, d_idata, subindx, subindy, ns,
      nvalid);

  carmaCheckMsg("subsum_kernel<<<>>> execution failed\n");
}

template void
pyr_subsum2<float>(float *d_odata, float *d_idata, int *subindx, int *subindy,
    int ns, int nvalid, carma_device *device);
template void
pyr_subsum2<double>(double *d_odata, double *d_idata, int *subindx, int *subindy,
    int ns, int nvalid, carma_device *device);

////////////////////////////////////////////////////////////
// ADDING PYR_SUBSUM MODIFIED FOR ROOF-PRISM: ROOF_SUBSUM //
////////////////////////////////////////////////////////////

template<class T>
__global__ void roof_subsum_krnl(T *g_odata, T *g_idata, int *subindx, int *subindy,
    unsigned int ns, unsigned int nvalid, unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nvalid) {
    int i2 = subindx[i] + subindy[i] * ns;
    g_odata[i] = 0;
    for (int cpt = 0; cpt < nim; cpt++) {
      g_odata[i] += g_idata[i2 + cpt * ns * ns];
    }
  }
}

template<class T>
void roof_subsum(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns,
    int nvalid, int nim, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, nvalid, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  roof_subsum_krnl<T> <<<grid, threads>>>(d_odata, d_idata, subindx, subindy, ns,
      nvalid, nim);

  carmaCheckMsg("subsum_kernel<<<>>> execution failed\n");
}

template void
roof_subsum<float>(float *d_odata, float *d_idata, int *subindx, int *subindy,
    int ns, int nvalid, int nim, carma_device *device);
template void
roof_subsum<double>(double *d_odata, double *d_idata, int *subindx, int *subindy,
    int ns, int nvalid, int nim, carma_device *device);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template<class T>
__global__ void pyrfact_krnl(T *g_data, T fact, unsigned int n,
    unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    for (int cpt = 0; cpt < nim; cpt++) {
      g_data[i + cpt * n * n] = g_data[i + cpt * n * n] * fact;
    }
    i += blockDim.x * gridDim.x;
  }
}

__global__ void pyrfact_krnl(cuFloatComplex *g_data, float fact, unsigned int n,
    unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    for (int cpt = 0; cpt < nim; cpt++) {
      g_data[i + cpt * n * n].x = g_data[i + cpt * n * n].x * fact;
      g_data[i + cpt * n * n].y = g_data[i + cpt * n * n].y * fact;
    }
    i += blockDim.x * gridDim.x;
  }
}

__global__ void pyrfact_krnl(float *g_data, float fact1, float *fact2,
    unsigned int n, unsigned int nim) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  while (i < n * n) {
    for (int cpt = 0; cpt < nim; cpt++) {
      g_data[i + cpt * n * n] = g_data[i + cpt * n * n] * fact1 / fact2[0];
    }
    i += blockDim.x * gridDim.x;
  }
}

template<class T>
void pyr_fact(T *d_data, T fact, int n, int nim, carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, n*n, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  pyrfact_krnl<T> <<<grid, threads>>>(d_data, fact, n, nim);

  carmaCheckMsg("pyrfact_kernel<<<>>> execution failed\n");
}

template void
pyr_fact<float>(float *d_data, float fact, int ns, int nim, carma_device *device);
template void
pyr_fact<double>(double *d_odata, double fact, int ns, int nim, carma_device *device);

void pyr_fact(cuFloatComplex *d_data, float fact, int n, int nim, carma_device *device) {
  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, n*n, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  pyrfact_krnl<<<grid, threads>>>(d_data, fact, n, nim);

  carmaCheckMsg("pyrfact_kernel<<<>>> execution failed\n");
}

void pyr_fact(float *d_data, float fact1, float *fact2, int n, int nim,
    carma_device *device) {

  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, n*n, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  pyrfact_krnl<<<grid, threads>>>(d_data, fact1, fact2, n, nim);

  carmaCheckMsg("pyrfact_kernel<<<>>> execution failed\n");
}

/*
 __global__ void fillcamplipup_krnl(cuFloatComplex *amplipup, float *phase,float *offset, float *mask, int *indx, int Nfft, 
 int Npup, int npup, int N)
 {
 int nim,idim,idimx,idimy,idx;
 int tid   = threadIdx.x + blockIdx.x * blockDim.x;

 while (tid < N) {
 nim   = tid / Npup;
 idim  = tid - nim * Npup;
 idimx = idim % npup;
 idimy = idim / npup;
 idx   = idimx + idimy * Nfft + nim * Nfft * Nfft;
 amplipup[idx].x = (cosf(phase[indx[tid]]-offset[idim]))*mask[indx[tid]];
 amplipup[idx].y = (sinf(phase[indx[tid]]-offset[idim]))*mask[indx[tid]];
 tid  += blockDim.x * gridDim.x;
 }
 }

 int fillcamplipup(cuFloatComplex *amplipup, float *phase, float *offset, float *mask, int *indx, int Nfft, int Npup, int Nsub, 
 int npup, carma_device *device)
 // here amplipup is a cube of data of size nfft x nfft x nsubap
 // phase is an array of size pupdiam x pupdiam
 // offset is an array of size pdiam x pdiam
 // mask is an array of size pupdiam x pupdiam
 // indx is an array of size pdiam x pdiam x nsubap
 // number of thread required : pdiam x pdiam x nsubap
 // Npup = pdiam x pdiam
 {

 struct cudaDeviceProp deviceProperties;
 cudaGetDeviceProperties(&deviceProperties, device);
 
  int nBlocks,nThreads;
  getNumBlocksAndThreads(device, Npup * Nsub, nBlocks, nThreads);
 dim3 grid(nBlocks), threads(nThreads);

 fillcamplipup_krnl<<<grid, threads>>>(amplipup,phase,offset,mask,indx,Nfft,Npup,npup,N);
 carmaCheckMsg("fillcamplipup_kernel<<<>>> execution failed\n");

 return EXIT_SUCCESS;
 }




 */
