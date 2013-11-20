#include <sutra_acquisim.h>
#include <sutra_ao_utils.h>

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T> struct SharedMemory
{
    __device__ inline operator       T*()
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }

    __device__ inline operator const T*() const
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }
};

// specialize for double to avoid unaligned memory
// access compile errors
template<> struct SharedMemory<double>
{
    __device__ inline operator       double*()
    {
        extern __shared__ double __smem_d[];
        return (double*)__smem_d;
    }

    __device__ inline operator const double*() const
    {
        extern __shared__ double __smem_d[];
        return (double*)__smem_d;
    }
};

__global__ void bcube_krnl(float *bimage, float *bcube, int npix, int npix2, int nsub, int *ivalid, int *jvalid, int N)
{
  /*
    indx is an array nrebin^2 * npix^2
    it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
    Npix = npix x npix
   */
  int tid     = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    int nim = tid / npix2;
    int tidim = tid - nim * npix2;
    int xim = tidim % npix;
    int yim = tidim / npix;

    int idbin = xim + yim * nsub + ivalid[nim] * npix + jvalid[nim] * npix * nsub;
    bcube[tid] = bimage[idbin];
    tid  += blockDim.x * gridDim.x;
  }
}

int fillbincube(float *bimage, float *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid, int device)
{
  int Npix = npix * npix;
  int N = Npix * nsub;
  int nthreads = 0,nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  bcube_krnl<<<grid, threads>>>(bimage,bcube,npix,Npix,Nsub,ivalid,jvalid,N);

  cutilCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void bcube_krnl_async(float *bimage, float *bcube, int npix, int npix2, int nsub, int *ivalid, int *jvalid, int N, int idstart)
{
  /*
    indx is an array nrebin^2 * npix^2
    it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
    Npix = npix x npix
   */
  int tid     = threadIdx.x + blockIdx.x * blockDim.x;
  tid += idstart;

  while (tid < N) {
    int nim = tid / npix2;
    int tidim = tid - nim * npix2;
    int xim = tidim % npix;
    int yim = tidim / npix;

    int idbin = xim + yim * nsub + ivalid[nim] * npix + jvalid[nim] * npix * nsub;
    bcube[tid] = bimage[idbin];
    tid  += blockDim.x * gridDim.x;
  }
}

int fillbincube_async(carma_host_obj<float> *image_telemetry, float *bimage, float *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid, int nim, int device)
{
  float *hdata = image_telemetry->getData();
  int nstreams = image_telemetry->get_nbStreams();

  int Npix = npix * npix;
  int N = Npix * nsub;
  int nthreads = 0,nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  // here nstreams should be : final image size / npix
  dim3 threads(nthreads);
  dim3 grid(N/(nstreams*threads.x));

  // asynchronously launch nstreams kernels, each operating on its own portion of data
  for(int i = 0; i < nstreams; i++){
	  cudaMemcpyAsync(&(bimage[i*nim/nstreams]), &(hdata[i*nim/nstreams]), sizeof(float) * nim / nstreams,
			  cudaMemcpyHostToDevice, image_telemetry->get_cudaStream_t(i));
    bcube_krnl_async<<<grid, threads, 0, image_telemetry->get_cudaStream_t(i)>>>(bimage,bcube,npix,Npix,Nsub,ivalid,jvalid,N,i*N/nstreams);
    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x will only
    //   commence executing when all previous CUDA calls in stream x have completed
  }
  //cudaStreamSynchronize(image_telemetry->get_cudaStream_t(nstreams-1));
  cutilCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

int fillbincube_async(carma_streams *streams, carma_obj<float> *bimage, carma_obj<float> *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid, int device)
{
  float *g_image = bimage->getData();
  float *g_cube = bcube->getData();
  int nstreams = streams->get_nbStreams();

  int Npix = npix * npix;
  int N = Npix * nsub;
  int nthreads = 0,nblocks = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  // here nstreams should be : final image size / npix
  dim3 threads(nthreads);
  dim3 grid(N/(nstreams*threads.x));

  // asynchronously launch nstreams kernels, each operating on its own portion of data
  for(int i = 0; i < nstreams; i++){

    bcube_krnl_async<<<grid, threads, 0, streams->get_stream(i)>>>(g_image,g_cube,npix,Npix,Nsub,ivalid,jvalid,N,i*N/nstreams);
    // asynchronously launch nstreams memcopies.  Note that memcopy in stream x will only
    //   commence executing when all previous CUDA calls in stream x have completed
  }

  cutilCheckMsg("binimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}


