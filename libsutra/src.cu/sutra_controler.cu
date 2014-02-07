#include <sutra_controller.h>


/*
 _  __                    _     
| |/ /___ _ __ _ __   ___| |___ 
| ' // _ \ '__| '_ \ / _ \ / __|
| . \  __/ |  | | | |  __/ \__ \
|_|\_\___|_|  |_| |_|\___|_|___/
                                
 */


__global__ void shift_krnl(float *data, int offset, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    data[tid] = data[tid+offset*N]; 
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void mult_krnl(float *i_data, float *scale, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    i_data[tid] = i_data[tid]*scale[tid]; 
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void mult_int_krnl(float *o_data, float *i_data, float *scale, float gain, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    o_data[tid] = gain *( i_data[tid] * scale[tid]) + o_data[tid]; 
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void mult_int_krnl(float *o_data, float *i_data, float *scale, float gain, int N, int istart)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  tid    += istart;

  if (tid < N) {
    o_data[tid] = gain *( i_data[tid] * scale[tid]) + o_data[tid]; 
  }
}

/*
 _                           _                   
| |    __ _ _   _ _ __   ___| |__   ___ _ __ ___ 
| |   / _` | | | | '_ \ / __| '_ \ / _ \ '__/ __|
| |__| (_| | |_| | | | | (__| | | |  __/ |  \__ \
|_____\__,_|\__,_|_| |_|\___|_| |_|\___|_|  |___/
                                                 
 */


int shift_buf(float *d_data,int offset, int N,int device)
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
    
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (N + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);

  shift_krnl<<<grid, threads>>>(d_data, offset, N);

  cutilCheckMsg("shift_kernel<<<>>> execution failed\n");
   return EXIT_SUCCESS;
}

int mult_vect(float *d_data,float *scale, int N,int device)
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
    
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (N + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);

  mult_krnl<<<grid, threads>>>(d_data, scale, N);

  cutilCheckMsg("mult_kernel<<<>>> execution failed\n");
   return EXIT_SUCCESS;
}

int mult_int(float *o_data,float *i_data,float *scale, float gain, int N,int device, carma_streams *streams)
{

  int nthreads = 0,nblocks = 0;

  int nstreams = streams->get_nbStreams();
  getNumBlocksAndThreads(device, N/nstreams, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  for(int i = 0; i < nstreams; i++) {
    mult_int_krnl<<<grid, threads, 0, streams->get_stream(i)>>>(o_data,i_data, scale, gain, N, i*nblocks*nthreads);
  cutilCheckMsg("multint_kernel<<<>>> execution failed\n");
  }

   return EXIT_SUCCESS;
}

int mult_int(float *o_data,float *i_data,float *scale, float gain, int N,int device)
{

  int nthreads = 0,nblocks = 0;

  getNumBlocksAndThreads(device, N, nblocks, nthreads);

  dim3 grid(nblocks), threads(nthreads);

  mult_int_krnl<<<grid, threads>>>(o_data,i_data, scale, gain, N);
  cutilCheckMsg("multint_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
