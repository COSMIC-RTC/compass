#include <sutra_roket.h>

__global__ void separate_modes_krnl(float *modes, float *filtmodes, int nmodes, int nfilt){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int bornemin = nmodes - nfilt - 2;
  int bornemax = nmodes - 2;
  while(tid < nmodes){
    if(tid >= bornemin && tid < bornemax){
      filtmodes[tid] = modes[tid];
      modes[tid] = 0.0f;
    }
    else filtmodes[tid] = 0.0f;

    tid += blockDim.x * gridDim.x;
  }
}

int separate_modes(float *modes, float *filtmodes, int nmodes, int nfilt, carma_device *device){
  int nthreads = 0, nblocks = 0;
	getNumBlocksAndThreads(device, nmodes, nblocks, nthreads);
	dim3 grid(nblocks), threads(nthreads);

  separate_modes_krnl<<< grid, threads >>>(modes,filtmodes,nmodes,nfilt);
  carmaCheckMsg("separate_modes_krnl<<<>>> execution failed\n");

	return EXIT_SUCCESS;
}
