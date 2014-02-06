#if 0
#include <carma_svd.h>

__global__ void kernel_setidd(double *d,int N)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) {
    d[idx+N*idx] = 1.0;
  }
}

int carma_setidd(double *d,int n)
{
  int blockSize = 8; 
  int nBlocks = n / blockSize + (n % blockSize == 0?0:1);
  
  kernel_setidd <<< nBlocks, blockSize >>> ((double *)d,n);
   
  return EXIT_SUCCESS;
}

__global__ void kernel_setidf(float *d,int N)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) {
    d[idx+N*idx] = 1.0;
  }
}

int carma_setidf(float *d,int n)
{
  int blockSize = 8; 
  int nBlocks = n / blockSize + (n % blockSize == 0?0:1);
  
  kernel_setidf <<< nBlocks, blockSize >>> ((float *)d,n);
   
  return EXIT_SUCCESS;
}
#endif
