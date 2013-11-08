#include <carma_obj.h>

/*
 _____                                           _  __                    _     
|_   _| __ __ _ _ __  ___ _ __   ___  ___  ___  | |/ /___ _ __ _ __   ___| |___ 
  | || '__/ _` | '_ \/ __| '_ \ / _ \/ __|/ _ \ | ' // _ \ '__| '_ \ / _ \ / __|
  | || | | (_| | | | \__ \ |_) | (_) \__ \  __/ | . \  __/ |  | | | |  __/ \__ \
  |_||_|  \__,_|_| |_|___/ .__/ \___/|___/\___| |_|\_\___|_|  |_| |_|\___|_|___/
                         |_|                                                    
*/

// Transpose that effectively reorders execution of thread blocks along diagonals of the 
// matrix (also coalesced and has no bank conflicts)
//
// Here blockIdx.x is interpreted as the distance along a diagonal and blockIdx.y as 
// corresponding to different diagonals
//
// blockIdx_x and blockIdx_y expressions map the diagonal coordinates to the more commonly 
// used cartesian coordinates so that the only changes to the code from the coalesced version 
// are the calculation of the blockIdx_x and blockIdx_y and replacement of blockIdx.x and 
// bloclIdx.y with the subscripted versions in the remaining code

#define TILE_DIM    32
#define BLOCK_ROWS  16

#define FLOOR(a,b) (a-(a%b))

template<class T> __global__ void transposeDiagonal(T *odata, T *idata, long width, long height, int nreps)
{
  __shared__ T tile[TILE_DIM][TILE_DIM+1];

  int blockIdx_x, blockIdx_y;

  // do diagonal reordering
  if (width == height) {
    blockIdx_y = blockIdx.x;
    blockIdx_x = (blockIdx.x+blockIdx.y)%gridDim.x;
  } else {
    int bid = blockIdx.x + gridDim.x*blockIdx.y;
    blockIdx_y = bid%gridDim.y;
    blockIdx_x = ((bid/gridDim.y)+blockIdx_y)%gridDim.x;
  }    

  // from here on the code is same as previous kernel except blockIdx_x replaces blockIdx.x
  // and similarly for y

  int xIndex = blockIdx_x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx_y * TILE_DIM + threadIdx.y;  
  int index_in = xIndex + (yIndex)*width;

  xIndex = blockIdx_y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx_x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;

  for (int r=0; r < nreps; r++) {
    for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
      tile[threadIdx.y+i][threadIdx.x] = idata[index_in+i*width];
    }
  
    __syncthreads();
  
    for (int i=0; i<TILE_DIM; i+=BLOCK_ROWS) {
      odata[index_out+i*height] = tile[threadIdx.x][threadIdx.y+i];
    }
  }
}

/*
template<class T> int get_tdim()
{
  return 32;
}
template<> int get_tdim<float>()
{
  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, 0);
  int shmemSize = deviceProperties.sharedMemPerBlock;
  return shmemSize/4;
}
template<> int get_tdim<double>()
{
  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, 0);
  int shmemSize = deviceProperties.sharedMemPerBlock;
  return shmemSize/8;
}
*/

template<class T> int transposeCU(T *d_idata,T *d_odata,long N1,long N2)
{
  /*
  int totTile = get_tdim<T>();
  fprintf(stderr,"tot Tiles : %d\n",totTile);
  fprintf(stderr,"Tile dim : %f\n",sqrt(totTile));
  //TILE_DIM = totTile / N1 / N2;

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, 0);

  int brows = deviceProperties.maxThreadsPerBlock;

  fprintf(stderr,"block rows : %f\n", brows/sqrt(totTile) );

  //BLOCK_ROWS = TILE_DIM / deviceProperties.maxThreadsPerBlock;
  */

  dim3 grid(N1/TILE_DIM,N2/TILE_DIM);
  dim3 threads(TILE_DIM,BLOCK_ROWS);

  transposeDiagonal<<<grid, threads>>>(d_odata, d_idata, N1 , N2 , 1);

   return EXIT_SUCCESS;
}

template int transposeCU<float>(float *d_idata,float *d_odata,long N1,long N2);

template int transposeCU<double>(double *d_idata,double *d_odata,long N1,long N2);


