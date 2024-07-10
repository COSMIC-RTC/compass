// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_transpose.cu
//! \ingroup   libcarma
//! \brief     this file provides transpose CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <carma_obj.hpp>

/*
 *  _____                                           _  __                    _
 * |_   _| __ __ _ _ __  ___ _ __   ___  ___  ___  | |/ /___ _ __ _ __   ___| |___
 *   | || '__/ _` | '_ \/ __| '_ \ / _ \/ __|/ _ \ | ' // _ \ '__| '_ \ / _ \ / __|
 *   | || | | (_| | | | \__ \ |_) | (_) \__ \  __/ | . \  __/ |  | | | |  __/ \__ \
 *   |_||_|  \__,_|_| |_|___/ .__/ \___/|___/\___| |_|\_\___|_|  |_| |_|\___|_|___/
 *                          |_|
 */

// Transpose that effectively reorders execution of thread blocks along
// diagonals of the matrix (also coalesced and has no bank conflicts)
//
// Here blockIdx.x is interpreted as the distance along a diagonal and
// blockIdx.y as corresponding to different diagonals
//
// blockIdx_x and blockIdx_y expressions map the diagonal coordinates to the
// more commonly used cartesian coordinates so that the only changes to the code
// from the coalesced version are the calculation of the blockIdx_x and
// blockIdx_y and replacement of blockIdx.x and bloclIdx.y with the subscripted
// versions in the remaining code
#define TILE_DIM 32
#define BLOCK_ROWS 16

#define FLOOR(a, b) (a - (a % b))

template <class T>
__global__ void transposeDiagonal(T *odata, T *idata, int64_t width, int64_t height,
                                  int32_t nreps) {
  __shared__ T tile[TILE_DIM][TILE_DIM + 1];

  int32_t blockIdx_x, blockIdx_y;

  // do diagonal reordering
  if (width == height) {
    blockIdx_y = blockIdx.x;
    blockIdx_x = (blockIdx.x + blockIdx.y) % gridDim.x;
  } else {
    int32_t bid = blockIdx.x + gridDim.x * blockIdx.y;
    blockIdx_y = bid % gridDim.y;
    blockIdx_x = ((bid / gridDim.y) + blockIdx_y) % gridDim.x;
  }

  // from here on the code is same as previous kernel except blockIdx_x replaces
  // blockIdx.x and similarly for y

  int32_t xIndex = blockIdx_x * TILE_DIM + threadIdx.x;
  int32_t yIndex = blockIdx_y * TILE_DIM + threadIdx.y;
  int32_t index_in = xIndex + (yIndex)*width;

  xIndex = blockIdx_y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx_x * TILE_DIM + threadIdx.y;
  int32_t index_out = xIndex + (yIndex)*height;

  for (int32_t r = 0; r < nreps; r++) {
    for (int32_t i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
      tile[threadIdx.y + i][threadIdx.x] = idata[index_in + i * width];
    }

    __syncthreads();

    for (int32_t i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
      odata[index_out + i * height] = tile[threadIdx.x][threadIdx.y + i];
    }
  }
}

/*
 template<class T> int32_t get_tdim()
 {
 return 32;
 }
 template<> int32_t get_tdim<float>()
 {
 struct cudaDeviceProp deviceProperties;
 cudaGetDeviceProperties(&deviceProperties, 0);
 int32_t shmemSize = deviceProperties.sharedMemPerBlock;
 return shmemSize/4;
 }
 template<> int32_t get_tdim<double>()
 {
 struct cudaDeviceProp deviceProperties;
 cudaGetDeviceProperties(&deviceProperties, 0);
 int32_t shmemSize = deviceProperties.sharedMemPerBlock;
 return shmemSize/8;
 }
 */

template <class T>
int32_t transposeCU(T *d_idata, T *d_odata, int64_t N1, int64_t N2) {
  /*
   int32_t totTile = get_tdim<T>();
   fprintf(stderr,"tot Tiles : %d\n",totTile);
   fprintf(stderr,"Tile dim : %f\n",sqrt(totTile));
   //TILE_DIM = totTile / N1 / N2;

   struct cudaDeviceProp deviceProperties;
   cudaGetDeviceProperties(&deviceProperties, 0);

   int32_t brows = deviceProperties.maxThreadsPerBlock;

   fprintf(stderr,"block rows : %f\n", brows/sqrt(totTile) );

   //BLOCK_ROWS = TILE_DIM / deviceProperties.maxThreadsPerBlock;
   */

  dim3 grid(N1 / TILE_DIM, N2 / TILE_DIM);
  dim3 threads(TILE_DIM, BLOCK_ROWS);

  transposeDiagonal<<<grid, threads>>>(d_odata, d_idata, N1, N2, 1);

  return EXIT_SUCCESS;
}

template int32_t transposeCU<int32_t>(int32_t *d_idata, int32_t *d_odata, int64_t N1, int64_t N2);

template int32_t transposeCU<uint32_t>(uint32_t *d_idata,
                                       uint32_t *d_odata, int64_t N1, int64_t N2);

template int32_t transposeCU<uint16_t>(uint16_t *d_idata, uint16_t *d_odata,
                                   int64_t N1, int64_t N2);
template int32_t transposeCU<float>(float *d_idata, float *d_odata, int64_t N1,
                                int64_t N2);

template int32_t transposeCU<double>(double *d_idata, double *d_odata, int64_t N1,
                                 int64_t N2);

template int32_t transposeCU<cuFloatComplex>(cuFloatComplex *d_idata,
                                         cuFloatComplex *d_odata, int64_t N1,
                                         int64_t N2);

template int32_t transposeCU<cuDoubleComplex>(cuDoubleComplex *d_idata,
                                          cuDoubleComplex *d_odata, int64_t N1,
                                          int64_t N2);
// template int32_t transposeCU<tuple_t<float>>(tuple_t<float> *d_idata,
//                                          tuple_t<float> *d_odata, int64_t N1,
//                                          int64_t N2);
#ifdef CAN_DO_HALF
template int32_t transposeCU<half>(half *d_idata, half *d_odata, int64_t N1, int64_t N2);
#endif
