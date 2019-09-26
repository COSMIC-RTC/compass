// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser 
//  General Public License as published by the Free Software Foundation, either version 3 of the License, 
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration 
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems. 
//  
//  The final product includes a software package for simulating all the critical subcomponents of AO, 
//  particularly in the context of the ELT and a real-time core based on several control approaches, 
//  with performances consistent with its integration into an instrument. Taking advantage of the specific 
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT. 
//  
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components 
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and 
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      carma_transpose.cu
//! \ingroup   libcarma
//! \brief     this file provides transpose CUDA kernels
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <carma_obj.h>

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
__global__ void transposeDiagonal(T *odata, T *idata, long width, long height,
                                  int nreps) {
  __shared__ T tile[TILE_DIM][TILE_DIM + 1];

  int blockIdx_x, blockIdx_y;

  // do diagonal reordering
  if (width == height) {
    blockIdx_y = blockIdx.x;
    blockIdx_x = (blockIdx.x + blockIdx.y) % gridDim.x;
  } else {
    int bid = blockIdx.x + gridDim.x * blockIdx.y;
    blockIdx_y = bid % gridDim.y;
    blockIdx_x = ((bid / gridDim.y) + blockIdx_y) % gridDim.x;
  }

  // from here on the code is same as previous kernel except blockIdx_x replaces
  // blockIdx.x and similarly for y

  int xIndex = blockIdx_x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx_y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*width;

  xIndex = blockIdx_y * TILE_DIM + threadIdx.x;
  yIndex = blockIdx_x * TILE_DIM + threadIdx.y;
  int index_out = xIndex + (yIndex)*height;

  for (int r = 0; r < nreps; r++) {
    for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
      tile[threadIdx.y + i][threadIdx.x] = idata[index_in + i * width];
    }

    __syncthreads();

    for (int i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
      odata[index_out + i * height] = tile[threadIdx.x][threadIdx.y + i];
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

template <class T>
int transposeCU(T *d_idata, T *d_odata, long N1, long N2) {
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

  dim3 grid(N1 / TILE_DIM, N2 / TILE_DIM);
  dim3 threads(TILE_DIM, BLOCK_ROWS);

  transposeDiagonal<<<grid, threads>>>(d_odata, d_idata, N1, N2, 1);

  return EXIT_SUCCESS;
}

template int transposeCU<int>(int *d_idata, int *d_odata, long N1, long N2);

template int transposeCU<unsigned int>(unsigned int *d_idata,
                                       unsigned int *d_odata, long N1, long N2);

template int transposeCU<uint16_t>(uint16_t *d_idata, uint16_t *d_odata,
                                   long N1, long N2);
template int transposeCU<float>(float *d_idata, float *d_odata, long N1,
                                long N2);

template int transposeCU<double>(double *d_idata, double *d_odata, long N1,
                                 long N2);

template int transposeCU<cuFloatComplex>(cuFloatComplex *d_idata,
                                         cuFloatComplex *d_odata, long N1,
                                         long N2);

template int transposeCU<cuDoubleComplex>(cuDoubleComplex *d_idata,
                                          cuDoubleComplex *d_odata, long N1,
                                          long N2);
// template int transposeCU<tuple_t<float>>(tuple_t<float> *d_idata,
//                                          tuple_t<float> *d_odata, long N1,
//                                          long N2);
#ifdef CAN_DO_HALF
template int transposeCU<half>(half *d_idata, half *d_odata, long N1, long N2);
#endif
