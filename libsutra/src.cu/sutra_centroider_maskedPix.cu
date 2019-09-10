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

//! \file      sutra_centroider_maskedPix.cu
//! \ingroup   libsutra
//! \class     sutra_centroider_maskedPix
//! \brief     this class provides the centroider_maskedPix features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#include <sutra_centroider_maskedPix.h>
#include <carma_utils.cuh>

template <class T>
__global__ void get_maskedPix_krnl(T *g_odata, T *ref, float *g_idata,
                                   int *subindx, int *subindy,
                                   float *intensities, int ns, int nslopes) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < nslopes) {
    int i2 = subindx[i] + subindy[i] * ns;
    g_odata[i] = T(g_idata[i2] / intensities[0]) - ref[i];
  }
}

template <class T>
void getMaskedPix(T *d_odata, T *ref, float *d_idata, int *subindx,
                  int *subindy, float *intensities, int ns, int nslopes,
                  carma_device *device) {
  // cout << "hello cu" << endl;

  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nslopes, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  get_maskedPix_krnl<T><<<grid, threads>>>(d_odata, ref, d_idata, subindx,
                                           subindy, intensities, ns, nslopes);

  carmaCheckMsg("get_maskedPix_kernel<<<>>> execution failed\n");
}

template void getMaskedPix<float>(float *d_odata, float *ref, float *d_idata,
                                  int *subindx, int *subindy,
                                  float *intensities, int ns, int nslopes,
                                  carma_device *device);

#ifdef CAN_DO_HALF
template void getMaskedPix<half>(half *d_odata, half *ref, float *d_idata,
                                 int *subindx, int *subindy, float *intensities,
                                 int ns, int nslopes, carma_device *device);
#endif

__global__ void fill_intensities_krnl(float *g_odata, float *g_idata,
                                      int *subindx, int *subindy, int ns,
                                      int nslopes) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < nslopes) {
    int i = subindx[tid] + subindy[tid] * ns;
    g_odata[tid] = g_idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

void fill_intensities(float *d_odata, float *d_idata, int *subindx,
                      int *subindy, int ns, int nslopes, carma_device *device) {
  int nBlocks, nThreads;
  getNumBlocksAndThreads(device, nslopes, nBlocks, nThreads);
  dim3 grid(nBlocks), threads(nThreads);

  fill_intensities_krnl<<<grid, threads>>>(d_odata, d_idata, subindx, subindy,
                                           ns, nslopes);

  carmaCheckMsg("fill_intensities_kernel<<<>>> execution failed\n");
}
