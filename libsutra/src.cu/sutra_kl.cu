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

//! \file      sutra_kl.cu
//! \ingroup   libsutra
//! \class     sutra_kl
//! \brief     this class provides the kl features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_kl.h>

__device__ float kl_sfi(float *rabas, float *azbas, int npix, int nrow) {
  return rabas[npix] * azbas[nrow];
}

__device__ void kl_interp(float alpha, float ampli, float *odata, float *rabas,
                          float *azbas, float xbi, float ybi, int nr, int np,
                          int tido) {
  if ((xbi >= 0) && (xbi <= nr - 1) && (ybi >= 0) && (ybi <= np - 1)) {
    int i0, i1, j0, j1;
    long ibi = (long)xbi;
    long jbi = (long)ybi;

    i0 = ibi;
    if (ibi < nr - 1)
      i1 = ibi + 1;
    else
      i1 = nr - 1;

    j0 = jbi;
    if (jbi < np - 1)
      j1 = jbi + 1;
    else
      j1 = np - 1;

    float wi, wj, w00, w01, w10, w11;

    wi = 1.0f - (xbi - ibi);
    wj = 1.0f - (ybi - jbi);

    w00 = wi * wj;
    w10 = (1 - wi) * wj;
    w01 = wi * (1 - wj);
    w11 = (1 - wi) * (1 - wj);

    odata[tido] =
        alpha * odata[tido] + ampli * (w00 * kl_sfi(rabas, azbas, i0, j0) +
                                       w10 * kl_sfi(rabas, azbas, i1, j0) +
                                       w01 * kl_sfi(rabas, azbas, i0, j1) +
                                       w11 * kl_sfi(rabas, azbas, i1, j1));
  } else
    odata[tido] = 0.0f;
}

__device__ void kl_interp(float alpha, float *ampli, int nkl, float *odata,
                          float *rabas, float *azbas, float xbi, float ybi,
                          int nr, int np, int tido) {
  if ((xbi >= 0) && (xbi <= nr - 1) && (ybi >= 0) && (ybi <= np - 1)) {
    int i0, i1, j0, j1;
    long ibi = (long)xbi;
    long jbi = (long)ybi;

    i0 = ibi;
    if (ibi < nr - 1)
      i1 = ibi + 1;
    else
      i1 = nr - 1;

    j0 = jbi;
    if (jbi < np - 1)
      j1 = jbi + 1;
    else
      j1 = np - 1;

    float wi, wj, w00, w01, w10, w11;

    wi = 1.0f - (xbi - ibi);
    wj = 1.0f - (ybi - jbi);

    w00 = wi * wj;
    w10 = (1 - wi) * wj;
    w01 = wi * (1 - wj);
    w11 = (1 - wi) * (1 - wj);

    odata[tido] =
        alpha * odata[tido] + ampli[nkl] * (w00 * kl_sfi(rabas, azbas, i0, j0) +
                                            w10 * kl_sfi(rabas, azbas, i1, j0) +
                                            w01 * kl_sfi(rabas, azbas, i0, j1) +
                                            w11 * kl_sfi(rabas, azbas, i1, j1));
  } else
    odata[tido] = 0.0f;
}

__global__ void getkl_krnl(float alpha, float ampli, float *odata, float *rabas,
                           float *azbas, float *cr, float *cp, int nr, int np,
                           int nx, int Nx, int xoff, int yoff) {
  int xid = threadIdx.x + blockIdx.x * blockDim.x;
  int yid = threadIdx.y + blockIdx.y * blockDim.y;

  if ((xid < nx) && (yid < nx)) {
    int xref = xid + xoff;
    int yref = yid + yoff;
    int tido = xref + yref * Nx;

    int tid = xid + yid * nx;
    float xbi = cr[tid];
    float ybi = cp[tid];

    kl_interp(alpha, ampli, odata, rabas, azbas, xbi, ybi, nr, np, tido);
  }
}

__global__ void combikl_krnl(float *com, int nkl, float *odata, float *rabas,
                             int *d_ord, float *azbas, float *cr, float *cp,
                             int nr, int np, int nx, int Nx, int xoff,
                             int yoff) {
  int xid = threadIdx.x + blockIdx.x * blockDim.x;
  int yid = threadIdx.y + blockIdx.y * blockDim.y;

  if ((xid < nx) && (yid < nx)) {
    int xref = xid + xoff;
    int yref = yid + yoff;
    int tido = xref + yref * Nx;

    int tid = xid + yid * nx;
    float xbi = cr[tid];
    float ybi = cp[tid];

    odata[tido] = 0.0f;
    int tmp;
    float *rabas_cc;
    float *azbas_cc;
    // int cc=10;
    for (int cc = 0; cc < nkl; cc++) {
      tmp = d_ord[cc] - 1;
      rabas_cc = &(rabas[cc * nr]);
      azbas_cc = &(azbas[tmp * np]);
      kl_interp(1.0f, com, cc, odata, rabas_cc, azbas_cc, xbi, ybi, nr, np,
                tido);
    }
    __syncthreads();
  }
}

int getkl(float alpha, float ampli, float *d_odata, float *rabas, float *azbas,
          float *cr, float *cp, int nr, int np, int nx, int Nx, int xoff,
          int yoff) {
  int block_size = 8;
  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  int nny = nx + block_size - nx % block_size;
  dim3 blocks(nnx / block_size, nny / block_size),
      threads(block_size, block_size);

  // int smemSize = (block_size +1) * (block_size +1) * sizeof(float);

  getkl_krnl<<<blocks, threads>>>(alpha, ampli, d_odata, rabas, azbas, cr, cp,
                                  nr, np, nx, Nx, xoff, yoff);

  carmaCheckMsg("get_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

int getkl(float ampli, float *d_odata, float *rabas, float *azbas, float *cr,
          float *cp, int nr, int np, int nx, int Nx, int xoff, int yoff) {
  return getkl(0.0f, ampli, d_odata, rabas, azbas, cr, cp, nr, np, nx, Nx, xoff,
               yoff);
}

int getkl(float *d_odata, float *rabas, float *azbas, float *cr, float *cp,
          int nr, int np, int nx, int Nx, int xoff, int yoff) {
  return getkl(0.0f, 1.0f, d_odata, rabas, azbas, cr, cp, nr, np, nx, Nx, xoff,
               yoff);
}

int combikl(float *com, int nkl, float *d_odata, float *rabas, int *d_ord,
            float *azbas, float *cr, float *cp, int nr, int np, int nx, int Nx,
            int xoff, int yoff) {
  int block_size = 8;
  int nnx =
      nx + block_size - nx % block_size;  // find next multiple of BLOCK_SZ
  int nny = nx + block_size - nx % block_size;
  dim3 blocks(nnx / block_size, nny / block_size),
      threads(block_size, block_size);

  // int smemSize = (block_size +1) * (block_size +1) * sizeof(float);

  // for (int cc=0;cc<nkl;cc++)
  combikl_krnl<<<blocks, threads>>>(com, nkl, d_odata, rabas, d_ord, azbas, cr,
                                    cp, nr, np, nx, Nx, xoff, yoff);

  carmaCheckMsg("get_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

// Florian features
__global__ void flokrnl(long dim, float *bas) {
  int tid = blockIdx.x;
  if (tid < dim) bas[tid * dim + tid] = tid;
}

int cget_flokl(long nkl, long dim, float *covmat, float *filter, float *bas) {
  // int i;
  printf("flag CUDA \n");
  // for (i=0;i<dim;i++) bas[i] = i;
  flokrnl<<<dim, 1>>>(dim, bas);
  printf("flag CUDA done \n");
  carmaCheckMsg("get_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}
