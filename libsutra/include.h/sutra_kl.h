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

//! \file      sutra_kl.h
//! \ingroup   libsutra
//! \class     sutra_kl
//! \brief     this class provides the kl features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_KL_H_
#define _SUTRA_KL_H_

#include <carma.h>
#include <carma_host_obj.h>
#include <carma_obj.h>

class sutra_kl {
 public:
  int device;  // # device
  long dim;    // dim of final array
  long nr;     // # radial points
  long np;     // # of elements
  long nkl;    // # of functions in the basis
  long nord;   // # number of radial orders

  carma_obj<float> *d_rabas;   // the radial array of the basis
  carma_obj<float> *d_azbas;   // the azimuthal array of the basis
  carma_host_obj<int> *h_ord;  // the radial orders of the basis
  carma_obj<int> *d_ord;       // the radial orders of the basis
  carma_obj<float> *d_cr;      //
  carma_obj<float> *d_cp;      //

  // Florian features
  carma_obj<float> *d_covmat;
  carma_obj<float> *d_filter;
  carma_obj<float> *d_bas;
  carma_obj<float> *d_evals;

  carma_context *current_context;  // the context in which it has been created
 public:
  sutra_kl(carma_context *context, long dim, long nr, long np, long nkl,
           long nord, int device);
  ~sutra_kl();

  int do_compute(float alpha, float ampli, float *odata, int nkl, int size,
                 int xoff, int yoff);
  int do_compute(float ampli, float *odata, int nkl, int size, int xoff,
                 int yoff);
  int do_compute(float *odata, int nkl, int size, int xoff, int yoff);
  int do_combi(float *com, float *odata, int size, int xoff, int yoff);

  // Florian features
  int get_flokl();
};
int getkl(float alpha, float ampli, float *d_odata, float *rabas, float *azbas,
          float *cr, float *cp, int nr, int np, int nx, int Nx, int xoff,
          int yoff);
int getkl(float ampli, float *d_odata, float *rabas, float *azbas, float *cr,
          float *cp, int nr, int np, int nx, int Nx, int xoff, int yoff);
int getkl(float *d_odata, float *rabas, float *azbas, float *cr, float *cp,
          int nr, int np, int nx, int Nx, int xoff, int yoff);
int combikl(float *com, int nkl, float *d_odata, float *rabas, int *h_ord,
            float *azbas, float *cr, float *cp, int nr, int np, int nx, int Nx,
            int xoff, int yoff);
int cget_flokl(long nkl, long dim, float *covmat, float *filter, float *bas);
// template <class T> void comp_kl(int threads, int blocks, T *d_idata, T
// *d_odata, int N);

#endif  // _SUTRA_KL_H_
