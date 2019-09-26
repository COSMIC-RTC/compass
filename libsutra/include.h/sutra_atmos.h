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

//! \file      sutra_atmos.h
//! \ingroup   libsutra
//! \class     sutra_atmos
//! \brief     this class provides the atmos features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_ATMOS_H_
#define _SUTRA_ATMOS_H_

#include <sutra_tscreen.h>

using std::pair;
using std::vector;

class sutra_atmos {
 public:
  int nscreens;
  vector<sutra_tscreen *> d_screens;
  float r0;
  carma_context *current_context;

 public:
  sutra_atmos(carma_context *context, int nscreens, float global_r0, float *r0,
              long *size, long *size2, float *altitude, float *windspeed,
              float *winddir, float *deltax, float *deltay, int device);
  ~sutra_atmos();

  int init_screen(int idx, float *h_A, float *h_B, unsigned int *h_istencilx,
                  unsigned int *h_istencily, int seed);

  int add_screen(float altitude, long size, long stencilSize, float amplitude,
                 float windspeed, float winddir, float deltax, float deltay,
                 int device);
  int del_screen(const int idx);
  int refresh_screen(int idx);

  int move_atmos();
  int set_global_r0(float r0);
  int set_frac(float *frac);
  int set_seed(int idx, float seed);
};

#endif  // _SUTRA_ATMOS_H_
