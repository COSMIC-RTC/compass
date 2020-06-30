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

//! \file      sutra_telescope.h
//! \ingroup   libsutra
//! \class     SutraTelescope
//! \brief     this class provides the telescope features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_TELESCOPE_H_
#define _SUTRA_TELESCOPE_H_

#include <carma.h>
#include <sutra_phase.h>
#include <map>
#include <string>
#include <vector>

using std::string;

class SutraTelescope {
 public:
  CarmaContext *current_context;
  int device;  // device #

  long pup_size;       // size of pupil
  long num_eleme_pup;  // number of points in the pupil

  CarmaObj<float> *d_pupil;        // the pupil mask
  CarmaObj<float> *d_phase_ab_M1;  // the phase aberration for M1

  long pup_size_m;  // size of pupil

  CarmaObj<float> *d_pupil_m;        // the pupil mask
  CarmaObj<float> *d_phase_ab_M1_m;  // the phase aberration for M1

 public:
  SutraTelescope(CarmaContext *context, long pup_size, long num_eleme_pup,
                  float *pupil, long pup_size_m, float *pupil_m);
  ~SutraTelescope();
  int set_phase_ab_M1(float *phase_ab_M1, int size);
  int set_phase_ab_M1_m(float *phase_ab_M1_m, int size);
};

#endif  // _SUTRA_TELESCOPE_H_
