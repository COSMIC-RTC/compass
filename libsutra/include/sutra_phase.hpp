// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_phase.hpp
//! \ingroup   libsutra
//! \class     SutraPhase
//! \brief     this class provides the phase features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_PHASE_H_
#define _SUTRA_PHASE_H_

// this is the generic class for a phase
// contains a yoga obj for the phase screen itself
// the screen size (assumed square)
// the following is only initialized on demand :
// an array of zernike coeffs on which the phase can be decomposed
// an array of zernike polynomials (carma_object)
// a matrix to decompose the phase on zernike coeffs

#include <carma.hpp>
#include <carma_obj.hpp>

class SutraPhase {
 public:
  CarmaContext *current_context;
  int32_t device;

  CarmaObj<float> *d_screen;
  int64_t screen_size;
  float *zer_coeff;
  CarmaObj<float> *zernikes;
  CarmaObj<float> *mat;

 public:
  SutraPhase(CarmaContext *current_context, int64_t size);
  SutraPhase(const SutraPhase &phase);
  ~SutraPhase();
};

#endif  // _SUTRA_PHASE_H_
