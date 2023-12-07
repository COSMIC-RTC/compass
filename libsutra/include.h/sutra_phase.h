// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_phase.h
//! \ingroup   libsutra
//! \class     SutraPhase
//! \brief     this class provides the phase features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
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

#include <carma.h>
#include <carma_obj.h>

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
