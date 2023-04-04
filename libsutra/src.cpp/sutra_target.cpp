// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_target.cpp
//! \ingroup   libsutra
//! \class     SutraTarget
//! \brief     this class provides the target features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24

#include <sutra_target.h>

SutraTarget::SutraTarget(CarmaContext *context, SutraTelescope *yTelescope,
                           int ntargets, float *xpos, float *ypos,
                           float *lambda, float *mag, float zerop, long *sizes,
                           int Npts, int device) {
  this->ntargets = ntargets;

  for (int i = 0; i < ntargets; i++) {
    d_targets.push_back(new SutraSource(context, xpos[i], ypos[i], lambda[i],
                                         mag[i], zerop, sizes[i], "target",
                                         yTelescope->d_pupil, Npts, device));
  }
}

SutraTarget::~SutraTarget() {
  //  for (size_t idx = 0; idx < (this->d_targets).size(); idx++) {
  while ((this->d_targets).size() > 0) {
    delete this->d_targets.back();
    d_targets.pop_back();
  }
}
