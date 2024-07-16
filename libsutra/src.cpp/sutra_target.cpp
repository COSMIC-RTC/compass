// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_target.cpp
//! \ingroup   libsutra
//! \class     SutraTarget
//! \brief     this class provides the target features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_target.hpp>

SutraTarget::SutraTarget(CarmaContext *context, SutraTelescope *yTelescope,
                           int32_t ntargets, float *xpos, float *ypos,
                           float *lambda, float *mag, float zerop, int64_t *sizes,
                           int32_t Npts, int32_t device) {
  this->ntargets = ntargets;

  for (int32_t i = 0; i < ntargets; i++) {
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
