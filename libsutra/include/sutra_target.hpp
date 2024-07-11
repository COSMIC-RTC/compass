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

//! \file      sutra_target.hpp
//! \ingroup   libsutra
//! \class     SutraTarget
//! \brief     this class provides the target features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_TARGET_H_
#define _SUTRA_TARGET_H_

#include <sutra_source.hpp>
#include <map>
#include <string>
#include <vector>

using std::vector;

typedef std::pair<std::string, int32_t> type_screen;

class SutraTarget {
 public:
  int32_t ntargets;
  vector<SutraSource *> d_targets;

 public:
  SutraTarget(CarmaContext *context, SutraTelescope *d_tel, int32_t ntargets,
               float *xpos, float *ypos, float *lambda, float *mag, float zerop,
               int64_t *sizes, int32_t Npts, int32_t device);
  ~SutraTarget();

  int32_t get_phase(int32_t ntarget, float *dest);
};

#endif  // _SUTRA_TARGET_H_
