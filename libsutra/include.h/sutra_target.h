// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_target.h
//! \ingroup   libsutra
//! \class     SutraTarget
//! \brief     this class provides the target features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_TARGET_H_
#define _SUTRA_TARGET_H_

#include <sutra_source.h>
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
