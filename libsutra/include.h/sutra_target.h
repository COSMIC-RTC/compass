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
//! \version   5.4.0
//! \date      2022/01/24

#ifndef _SUTRA_TARGET_H_
#define _SUTRA_TARGET_H_

#include <sutra_source.h>
#include <map>
#include <string>
#include <vector>

using std::vector;

typedef std::pair<std::string, int> type_screen;

class SutraTarget {
 public:
  int ntargets;
  vector<SutraSource *> d_targets;

 public:
  SutraTarget(CarmaContext *context, SutraTelescope *d_tel, int ntargets,
               float *xpos, float *ypos, float *lambda, float *mag, float zerop,
               long *sizes, int Npts, int device);
  ~SutraTarget();

  int get_phase(int ntarget, float *dest);
};

#endif  // _SUTRA_TARGET_H_
