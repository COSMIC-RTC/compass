// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_atmos.h
//! \ingroup   libsutra
//! \class     SutraAtmos
//! \brief     this class provides the atmos features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#ifndef _SUTRA_ATMOS_H_
#define _SUTRA_ATMOS_H_

#include <sutra_tscreen.h>

using std::pair;
using std::vector;

class SutraAtmos {
 public:
  int nscreens;
  vector<SutraTurbuScreen *> d_screens;
  float r0;
  CarmaContext *current_context;

 public:
  SutraAtmos(CarmaContext *context, int nscreens, float global_r0, float *r0,
              long *size, long *size2, float *altitude, float *windspeed,
              float *winddir, float *deltax, float *deltay, int device);
  ~SutraAtmos();

  int init_screen(int idx, float *h_A, float *h_B, unsigned int *h_istencilx,
                  unsigned int *h_istencily, int seed);

  int add_screen(float altitude, long size, long stencilSize, float amplitude,
                 float windspeed, float winddir, float deltax, float deltay,
                 int device);
  int del_screen(const int idx);
  int refresh_screen(int idx);

  int move_atmos();
  int set_r0(float r0);
  int set_seed(int idx, float seed);
};

#endif  // _SUTRA_ATMOS_H_
