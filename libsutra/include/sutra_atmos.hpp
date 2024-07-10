// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_atmos.hpp
//! \ingroup   libsutra
//! \class     SutraAtmos
//! \brief     this class provides the atmos features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_ATMOS_H_
#define _SUTRA_ATMOS_H_

#include <sutra_tscreen.hpp>

using std::pair;
using std::vector;

class SutraAtmos {
 public:
  int32_t nscreens;
  vector<SutraTurbuScreen *> d_screens;
  float r0;
  CarmaContext *current_context;

 public:
  SutraAtmos(CarmaContext *context, int32_t nscreens, float global_r0, float *r0,
              int64_t *size, int64_t *size2, float *altitude, float *windspeed,
              float *winddir, float *deltax, float *deltay, int32_t device);
  ~SutraAtmos();

  int32_t init_screen(int32_t idx, float *h_A, float *h_B, uint32_t *h_istencilx,
                  uint32_t *h_istencily, int32_t seed);

  int32_t add_screen(float altitude, int64_t size, int64_t stencilSize, float amplitude,
                 float windspeed, float winddir, float deltax, float deltay,
                 int32_t device);
  int32_t del_screen(const int32_t idx);
  int32_t refresh_screen(int32_t idx);

  int32_t move_atmos();
  int32_t set_r0(float r0);
  int32_t set_seed(int32_t idx, float seed);
};

#endif  // _SUTRA_ATMOS_H_
