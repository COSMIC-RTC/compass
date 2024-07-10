// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_atmos.cpp
//! \ingroup   libsutra
//! \class     SutraAtmos
//! \brief     this class provides the atmos features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_atmos.hpp>
#include <algorithm>

SutraAtmos::SutraAtmos(CarmaContext *context, int32_t nscreens, float global_r0,
                         float *r0, int64_t *size, int64_t *stencilSize,
                         float *altitude, float *windspeed, float *winddir,
                         float *deltax, float *deltay, int32_t device) {
  this->nscreens = nscreens;
  // this->r0       = r0;
  this->current_context = context;
  this->r0 = global_r0;

  for (int32_t i = 0; i < nscreens; i++) {
    d_screens.push_back(new SutraTurbuScreen(
        context, size[i], stencilSize[i], r0[i], altitude[i], windspeed[i],
        winddir[i], deltax[i], deltay[i], device));
  }
}

SutraAtmos::~SutraAtmos() {
  for (vector<SutraTurbuScreen *>::iterator it = this->d_screens.begin();
       this->d_screens.end() != it; it++) {
    delete *it;
  }

  this->d_screens.clear();
  // d_screens.erase(d_screens.begin(),d_screens.end());
}

int32_t SutraAtmos::init_screen(int32_t idx, float *h_A, float *h_B,
                             uint32_t *h_istencilx,
                             uint32_t *h_istencily, int32_t seed) {
  if (idx < this->d_screens.size()) {
    d_screens[idx]->init_screen(h_A, h_B, h_istencilx, h_istencily, seed);
    d_screens[idx]->refresh_screen();
  } else
    DEBUG_TRACE("Index exceed vector size");

  return EXIT_SUCCESS;
}

int32_t SutraAtmos::refresh_screen(int32_t idx) {
  if (idx < this->d_screens.size())
    this->d_screens[idx]->refresh_screen();
  else
    DEBUG_TRACE("Index exceed vector size");
  return EXIT_SUCCESS;
}

int32_t SutraAtmos::add_screen(float alt, int64_t size, int64_t stencilSize,
                            float r0_thislayer, float windspeed, float winddir,
                            float deltax, float deltay, int32_t device) {
  this->d_screens.push_back(
      new SutraTurbuScreen(current_context, size, stencilSize, r0_thislayer, alt,
                        windspeed, winddir, deltax, deltay, device));
  this->r0 = powf(powf(r0, -5.0f / 3.0f) + powf(r0_thislayer, -5.0f / 3.0f),
                  -3.0f / 5.0f);
  this->nscreens++;

  return EXIT_SUCCESS;
}

int32_t SutraAtmos::del_screen(int32_t idx) {
  if (idx < this->d_screens.size()) {
    this->nscreens--;
    this->r0 =
        powf(powf(r0, -5.0f / 3.0f) - powf(d_screens[idx]->r0, -5.0f / 3.0f),
             -3.0f / 5.0f);
    delete this->d_screens[idx];
    this->d_screens.erase(this->d_screens.begin() + idx);
  }

  else
    DEBUG_TRACE("Index exceed vector size");

  return EXIT_SUCCESS;
}

int32_t SutraAtmos::move_atmos() {
  vector<SutraTurbuScreen *>::iterator p;
  p = this->d_screens.begin();

  while (p != this->d_screens.end()) {
    (*p)->accumx += (*p)->deltax;
    (*p)->accumy += (*p)->deltay;

    int32_t deltax = (int32_t)(*p)->accumx;
    int32_t deltay = (int32_t)(*p)->accumy;
    int32_t cx = deltax > 0 ? 1 : -1;
    int32_t cy = deltay > 0 ? 1 : -1;
    for (int32_t cc = 0; cc < cx * deltax; cc++) (*p)->extrude(1 * cx);
    (*p)->accumx -= deltax;
    for (int32_t cc = 0; cc < cy * deltay; cc++) (*p)->extrude(2 * cy);
    (*p)->accumy -= deltay;
    ++p;
  }
  return EXIT_SUCCESS;
}

int32_t SutraAtmos::set_r0(float r0) {
  // this->amplitude = powf(r0, -5.0f / 6.0f)
  float scaling = powf(r0 / this->r0, -5.0f / 6.0f);
  for (vector<SutraTurbuScreen *>::iterator it = this->d_screens.begin();
       this->d_screens.end() != it; it++) {
    (*it)->r0 *= this->r0 / r0;
    (*it)->amplitude *= scaling;
  }
  this->r0 = r0;
  return EXIT_SUCCESS;
}

int32_t SutraAtmos::set_seed(int32_t idx, float seed) {
  if (idx < this->d_screens.size())
    this->d_screens[idx]->set_seed(seed);
  else
    DEBUG_TRACE("Index exceed vector size");

  return EXIT_SUCCESS;
}
