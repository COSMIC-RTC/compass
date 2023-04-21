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
//! \version   5.4.3
//! \date      2022/01/24

#include <sutra_atmos.h>
#include <algorithm>

SutraAtmos::SutraAtmos(CarmaContext *context, int nscreens, float global_r0,
                         float *r0, long *size, long *stencilSize,
                         float *altitude, float *windspeed, float *winddir,
                         float *deltax, float *deltay, int device) {
  this->nscreens = nscreens;
  // this->r0       = r0;
  this->current_context = context;
  this->r0 = global_r0;

  for (int i = 0; i < nscreens; i++) {
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

int SutraAtmos::init_screen(int idx, float *h_A, float *h_B,
                             unsigned int *h_istencilx,
                             unsigned int *h_istencily, int seed) {
  if (idx < this->d_screens.size()) {
    d_screens[idx]->init_screen(h_A, h_B, h_istencilx, h_istencily, seed);
    d_screens[idx]->refresh_screen();
  } else
    DEBUG_TRACE("Index exceed vector size");

  return EXIT_SUCCESS;
}

int SutraAtmos::refresh_screen(int idx) {
  if (idx < this->d_screens.size())
    this->d_screens[idx]->refresh_screen();
  else
    DEBUG_TRACE("Index exceed vector size");
  return EXIT_SUCCESS;
}

int SutraAtmos::add_screen(float alt, long size, long stencilSize,
                            float r0_thislayer, float windspeed, float winddir,
                            float deltax, float deltay, int device) {
  this->d_screens.push_back(
      new SutraTurbuScreen(current_context, size, stencilSize, r0_thislayer, alt,
                        windspeed, winddir, deltax, deltay, device));
  this->r0 = powf(powf(r0, -5.0f / 3.0f) + powf(r0_thislayer, -5.0f / 3.0f),
                  -3.0f / 5.0f);
  this->nscreens++;

  return EXIT_SUCCESS;
}

int SutraAtmos::del_screen(int idx) {
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

int SutraAtmos::move_atmos() {
  vector<SutraTurbuScreen *>::iterator p;
  p = this->d_screens.begin();

  while (p != this->d_screens.end()) {
    (*p)->accumx += (*p)->deltax;
    (*p)->accumy += (*p)->deltay;

    int deltax = (int)(*p)->accumx;
    int deltay = (int)(*p)->accumy;
    int cx = deltax > 0 ? 1 : -1;
    int cy = deltay > 0 ? 1 : -1;
    for (int cc = 0; cc < cx * deltax; cc++) (*p)->extrude(1 * cx);
    (*p)->accumx -= deltax;
    for (int cc = 0; cc < cy * deltay; cc++) (*p)->extrude(2 * cy);
    (*p)->accumy -= deltay;
    ++p;
  }
  return EXIT_SUCCESS;
}

int SutraAtmos::set_r0(float r0) {
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

int SutraAtmos::set_seed(int idx, float seed) {
  if (idx < this->d_screens.size())
    this->d_screens[idx]->set_seed(seed);
  else
    DEBUG_TRACE("Index exceed vector size");

  return EXIT_SUCCESS;
}
