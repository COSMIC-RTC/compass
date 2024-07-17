// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_tscreen.hpp
//! \ingroup   libsutra
//! \class     SutraTurbuScreen
//! \brief     this class provides the turbulent screen features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_TSCREEN_H_
#define _SUTRA_TSCREEN_H_

#include <sutra_phase.hpp>
#include <sutra_utils.hpp>

class SutraTurbuScreen {
 public:
  int32_t device;              // The device #
  SutraPhase *d_tscreen;  // The phase screen
  CarmaObj<float>
      *d_tscreen_o;  // Additional space of the same size as the phase screen
  CarmaObj<float> *d_mat_a;                 // A matrix for extrusion
  CarmaObj<float> *d_mat_b;                 // B matrix for extrusion
  CarmaObj<uint32_t> *d_istencilx;  // stencil for column extrusion
  CarmaObj<uint32_t> *d_istencily;  // stencil for row extrusion
  CarmaObj<float> *d_z;                 // tmp array for extrusion process
  CarmaObj<float> *d_noise;             // tmp array containing random numbers
  CarmaObj<float> *d_ytmp;  // contains the extrude update (row or column)
  int64_t screen_size;          // size of phase screens
  float r0;                  // layer r0 (pixel units)
  float amplitude;           // amplitude for extrusion (r0^-5/6)
  float altitude;
  float windspeed;
  float winddir;
  float deltax;  // number of rows to extrude per iteration
  float deltay;  // number of lines to extrude per iteration
  // internal
  float accumx;
  float accumy;
  cudaChannelFormatDesc
      channel_desc;  // Channel descriptor for texture memory access

  CarmaObj<cuFloatComplex>
      *d_tscreen_c;  // Additional space for von karman screen generation
  float norm_vk;
  bool vk_on;
  CarmaContext *current_context;

 public:
  SutraTurbuScreen(CarmaContext *context, int64_t size, int64_t size2, float amplitude,
                float altitude, float windspeed, float winddir, float deltax,
                float deltay, int32_t device);
  // SutraTurbuScreen(const SutraTurbuScreen &tscreen);
  ~SutraTurbuScreen();

  int32_t init_screen(float *h_A, float *h_B, uint32_t *h_istencilx,
                  uint32_t *h_istencily, int32_t seed);
  int32_t extrude(int32_t dir);
  int32_t init_vk(int32_t seed, int32_t pupd);
  int32_t generate_vk(float l0, int32_t nalias);
  int32_t refresh_screen();
  int32_t set_seed(int32_t seed);
  int32_t set_deltax(float deltax);
  int32_t set_deltay(float deltay);
  int32_t set_istencilx(uint32_t* istencil);
  int32_t set_istencily(uint32_t* istencil);
};

int32_t gene_vonkarman(cuFloatComplex *d_odata, float *d_idata, float k0,
                   int32_t nalias, int32_t nx, int32_t ny, int32_t block_size);
int32_t norm_pscreen(float *d_odata, float *d_idata, int32_t nx, int32_t ny,
                 float norm_fact, CarmaDevice *device);

#endif  // _SUTRA_TSCREEN_H_
