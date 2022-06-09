// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_tscreen.h
//! \ingroup   libsutra
//! \class     SutraTurbuScreen
//! \brief     this class provides the turbulent screen features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_TSCREEN_H_
#define _SUTRA_TSCREEN_H_

#include <sutra_phase.h>
#include <sutra_utils.h>

class SutraTurbuScreen {
 public:
  int device;              // The device #
  SutraPhase *d_tscreen;  // The phase screen
  CarmaObj<float>
      *d_tscreen_o;  // Additional space of the same size as the phase screen
  CarmaObj<float> *d_mat_a;                 // A matrix for extrusion
  CarmaObj<float> *d_mat_b;                 // B matrix for extrusion
  CarmaObj<unsigned int> *d_istencilx;  // stencil for column extrusion
  CarmaObj<unsigned int> *d_istencily;  // stencil for row extrusion
  CarmaObj<float> *d_z;                 // tmp array for extrusion process
  CarmaObj<float> *d_noise;             // tmp array containing random numbers
  CarmaObj<float> *d_ytmp;  // contains the extrude update (row or column)
  long screen_size;          // size of phase screens
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
  SutraTurbuScreen(CarmaContext *context, long size, long size2, float amplitude,
                float altitude, float windspeed, float winddir, float deltax,
                float deltay, int device);
  // SutraTurbuScreen(const SutraTurbuScreen &tscreen);
  ~SutraTurbuScreen();

  int init_screen(float *h_A, float *h_B, unsigned int *h_istencilx,
                  unsigned int *h_istencily, int seed);
  int extrude(int dir);
  int init_vk(int seed, int pupd);
  int generate_vk(float l0, int nalias);
  int refresh_screen();
  int set_seed(int seed);
  int set_deltax(float deltax);
  int set_deltay(float deltay);
  int set_istencilx(unsigned int* istencil);
  int set_istencily(unsigned int* istencil);
};

int gene_vonkarman(cuFloatComplex *d_odata, float *d_idata, float k0,
                   int nalias, int nx, int ny, int block_size);
int norm_pscreen(float *d_odata, float *d_idata, int nx, int ny,
                 float norm_fact, CarmaDevice *device);

#endif  // _SUTRA_TSCREEN_H_
