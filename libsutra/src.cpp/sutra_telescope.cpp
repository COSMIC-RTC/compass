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

//! \file      sutra_telescope.cpp
//! \ingroup   libsutra
//! \class     SutraTelescope
//! \brief     this class provides the telescope features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.1
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <sutra_phase.h>
#include <sutra_telescope.h>
#include <sutra_utils.h>

SutraTelescope::SutraTelescope(CarmaContext *current_context, long n_pup,
                                 long npos, float *pupil, long n_pup_m,
                                 float *pupil_m) {
  this->current_context = current_context;
  this->device = current_context->get_active_device();

  this->pup_size = n_pup;
  this->pup_size_m = n_pup_m;

  this->num_eleme_pup = npos;
  this->d_phase_ab_M1 = nullptr;
  this->d_phase_ab_M1_m = nullptr;

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->pup_size;
  dims_data2[2] = this->pup_size;

  this->d_pupil = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_pupil->host2device(pupil);

  long *dims_data3 = new long[3];
  dims_data3[0] = 2;
  dims_data3[1] = this->pup_size_m;
  dims_data3[2] = this->pup_size_m;

  this->d_pupil_m = new CarmaObj<float>(this->current_context, dims_data3);
  this->d_pupil_m->host2device(pupil_m);

  this->d_input_phase = nullptr;
  this->input_phase_counter = 0;

  delete[] dims_data2;
  delete[] dims_data3;
}

SutraTelescope::~SutraTelescope() {
  // delete this->current_context;
  current_context->set_active_device(device, 1);
  delete this->d_pupil;
  delete this->d_pupil_m;
  if (this->d_phase_ab_M1 != nullptr) delete this->d_phase_ab_M1;
  if (this->d_phase_ab_M1_m != nullptr) delete this->d_phase_ab_M1_m;
}

int SutraTelescope::set_phase_ab_M1(float *phase_ab_M1, int size) {
  current_context->set_active_device(device, 1);
  if (size == this->pup_size * this->pup_size) {
    long *dims_data2 = new long[3];
    dims_data2[0] = 2;
    dims_data2[1] = this->pup_size;
    dims_data2[2] = this->pup_size;
    if (this->d_phase_ab_M1 == nullptr)
      this->d_phase_ab_M1 =
          new CarmaObj<float>(this->current_context, dims_data2, phase_ab_M1);
    else
      this->d_phase_ab_M1->host2device(phase_ab_M1);
  } else
    DEBUG_TRACE("Wrong dimensions");

  return EXIT_SUCCESS;
}

int SutraTelescope::set_phase_ab_M1_m(float *phase_ab_M1_m, int size) {
  current_context->set_active_device(device, 1);
  if (size == this->pup_size_m * this->pup_size_m) {
    long *dims_data2 = new long[3];
    dims_data2[0] = 2;
    dims_data2[1] = this->pup_size_m;
    dims_data2[2] = this->pup_size_m;
    if (this->d_phase_ab_M1_m == nullptr)
      this->d_phase_ab_M1_m = new CarmaObj<float>(this->current_context,
                                                   dims_data2, phase_ab_M1_m);
    else
      this->d_phase_ab_M1_m->host2device(phase_ab_M1_m);
  } else
    DEBUG_TRACE("Wrong dimensions");

  return EXIT_SUCCESS;
}

int SutraTelescope::set_input_phase(float *input_phase, int size, int N) {
  current_context->set_active_device(device, 1);
  if (size == this->pup_size_m * this->pup_size_m) {
    long *dims_data3 = new long[4];
    dims_data3[0] = 3;
    dims_data3[1] = this->pup_size_m;
    dims_data3[2] = this->pup_size_m;
    dims_data3[3] = N;
    if (this->d_input_phase != nullptr)
      delete this->d_input_phase;
    
    this->d_input_phase = new CarmaObj<float>(this->current_context,
                                                   dims_data3, input_phase);
    this->input_phase_counter = 0;

  } else {
      DEBUG_TRACE("Wrong dimensions");
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int SutraTelescope::update_input_phase() {
  if (this->d_input_phase != nullptr) {
    this->input_phase_counter++;
    if(this->input_phase_counter >= this->d_input_phase->get_dims(3))
      this->input_phase_counter = 0;
  }

  return EXIT_SUCCESS;
}

int SutraTelescope::reset_input_phase() {
  if (this->d_input_phase != nullptr) {
    delete this->d_input_phase;
    this->d_input_phase = nullptr;
    this->input_phase_counter = 0;
  }

  return EXIT_SUCCESS;
}