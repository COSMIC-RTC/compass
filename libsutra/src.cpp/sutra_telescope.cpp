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

//! \file      sutra_telescope.cpp
//! \ingroup   libsutra
//! \class     SutraTelescope
//! \brief     this class provides the telescope features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_phase.hpp>
#include <sutra_telescope.hpp>
#include <sutra_utils.hpp>

SutraTelescope::SutraTelescope(CarmaContext *current_context, int64_t n_pup,
                                 int64_t npos, float *pupil, int64_t n_pup_m,
                                 float *pupil_m) {
  this->current_context = current_context;
  this->device = current_context->get_active_device();

  this->pup_size = n_pup;
  this->pup_size_m = n_pup_m;

  this->num_eleme_pup = npos;
  this->d_phase_ab_M1 = nullptr;
  this->d_phase_ab_M1_m = nullptr;

  int64_t *dims_data2 = new int64_t[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->pup_size;
  dims_data2[2] = this->pup_size;

  this->d_pupil = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_pupil->host2device(pupil);

  int64_t *dims_data3 = new int64_t[3];
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

int32_t SutraTelescope::set_phase_ab_M1(float *phase_ab_M1, int32_t size) {
  current_context->set_active_device(device, 1);
  if (size == this->pup_size * this->pup_size) {
    int64_t *dims_data2 = new int64_t[3];
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

int32_t SutraTelescope::set_phase_ab_M1_m(float *phase_ab_M1_m, int32_t size) {
  current_context->set_active_device(device, 1);
  if (size == this->pup_size_m * this->pup_size_m) {
    int64_t *dims_data2 = new int64_t[3];
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

int32_t SutraTelescope::set_input_phase(float *input_phase, int32_t size, int32_t N) {
  current_context->set_active_device(device, 1);
  if (size == this->pup_size_m * this->pup_size_m) {
    int64_t *dims_data3 = new int64_t[4];
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

int32_t SutraTelescope::update_input_phase() {
  if (this->d_input_phase != nullptr) {
    this->input_phase_counter++;
    if(this->input_phase_counter >= this->d_input_phase->get_dims(3))
      this->input_phase_counter = 0;
  }

  return EXIT_SUCCESS;
}

int32_t SutraTelescope::reset_input_phase() {
  if (this->d_input_phase != nullptr) {
    delete this->d_input_phase;
    this->d_input_phase = nullptr;
    this->input_phase_counter = 0;
  }

  return EXIT_SUCCESS;
}