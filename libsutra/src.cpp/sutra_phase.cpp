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

//! \file      sutra_phase.cpp
//! \ingroup   libsutra
//! \class     SutraPhase
//! \brief     this class provides the phase features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_phase.hpp>

SutraPhase::SutraPhase(CarmaContext *current_context, int64_t size) {
  this->current_context = current_context;
  this->screen_size = size;
  this->mat = 0;
  this->zernikes = 0;
  this->zer_coeff = 0;
  this->device = current_context->get_active_device();

  int64_t *dims_data2 = new int64_t[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->screen_size;
  dims_data2[2] = this->screen_size;

  this->d_screen = new CarmaObj<float>(current_context, dims_data2);
  delete[] dims_data2;
}

SutraPhase::~SutraPhase() {
  current_context->set_active_device(device, 1);
  delete this->d_screen;
}
