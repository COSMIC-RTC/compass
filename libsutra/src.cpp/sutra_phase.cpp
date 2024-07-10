// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_phase.cpp
//! \ingroup   libsutra
//! \class     SutraPhase
//! \brief     this class provides the phase features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
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
