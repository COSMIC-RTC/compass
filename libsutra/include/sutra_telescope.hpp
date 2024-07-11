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

//! \file      sutra_telescope.hpp
//! \ingroup   libsutra
//! \class     SutraTelescope
//! \brief     this class provides the telescope features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_TELESCOPE_H_
#define _SUTRA_TELESCOPE_H_

#include <carma.hpp>
#include <sutra_phase.hpp>
#include <map>
#include <string>
#include <vector>

using std::string;

class SutraTelescope {
 public:
  CarmaContext *current_context;
  int32_t device;  // device #

  int64_t pup_size;       // size of pupil
  int64_t num_eleme_pup;  // number of points in the pupil

  CarmaObj<float> *d_pupil;        // the pupil mask
  CarmaObj<float> *d_phase_ab_M1;  // the phase aberration for M1

  int64_t pup_size_m;  // size of pupil

  CarmaObj<float> *d_pupil_m;        // the pupil mask
  CarmaObj<float> *d_phase_ab_M1_m;  // the phase aberration for M1

  CarmaObj<float> *d_input_phase; // Circular buffer of phase screens to be played
  int32_t input_phase_counter; // Current input phase screen in the circular buffer


 public:
  SutraTelescope(CarmaContext *context, int64_t pup_size, int64_t num_eleme_pup,
                  float *pupil, int64_t pup_size_m, float *pupil_m);
  ~SutraTelescope();
  int32_t set_phase_ab_M1(float *phase_ab_M1, int32_t size);
  int32_t set_phase_ab_M1_m(float *phase_ab_M1_m, int32_t size);
  /**
   * @brief Set a 3D cube of phase screens to be played. Each phase screen
   * is shown to sources as an additional layer to be raytraced.
   *
   * @param input_phase Cube of phase screens
   * @param size 1 phase screen size. Must be equal to d_pupil_m size
   * @param N Number of phase screens in the cube
   * @return int32_t Success status
   */
  int32_t set_input_phase(float *input_phase, int32_t size, int32_t N);
  /**
   * @brief Update input_phase_counter to take the next phase screen in
   * the circular buffer d_input_phase
   *
   * @return int32_t Success status
   */
  int32_t update_input_phase();
  /**
   * @brief Reset circular buffer d_input_phase
   *
   * @return int32_t Success status
   */
  int32_t reset_input_phase();
};

#endif  // _SUTRA_TELESCOPE_H_
