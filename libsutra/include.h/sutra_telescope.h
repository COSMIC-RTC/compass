// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_telescope.h
//! \ingroup   libsutra
//! \class     SutraTelescope
//! \brief     this class provides the telescope features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#ifndef _SUTRA_TELESCOPE_H_
#define _SUTRA_TELESCOPE_H_

#include <carma.h>
#include <sutra_phase.h>
#include <map>
#include <string>
#include <vector>

using std::string;

class SutraTelescope {
 public:
  CarmaContext *current_context;
  int device;  // device #

  long pup_size;       // size of pupil
  long num_eleme_pup;  // number of points in the pupil

  CarmaObj<float> *d_pupil;        // the pupil mask
  CarmaObj<float> *d_phase_ab_M1;  // the phase aberration for M1

  long pup_size_m;  // size of pupil

  CarmaObj<float> *d_pupil_m;        // the pupil mask
  CarmaObj<float> *d_phase_ab_M1_m;  // the phase aberration for M1

  CarmaObj<float> *d_input_phase; // Circular buffer of phase screens to be played
  int input_phase_counter; // Current input phase screen in the circular buffer


 public:
  SutraTelescope(CarmaContext *context, long pup_size, long num_eleme_pup,
                  float *pupil, long pup_size_m, float *pupil_m);
  ~SutraTelescope();
  int set_phase_ab_M1(float *phase_ab_M1, int size);
  int set_phase_ab_M1_m(float *phase_ab_M1_m, int size);
  /**
   * @brief Set a 3D cube of phase screens to be played. Each phase screen 
   * is shown to sources as an additional layer to be raytraced. 
   * 
   * @param input_phase Cube of phase screens
   * @param size 1 phase screen size. Must be equal to d_pupil_m size
   * @param N Number of phase screens in the cube
   * @return int Success status
   */
  int set_input_phase(float *input_phase, int size, int N);
  /**
   * @brief Update input_phase_counter to take the next phase screen in
   * the circular buffer d_input_phase
   * 
   * @return int Success status
   */
  int update_input_phase();
  /**
   * @brief Reset circular buffer d_input_phase
   * 
   * @return int Success status
   */
  int reset_input_phase();
};

#endif  // _SUTRA_TELESCOPE_H_
