// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_roket.h
//! \ingroup   libsutra
//! \class     SutraRoket
//! \brief     this class provides the roket features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_ROKET_H_
#define _SUTRA_ROKET_H_

#include <carma.h>
#include <carma_obj.h>
#include <sutra_atmos.h>
#include <sutra_rtc.h>
#include <sutra_sensors.h>
#include <sutra_target.h>
#include <sutra_telescope.h>
#include <sutra_utils.h>

class SutraRoket {
 public:
  CarmaContext *current_context;
  int32_t device;
  float gain;  // Loop gain
  int32_t nfilt;   // Number of filtered modes
  int32_t nactus;  // number of actuators
  int32_t nmodes;  // number of modes
  int32_t iterk;   // current iteration number
  int32_t niter;
  int32_t loopcontroller;  // index of the loop controller
  int32_t geocontroller;   // index of the geo controller
  int32_t nslopes;
  // sutra objects to supervise
  SutraRtc *rtc;
  SutraSensors *sensors;
  SutraTarget *target;
  SutraTelescope *tel;
  SutraAtmos *atm;
  SutraDms *dms;
  SutraControllerLs *loopcontrol;
  SutraControllerGeo *geocontrol;

  // Projection matrices
  CarmaObj<float> *d_P;
  CarmaObj<float> *d_Btt;

  // Error contributors buffers
  CarmaObj<float> *d_noise;
  CarmaObj<float> *d_nonlinear;
  CarmaObj<float> *d_tomo;
  CarmaObj<float> *d_filtered;
  CarmaObj<float> *d_alias;
  CarmaObj<float> *d_bandwidth;
  float fitting;

  // Residual error buffers
  CarmaObj<float> *d_fullErr;
  CarmaObj<float> *d_err1;
  CarmaObj<float> *d_err2;
  // Command loop backup
  CarmaObj<float> *d_bkup_com;
  // Target screen backup
  CarmaObj<float> *d_bkup_screen;
  // Additional buffers
  CarmaObj<float> *d_commanded;
  CarmaObj<float> *d_modes;
  CarmaObj<float> *d_filtmodes;
  CarmaObj<float> *d_tmpdiff;
  // Loop filter matrix
  CarmaObj<float> *d_gRD;
  // R*D matrix
  CarmaObj<float> *d_RD;
  // PSF ortho
  CarmaObj<float> *d_psfortho;
  // Loop output covariance matrices
  CarmaObj<float> *d_covv;
  CarmaObj<float> *d_covm;

 public:
  SutraRoket(CarmaContext *context, int32_t device, SutraRtc *rtc,
              SutraSensors *sensors, SutraTarget *target, SutraDms *dms,
              SutraTelescope *tel, SutraAtmos *atm, int32_t loopcontroller,
              int32_t geocontroller, int32_t nactus, int32_t nmodes, int32_t nfilt, int32_t niter,
              float *Btt, float *P, float *gRD, float *RD);
  ~SutraRoket();

  int32_t compute_breakdown();
  int32_t save_loop_state();
  int32_t restore_loop_state();
  int32_t apply_loop_filter(CarmaObj<float> *d_odata, CarmaObj<float> *d_idata1,
                        CarmaObj<float> *d_idata2, float gain, int32_t k);
};

int32_t separate_modes(float *modes, float *filtmodes, int32_t nmodes, int32_t nfilt,
                   CarmaDevice *device);

#endif
