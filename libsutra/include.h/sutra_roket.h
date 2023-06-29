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
//! \version   5.4.4
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
  int device;
  float gain;  // Loop gain
  int nfilt;   // Number of filtered modes
  int nactus;  // number of actuators
  int nmodes;  // number of modes
  int iterk;   // current iteration number
  int niter;
  int loopcontroller;  // index of the loop controller
  int geocontroller;   // index of the geo controller
  int nslopes;
  // sutra objects to supervise
  SutraRtc *rtc;
  SutraSensors *sensors;
  SutraTarget *target;
  SutraTelescope *tel;
  SutraAtmos *atm;
  SutraDms *dms;
  sutra_controller_ls *loopcontrol;
  sutra_controller_geo *geocontrol;

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
  SutraRoket(CarmaContext *context, int device, SutraRtc *rtc,
              SutraSensors *sensors, SutraTarget *target, SutraDms *dms,
              SutraTelescope *tel, SutraAtmos *atm, int loopcontroller,
              int geocontroller, int nactus, int nmodes, int nfilt, int niter,
              float *Btt, float *P, float *gRD, float *RD);
  ~SutraRoket();

  int compute_breakdown();
  int save_loop_state();
  int restore_loop_state();
  int apply_loop_filter(CarmaObj<float> *d_odata, CarmaObj<float> *d_idata1,
                        CarmaObj<float> *d_idata2, float gain, int k);
};

int separate_modes(float *modes, float *filtmodes, int nmodes, int nfilt,
                   CarmaDevice *device);

#endif
