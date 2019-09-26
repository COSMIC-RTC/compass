// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      sutra_roket.h
//! \ingroup   libsutra
//! \class     sutra_roket
//! \brief     this class provides the roket features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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

class sutra_roket {
 public:
  carma_context *current_context;
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
  sutra_rtc *rtc;
  sutra_sensors *sensors;
  sutra_target *target;
  sutra_telescope *tel;
  sutra_atmos *atm;
  sutra_dms *dms;
  sutra_controller_ls *loopcontrol;
  sutra_controller_geo *geocontrol;

  // Projection matrices
  carma_obj<float> *d_P;
  carma_obj<float> *d_Btt;

  // Error contributors buffers
  carma_obj<float> *d_noise;
  carma_obj<float> *d_nonlinear;
  carma_obj<float> *d_tomo;
  carma_obj<float> *d_filtered;
  carma_obj<float> *d_alias;
  carma_obj<float> *d_bandwidth;
  float fitting;

  // Residual error buffers
  carma_obj<float> *d_fullErr;
  carma_obj<float> *d_err1;
  carma_obj<float> *d_err2;
  // Command loop backup
  carma_obj<float> *d_bkup_com;
  // Target screen backup
  carma_obj<float> *d_bkup_screen;
  // Additional buffers
  carma_obj<float> *d_commanded;
  carma_obj<float> *d_modes;
  carma_obj<float> *d_filtmodes;
  carma_obj<float> *d_tmpdiff;
  // Loop filter matrix
  carma_obj<float> *d_gRD;
  // R*D matrix
  carma_obj<float> *d_RD;
  // PSF ortho
  carma_obj<float> *d_psfortho;
  // Loop output covariance matrices
  carma_obj<float> *d_covv;
  carma_obj<float> *d_covm;

 public:
  sutra_roket(carma_context *context, int device, sutra_rtc *rtc,
              sutra_sensors *sensors, sutra_target *target, sutra_dms *dms,
              sutra_telescope *tel, sutra_atmos *atm, int loopcontroller,
              int geocontroller, int nactus, int nmodes, int nfilt, int niter,
              float *Btt, float *P, float *gRD, float *RD);
  ~sutra_roket();

  int compute_breakdown();
  int save_loop_state();
  int restore_loop_state();
  int apply_loop_filter(carma_obj<float> *d_odata, carma_obj<float> *d_idata1,
                        carma_obj<float> *d_idata2, float gain, int k);
};

int separate_modes(float *modes, float *filtmodes, int nmodes, int nfilt,
                   carma_device *device);

#endif
