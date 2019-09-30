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

//! \file      sutra_sensors.h
//! \ingroup   libsutra
//! \class     sutra_sensors
//! \brief     this class provides the sensors features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_SENSORS_H_
#define _SUTRA_SENSORS_H_

#include <carma_utils.h>
#include <sutra_telescope.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <sutra_wfs_sh.h>

#include <map>
#include <vector>

using std::map;
using std::string;
using std::vector;

class sutra_sensors {
 public:
  int device;
  bool roket;
  carma_context *current_context;
  size_t nsensors() { return d_wfs.size(); }
  vector<sutra_wfs *> d_wfs;
  map<vector<int>, cufftHandle *> campli_plans;
  map<vector<int>, cufftHandle *> fttotim_plans;
  map<vector<int>, cufftHandle *> ftlgskern_plans;

  carma_obj<cuFloatComplex> *d_camplipup;
  carma_obj<cuFloatComplex> *d_camplifoc;
  carma_obj<cuFloatComplex> *d_fttotim;
  carma_obj<cuFloatComplex> *d_ftlgskern;
  carma_obj<float> *d_lgskern;

 public:
  sutra_sensors(carma_context *context, sutra_telescope *d_tel,
                vector<string> type, int nwfs, long *nxsub, long *nvalid,
                long *npupils, long *npix, long *nphase, long *nrebin,
                long *nfft, long *ntot, long *npup, float *pdiam, float *nphot,
                float *nphot4imat, int *lgs, bool *fakecam, int *maxFluxPerPix,
                int *maxPixValue, int device, bool roket);
  ~sutra_sensors();

  int allocate_buffers();
  int define_mpi_rank(int rank, int size);
  int set_noise(int nwfs, float noise, long seed);

  int initgs(float *xpos, float *ypos, float *lambda, float *mag, float zerop,
             long *size, float *noise, long *seed, float *G, float *thetaML,
             float *dx, float *dy);
  int initgs(float *xpos, float *ypos, float *lambda, float *mag, float zerop,
             long *size, float *noise, float *G, float *thetaML, float *dx,
             float *dy);
  int initgs(float *xpos, float *ypos, float *lambda, float *mag, float zerop,
             long *size, float *G, float *thetaML, float *dx, float *dy);
};

// General utilities

#endif  // _SUTRA_SENSORS_H_
