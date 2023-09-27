// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_sensors.h
//! \ingroup   libsutra
//! \class     SutraSensors
//! \brief     this class provides the sensors features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

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

class SutraSensors {
 public:
  int device;
  bool roket;
  CarmaContext *current_context;
  size_t nsensors() { return d_wfs.size(); }
  vector<SutraWfs *> d_wfs;
  map<vector<int>, cufftHandle *> campli_plans;
  map<vector<int>, cufftHandle *> fttotim_plans;
  map<vector<int>, cufftHandle *> ftlgskern_plans;
  map<vector<int>, cufftHandle *> field_stop_plans;

  CarmaObj<cuFloatComplex> *d_camplipup;
  CarmaObj<cuFloatComplex> *d_camplifoc;
  CarmaObj<cuFloatComplex> *d_fttotim;
  CarmaObj<cuFloatComplex> *d_ftlgskern;
  CarmaObj<float> *d_lgskern;

 public:
  SutraSensors(CarmaContext *context, SutraTelescope *d_tel,
                vector<string> type, int nwfs, long *nxsub, long *nvalid,
                long *npupils, long *npix, long *nphase, long *nrebin,
                long *nfft, long *ntot, long *npup, float *pdiam, float *nphot,
                float *nphot4imat, int *lgs, bool *fakecam, int *max_flux_per_pix,
                int *max_pix_value, int device, bool roket);
  ~SutraSensors();

  int allocate_buffers();
  int define_mpi_rank(int rank, int size);
  int set_noise(int nwfs, float noise, long seed);
  int set_field_stop(int nwfs, float* field_stop, int N);

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
