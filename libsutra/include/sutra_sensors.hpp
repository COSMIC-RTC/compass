// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_sensors.hpp
//! \ingroup   libsutra
//! \class     SutraSensors
//! \brief     this class provides the sensors features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_SENSORS_H_
#define _SUTRA_SENSORS_H_

#include <carma_utils.hpp>
#include <sutra_telescope.hpp>
#include <sutra_utils.hpp>
#include <sutra_wfs.hpp>
#include <sutra_wfs_pyr_pyrhr.hpp>
#include <sutra_wfs_sh.hpp>

#include <map>
#include <vector>

using std::map;
using std::string;
using std::vector;

class SutraSensors {
 public:
  int32_t device;
  bool roket;
  CarmaContext *current_context;
  size_t nsensors() { return d_wfs.size(); }
  vector<SutraWfs *> d_wfs;
  map<vector<int32_t>, cufftHandle *> campli_plans;
  map<vector<int32_t>, cufftHandle *> fttotim_plans;
  map<vector<int32_t>, cufftHandle *> ftlgskern_plans;
  map<vector<int32_t>, cufftHandle *> field_stop_plans;

  CarmaObj<cuFloatComplex> *d_camplipup;
  CarmaObj<cuFloatComplex> *d_camplifoc;
  CarmaObj<cuFloatComplex> *d_fttotim;
  CarmaObj<cuFloatComplex> *d_ftlgskern;
  CarmaObj<float> *d_lgskern;

 public:
  SutraSensors(CarmaContext *context, SutraTelescope *d_tel,
                vector<string> type, int32_t nwfs, int64_t *nxsub, int64_t *nvalid,
                int64_t *npupils, int64_t *npix, int64_t *nphase, int64_t *nrebin,
                int64_t *nfft, int64_t *ntot, int64_t *npup, float *pdiam, float *nphot,
                float *nphot4imat, int32_t *lgs, bool *fakecam, int32_t *max_flux_per_pix,
                int32_t *max_pix_value, int32_t device, bool roket);
  ~SutraSensors();

  int32_t allocate_buffers();
  int32_t define_mpi_rank(int32_t rank, int32_t size);
  int32_t set_noise(int32_t nwfs, float noise, int64_t seed);
  int32_t set_field_stop(int32_t nwfs, float* field_stop, int32_t N);

  int32_t initgs(float *xpos, float *ypos, float *lambda, float *mag, float zerop,
             int64_t *size, float *noise, int64_t *seed, float *G, float *thetaML,
             float *dx, float *dy);
  int32_t initgs(float *xpos, float *ypos, float *lambda, float *mag, float zerop,
             int64_t *size, float *noise, float *G, float *thetaML, float *dx,
             float *dy);
  int32_t initgs(float *xpos, float *ypos, float *lambda, float *mag, float zerop,
             int64_t *size, float *G, float *thetaML, float *dx, float *dy);
};

// General utilities

#endif  // _SUTRA_SENSORS_H_
