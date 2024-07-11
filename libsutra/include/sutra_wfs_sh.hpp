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

//! \file      sutra_wfs_sh.hpp
//! \ingroup   libsutra
//! \class     SutraWfsSH
//! \brief     this class provides the wfs_sh features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_WFS_SH_H_
#define _SUTRA_WFS_SH_H_

#include <sutra_lgs.hpp>
#include <sutra_phase.hpp>
#include <sutra_target.hpp>
#include <sutra_telemetry.hpp>
#include <sutra_utils.hpp>
#include <sutra_wfs.hpp>
#include <map>
#include <vector>

class SutraWfsSH : public SutraWfs {
 public:
  // sh only
  CarmaObj<int32_t> *d_binmap;
  CarmaObj<int32_t> *d_validpuppixx;  // nxsub
  CarmaObj<int32_t> *d_validpuppixy;  // nxsub
  CarmaObj<cuFloatComplex> *d_fsamplipup; // Field stop computation arrays
  CarmaObj<cuFloatComplex> *d_fsamplifoc; // Field stop computation arrays
  cufftHandle *fsampli_plan;

 public:
  SutraWfsSH(CarmaContext *context, SutraTelescope *d_tel,
               CarmaObj<cuFloatComplex> *d_camplipup,
               CarmaObj<cuFloatComplex> *d_camplifoc,
               CarmaObj<cuFloatComplex> *d_fttotim, int64_t nxsub, int64_t nvalid,
               int64_t npix, int64_t nphase, int64_t nrebin, int64_t nfft, int64_t ntot,
               int64_t npup, float pdiam, float nphotons, float nphot4imat,
               int32_t lgs, bool fakecam, int32_t max_flux_per_pix, int32_t max_pix_value,
               bool is_low_order, bool roket, int32_t device);
  SutraWfsSH(const SutraWfsSH &wfs);
  ~SutraWfsSH();

  int32_t define_mpi_rank(int32_t rank, int32_t size);
  int32_t allocate_buffers(map<vector<int32_t>, cufftHandle *> campli_plans,
                       map<vector<int32_t>, cufftHandle *> fttotim_plans);

  int32_t load_arrays(int32_t *phasemap, int32_t *hrmap, int32_t *binmap, float *offsets,
                 float *fluxPerSub, int32_t *validsubsx, int32_t *validsubsy,
                 int32_t *istart, int32_t *jstart, float *ttprojmat,
                cuFloatComplex *kernel);

  int32_t fill_binimage(int32_t async);
  int32_t comp_image(bool noise = true);
  int32_t comp_nphot(float ittime, float optthroughput, float diam, int32_t nxsub,
                 float zerop = 0, float gsmag = 0, float lgsreturnperwatt = 0,
                 float laserpower = 0);
  int32_t set_bincube(float *bincube, int32_t nElem);
  int32_t set_field_stop(map<vector<int32_t>, cufftHandle *> campli_plans, float* field_stop, int32_t N);

 private:
  int32_t comp_generic();
};

#endif  // _SUTRA_WFS_SH_H_
