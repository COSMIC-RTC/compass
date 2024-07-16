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

//! \file      sutra_wfs_geom.hpp
//! \ingroup   libsutra
//! \class     SutraWfsGeom
//! \brief     this class provides the wfs_geom features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_WFS_GEOM_H_
#define _SUTRA_WFS_GEOM_H_

#include <sutra_lgs.hpp>
#include <sutra_phase.hpp>
#include <sutra_target.hpp>
#include <sutra_telemetry.hpp>
#include <sutra_telescope.hpp>
#include <sutra_wfs.hpp>
#include <map>
#include <vector>

using std::string;
class SutraWfsGeom : public SutraWfs {
 public:
 public:
  SutraWfsGeom(CarmaContext *context, SutraTelescope *d_tel, int64_t nxsub,
                 int64_t nvalid, int64_t nphase, int64_t npup, float pdiam, int32_t device);
  SutraWfsGeom(const SutraWfsGeom &wfs);
  ~SutraWfsGeom();

  int32_t wfs_initarrays(int32_t *phasemap, float *offsets, float *fluxPerSub,
                     int32_t *validsubsx, int32_t *validsubsy);
  int32_t slopes_geom(int32_t type, float *slopes);
  int32_t slopes_geom(int32_t type);

  int32_t define_mpi_rank(int32_t rank, int32_t size) { return EXIT_SUCCESS; }
  int32_t allocate_buffers(map<vector<int32_t>, cufftHandle *> campli_plans,
                       map<vector<int32_t>, cufftHandle *> fttotim_plans) {
    return EXIT_SUCCESS;
  }

  int32_t fill_binimage(int32_t async) { return 0; }
  int32_t comp_image() { return 0; }

 private:
  int32_t comp_generic() { return 0; }
};

#endif  // _SUTRA_WFS_GEOM_H_
