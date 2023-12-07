// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_wfs_geom.h
//! \ingroup   libsutra
//! \class     SutraWfsGeom
//! \brief     this class provides the wfs_geom features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_WFS_GEOM_H_
#define _SUTRA_WFS_GEOM_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_telescope.h>
#include <sutra_wfs.h>
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
