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
//! \version   5.4.0
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
  SutraWfsGeom(CarmaContext *context, SutraTelescope *d_tel, long nxsub,
                 long nvalid, long nphase, long npup, float pdiam, int device);
  SutraWfsGeom(const SutraWfsGeom &wfs);
  ~SutraWfsGeom();

  int wfs_initarrays(int *phasemap, float *offsets, float *fluxPerSub,
                     int *validsubsx, int *validsubsy);
  int slopes_geom(int type, float *slopes);
  int slopes_geom(int type);

  int define_mpi_rank(int rank, int size) { return EXIT_SUCCESS; }
  int allocate_buffers(map<vector<int>, cufftHandle *> campli_plans,
                       map<vector<int>, cufftHandle *> fttotim_plans) {
    return EXIT_SUCCESS;
  }

  int fill_binimage(int async) { return 0; }
  int comp_image() { return 0; }

 private:
  int comp_generic() { return 0; }
};

#endif  // _SUTRA_WFS_GEOM_H_
