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

//! \file      sutra_wfs_geom.h
//! \ingroup   libsutra
//! \class     SutraWfsGeom
//! \brief     this class provides the wfs_geom features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
