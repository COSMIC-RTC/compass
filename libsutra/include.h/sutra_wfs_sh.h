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

//! \file      sutra_wfs_sh.h
//! \ingroup   libsutra
//! \class     SutraWfsSH
//! \brief     this class provides the wfs_sh features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_WFS_SH_H_
#define _SUTRA_WFS_SH_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>
#include <map>
#include <vector>

class SutraWfsSH : public SutraWfs {
 public:
  // sh only
  CarmaObj<int> *d_binmap;
  CarmaObj<int> *d_validpuppixx;  // nxsub
  CarmaObj<int> *d_validpuppixy;  // nxsub

 public:
  SutraWfsSH(CarmaContext *context, SutraTelescope *d_tel,
               CarmaObj<cuFloatComplex> *d_camplipup,
               CarmaObj<cuFloatComplex> *d_camplifoc,
               CarmaObj<cuFloatComplex> *d_fttotim, long nxsub, long nvalid,
               long npix, long nphase, long nrebin, long nfft, long ntot,
               long npup, float pdiam, float nphotons, float nphot4imat,
               int lgs, bool fakecam, int max_flux_per_pix, int max_pix_value,
               bool is_low_order, bool roket, int device);
  SutraWfsSH(const SutraWfsSH &wfs);
  ~SutraWfsSH();

  int define_mpi_rank(int rank, int size);
  int allocate_buffers(map<vector<int>, cufftHandle *> campli_plans,
                       map<vector<int>, cufftHandle *> fttotim_plans);

  int load_arrays(int *phasemap, int *hrmap, int *binmap, float *offsets,
                 float *fluxPerSub, int *validsubsx, int *validsubsy,
                 int *istart, int *jstart, cuFloatComplex *kernel);

  int fill_binimage(int async);
  int comp_image(bool noise = true);
  int comp_nphot(float ittime, float optthroughput, float diam, int nxsub,
                 float zerop = 0, float gsmag = 0, float lgsreturnperwatt = 0,
                 float laserpower = 0);
  int set_bincube(float *bincube, int nElem);

 private:
  int comp_generic();
};

#endif  // _SUTRA_WFS_SH_H_
