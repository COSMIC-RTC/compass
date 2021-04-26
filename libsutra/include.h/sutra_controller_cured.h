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

//! \file      sutra_controller_cured.h
//! \ingroup   libsutra
//! \class     SutraControllerCured
//! \brief     this class provides the controller_cured features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_CONTROLLER_CURED_H_
#define _SUTRA_CONTROLLER_CURED_H_

#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class SutraControllerCured : public SutraController<Tcomp, Tout> {
 public:
  int ndivs;     // number of subdivision levels for cured
  bool tt_flag;  // flag for separate tt

  // data for CuReD */
  CarmaHostObj<Tcomp> *h_centroids;
  CarmaHostObj<Tcomp> *h_err;
  CarmaObj<Tcomp> *d_err;      // current error
  CarmaObj<Tcomp> *d_cenbuff;  // centroids circular buffer

  // data for CuReD */
  CarmaObj<Tcomp> *d_imat;

  // structures needed to run CuReD */
  // sysCure* h_syscure;
  void *h_syscure;
  // parCure* h_parcure;
  void *h_parcure;

 public:
  SutraControllerCured(CarmaContext *context, long nvalid, long nslope,
                         long nactu, float delay, SutraDms *dms, int *idx_dms,
                         int ndm, int *idx_centro, int ncentro);
  SutraControllerCured(const SutraControllerCured &controller);
  ~SutraControllerCured();

  string get_type() { return "cured"; }

  int comp_com();

  int init_cured(int nxsubs, int *isvalid, int ndivs, int tt);
  int frame_delay();
};

#endif  // _SUTRA_CONTROLLER_CURED_H_
