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

//! \file      sutra_controller_geo.h
//! \ingroup   libsutra
//! \class     sutra_controller_geo
//! \brief     this class provides the controller_geo features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef SUTRA_CONTROLLER_GEO_H_
#define SUTRA_CONTROLLER_GEO_H_

#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class sutra_controller_geo : public SutraController<Tcomp, Tout> {
 public:
  long Nphi;
  int Ntt;

  CarmaObj<Tcomp> *d_gain;
  CarmaObj<Tcomp> *d_proj;
  CarmaObj<double> *d_phi;
  CarmaObj<Tcomp> *d_phif;
  CarmaObj<int> *d_indx_pup;
  CarmaObj<int> *d_indx_mpup;
  CarmaSparseObj<double> *d_IFsparse;
  CarmaObj<Tcomp> *d_geocov;
  CarmaObj<double> *d_compdouble;
  CarmaObj<float> *d_compfloat;
  CarmaObj<Tcomp> *d_TT;
  CarmaObj<Tcomp> *d_geocovTT;
  //  CarmaObj<T> *d_Btt;
  // CarmaObj<T> *d_cenbuff; // centroids circular buffer

 public:
  sutra_controller_geo(CarmaContext *context, long nactu, long Nphi,
                       float delay, SutraDms *dms, int *idx_dms, int ndm,
                       int *idx_centro, int ncentro, bool wfs_direction);
  sutra_controller_geo(const sutra_controller_geo &controller);
  ~sutra_controller_geo();

  string get_type();

  cusparseHandle_t cusparse_handle() {
    return this->current_context->get_cusparse_handle();
  }
  int load_Btt(Tcomp *Btt_pzt, Tcomp *Btt_TT);
  int load_mgain(Tcomp *mgain);
  int comp_dphi(SutraSource *target, bool wfs_direction);
  int comp_com();
  int init_proj(SutraDms *dms, int *indx_dm, Tcomp *unitpervolt,
                int *indx_pup);
  int init_proj_sparse(SutraDms *dms, int *indx_dm, Tcomp *unitpervolt,
                       int *indx_pup, int *indx_mpup, bool roket);
};

#endif /* SUTRA_CONTROLLER_GEO_H_ */
