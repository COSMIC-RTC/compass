// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_geo.h
//! \ingroup   libsutra
//! \class     sutra_controller_geo
//! \brief     this class provides the controller_geo features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

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
