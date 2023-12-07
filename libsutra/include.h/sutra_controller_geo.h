// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_geo.h
//! \ingroup   libsutra
//! \class     SutraControllerGeo
//! \brief     this class provides the controller_geo features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef SUTRA_CONTROLLER_GEO_H_
#define SUTRA_CONTROLLER_GEO_H_

#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class SutraControllerGeo : public SutraController<Tcomp, Tout> {
 public:
  int64_t Nphi;
  int32_t Ntt;

  CarmaObj<Tcomp> *d_gain;
  CarmaObj<Tcomp> *d_proj;
  CarmaObj<double> *d_phi;
  CarmaObj<Tcomp> *d_phif;
  CarmaObj<int32_t> *d_indx_pup;
  CarmaObj<int32_t> *d_indx_mpup;
  CarmaSparseObj<double> *d_IFsparse;
  CarmaObj<Tcomp> *d_geocov;
  CarmaObj<double> *d_compdouble;
  CarmaObj<float> *d_compfloat;
  CarmaObj<Tcomp> *d_TT;
  CarmaObj<Tcomp> *d_geocovTT;
  //  CarmaObj<T> *d_Btt;
  // CarmaObj<T> *d_cenbuff; // centroids circular buffer

 public:
  SutraControllerGeo(CarmaContext *context, int64_t nactu, int64_t Nphi,
                       float delay, SutraDms *dms, int32_t *idx_dms, int32_t ndm,
                       int32_t *idx_centro, int32_t ncentro, bool wfs_direction);
  SutraControllerGeo(const SutraControllerGeo &controller);
  ~SutraControllerGeo();

  string get_type();

  cusparseHandle_t cusparse_handle() {
    return this->current_context->get_cusparse_handle();
  }
  int32_t load_Btt(Tcomp *Btt_pzt, Tcomp *Btt_TT);
  int32_t load_mgain(Tcomp *mgain);
  int32_t comp_dphi(SutraSource *target, bool wfs_direction);
  int32_t comp_com();
  int32_t init_proj(SutraDms *dms, int32_t *indx_dm, Tcomp *unitpervolt,
                int32_t *indx_pup);
  int32_t init_proj_sparse(SutraDms *dms, int32_t *indx_dm, Tcomp *unitpervolt,
                       int32_t *indx_pup, int32_t *indx_mpup, bool roket);
};

#endif /* SUTRA_CONTROLLER_GEO_H_ */
