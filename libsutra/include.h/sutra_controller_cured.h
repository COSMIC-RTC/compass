// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_cured.h
//! \ingroup   libsutra
//! \class     SutraControllerCured
//! \brief     this class provides the controller_cured features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

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
  SutraControllerCured(CarmaContext *context, long nslope,
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
