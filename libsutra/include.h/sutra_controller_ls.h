// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_ls.h
//! \ingroup   libsutra
//! \class     SutraControllerLs
//! \brief     this class provides the controller_ls features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CONTROLLER_LS_H_
#define _SUTRA_CONTROLLER_LS_H_

#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class SutraControllerLs : public SutraController<Tcomp, Tout> {
 public:
  CarmaObj<Tcomp> *d_imat;
  CarmaObj<Tcomp> *d_cmat;
  CarmaObj<Tcomp> *d_gain;

  // svd computations
  CarmaObj<Tcomp> *d_eigenvals;
  CarmaHostObj<Tcomp> *h_eigenvals;
  CarmaObj<Tcomp> *d_U;

  // loop components
  CarmaObj<Tcomp> *d_cenbuff;  // centroids circular buffer
  CarmaObj<Tcomp> *d_err;      // current error

  // Modal optimization components
  int32_t is_modopti;                // Flag for using modal optimization
  int32_t nrec;                      // Number of recorded open slopes measurements
  int32_t nmodes;                    // Number of modes
  Tcomp gmin;                    // Gain min
  Tcomp gmax;                    // Gain max
  int32_t ngain;                     // Number of tested gains between gmin and gmax
  Tcomp Fs;                      // Sampling frequency
  int32_t cpt_rec;                   // Counter for modal gains refresh
  CarmaObj<Tcomp> *d_M2V;       // Modes to Volts matrix
  CarmaObj<Tcomp> *d_S2M;       // Slopes to Modes matrix
  CarmaObj<Tcomp> *d_slpol;     // Open-loop measurements buffer, recorded and
                                 // loaded from Yorick
  CarmaObj<Tcomp> *d_Hcor;      // Transfer function
  CarmaObj<Tcomp> *d_compbuff;  // Buffer for POLC computation
  CarmaObj<Tcomp> *d_compbuff2;  // Buffer for POLC computation

 public:
  SutraControllerLs(CarmaContext *context, int64_t nslope,
                      int64_t nactu, float delay, SutraDms *dms, int32_t *idx_dms,
                      int32_t ndm, int32_t *idx_centro, int32_t ncentro);
  ~SutraControllerLs();

  string get_type();

  int32_t svdec_imat();
  int32_t build_cmat(int32_t nfilt, bool filt_tt);
  int32_t build_cmat(int32_t nfilt);
  int32_t build_cmat_modopti();
  int32_t frame_delay();
  int32_t comp_com();
  int32_t set_modal_gains(Tcomp *mgain);
  int32_t set_cmat(Tcomp *cmat);
  int32_t set_imat(Tcomp *imat);
  int32_t init_modalOpti(int32_t nmodes, int32_t nrec, Tcomp *M2V, Tcomp gmin, Tcomp gmax,
                     int32_t ngain, Tcomp Fs);
  int32_t loadopen_loopSlp(Tcomp *ol_slopes);
  int32_t modalControlOptimization();
  int32_t compute_Hcor();
};

#endif  // _SUTRA_CONTROLLER_LS_H_
