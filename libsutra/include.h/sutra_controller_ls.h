// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_ls.h
//! \ingroup   libsutra
//! \class     sutra_controller_ls
//! \brief     this class provides the controller_ls features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#ifndef _SUTRA_CONTROLLER_LS_H_
#define _SUTRA_CONTROLLER_LS_H_

#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class sutra_controller_ls : public SutraController<Tcomp, Tout> {
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
  int is_modopti;                // Flag for using modal optimization
  int nrec;                      // Number of recorded open slopes measurements
  int nmodes;                    // Number of modes
  Tcomp gmin;                    // Gain min
  Tcomp gmax;                    // Gain max
  int ngain;                     // Number of tested gains between gmin and gmax
  Tcomp Fs;                      // Sampling frequency
  int cpt_rec;                   // Counter for modal gains refresh
  CarmaObj<Tcomp> *d_M2V;       // Modes to Volts matrix
  CarmaObj<Tcomp> *d_S2M;       // Slopes to Modes matrix
  CarmaObj<Tcomp> *d_slpol;     // Open-loop measurements buffer, recorded and
                                 // loaded from Yorick
  CarmaObj<Tcomp> *d_Hcor;      // Transfer function
  CarmaObj<Tcomp> *d_compbuff;  // Buffer for POLC computation
  CarmaObj<Tcomp> *d_compbuff2;  // Buffer for POLC computation

 public:
  sutra_controller_ls(CarmaContext *context, long nslope,
                      long nactu, float delay, SutraDms *dms, int *idx_dms,
                      int ndm, int *idx_centro, int ncentro);
  sutra_controller_ls(const sutra_controller_ls &controller);
  ~sutra_controller_ls();

  string get_type();

  int svdec_imat();
  int build_cmat(int nfilt, bool filt_tt);
  int build_cmat(int nfilt);
  int build_cmat_modopti();
  int frame_delay();
  int comp_com();
  int set_modal_gains(Tcomp *mgain);
  int set_cmat(Tcomp *cmat);
  int set_imat(Tcomp *imat);
  int init_modalOpti(int nmodes, int nrec, Tcomp *M2V, Tcomp gmin, Tcomp gmax,
                     int ngain, Tcomp Fs);
  int loadopen_loopSlp(Tcomp *ol_slopes);
  int modalControlOptimization();
  int compute_Hcor();
};

#endif  // _SUTRA_CONTROLLER_LS_H_
