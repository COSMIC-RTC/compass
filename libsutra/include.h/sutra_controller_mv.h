// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_mv.h
//! \ingroup   libsutra
//! \class     SutraControllerMv
//! \brief     this class provides the controller_mv features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CONTROLLER_MV_H_
#define _SUTRA_CONTROLLER_MV_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_centroider.h>
#include <sutra_controller.h>
#include <sutra_dm.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>

template <typename Tcomp, typename Tout>
class SutraControllerMv : public SutraController<Tcomp, Tout> {
 public:
  CarmaObj<Tcomp> *d_imat;
  CarmaObj<Tcomp> *d_cmat;
  CarmaObj<Tcomp> *d_gain;

  // Cphim & Cmm features
  CarmaObj<Tcomp> *d_covmat;
  CarmaObj<Tcomp> *d_KLbasis;
  CarmaObj<Tcomp> *d_noisemat;
  CarmaObj<Tcomp> *d_Cmm;
  CarmaObj<Tcomp> *d_Cphim;
  // svd computations
  CarmaHostObj<Tcomp> *h_Cmmeigenvals;
  CarmaHostObj<Tcomp> *h_eigenvals;
  // CarmaObj<Tcomp> *d_U;

  // loop components
  CarmaObj<Tcomp> *d_cenbuff;   // centroids circular buffer
  CarmaObj<Tcomp> *d_compbuff;  // Buffer for computations
  CarmaObj<Tcomp> *d_compbuff2;
  CarmaObj<Tcomp> *d_olmeas;  // Open-loop measurements for POLC
  CarmaObj<Tcomp> *d_err;     // current error

  cublasHandle_t cublas_handle;

 public:
  SutraControllerMv(CarmaContext *context, int64_t nslope,
                      int64_t nactu, float delay, SutraDms *dms, int32_t *idx_dms,
                      int32_t ndm, int32_t *idx_centro, int32_t ncentro);
  SutraControllerMv(const SutraControllerMv &controller);
  ~SutraControllerMv();

  string get_type();

  int32_t svdec_imat();
  int32_t build_cmat(const char *dmtype, char *method);
  int32_t build_cmat(Tcomp cond);
  int32_t frame_delay();
  int32_t comp_com();
  int32_t set_modal_gains(Tcomp *mgain);
  int32_t set_cmat(Tcomp *cmat);
  int32_t set_imat(Tcomp *imat);
  // Florian features
  int32_t load_noisemat(Tcomp *noise);
  int32_t do_covmat(SutraDm *ydm, char *method, int32_t *indx_pup, int64_t dim,
                Tcomp *xpos, Tcomp *ypos, int64_t Nkl, Tcomp norm, Tcomp ampli);
  int32_t do_geomat(CarmaObj<Tcomp> *d_geocov, CarmaObj<Tcomp> *d_IF, int64_t n_pts,
                Tcomp ampli);
  int32_t piston_filt(CarmaObj<Tcomp> *d_statcov);
  int32_t piston_filt_cphim(CarmaObj<Tcomp> *d_cphim, Tcomp *F);
  int32_t filter_cphim(Tcomp *F, Tcomp *Nact);
  int32_t filter_cmat(Tcomp cond);
  int32_t invgen(CarmaObj<Tcomp> *d_mat, Tcomp cond, int32_t job);
  int32_t invgen(CarmaObj<Tcomp> *d_mat, CarmaHostObj<Tcomp> *h_eigen,
             Tcomp cond);
  int32_t invgen_cpu(CarmaObj<Tcomp> *d_mat, CarmaHostObj<Tcomp> *h_eigen,
                 Tcomp cond);
  // int32_t
  // do_statmat(T *statcov,int64_t dim, T *xpos, T *ypos, T norm,
  // CarmaDevice *device);
  int32_t DDiago(CarmaObj<Tcomp> *d_statcov, CarmaObj<Tcomp> *d_geocov);
  int32_t load_covmat(Tcomp *covmat);
  int32_t load_klbasis(Tcomp *klbasis);
  int32_t compute_Cmm(SutraAtmos *atmos, SutraSensors *sensors, double *L0,
                  double *cn2, double *alphaX, double *alphaY, double diamTel,
                  double cobs);
  int32_t compute_Cphim(SutraAtmos *atmos, SutraSensors *sensors, SutraDms *dms,
                    double *L0, double *cn2, double *alphaX, double *alphaY,
                    double *X, double *Y, double *xactu, double *yactu,
                    double diamTel, double *k2, int64_t *NlayerDm,
                    int64_t *indLayerDm, double FoV, double *pitch,
                    double *alt_dm);
};

#endif  // _SUTRA_CONTROLLER_MV_H_
