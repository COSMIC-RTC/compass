// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_controller_utils.hpp
//! \ingroup   libsutra
//! \brief     this file provides the controller utilities to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef SUTRA_CONTROLLER_UTILS_H_
#define SUTRA_CONTROLLER_UTILS_H_

#include <sutra_atmos.hpp>
#include <sutra_dm.hpp>
#include <sutra_sensors.hpp>
#include <cstdio>

struct gtomo_struct {
  double DiamTel;   // Telescope diameter [m]
  double obs;       // Central obstruction
  double pasDPHI;   // DPHI step
  int64_t Nw;          // Number of WFS
  int64_t Nx;          // Total number of valid ssp
  int64_t *Nssp;       // Number of nxsub for each wfs
  int64_t *Nsubap;     // Number of valid subap for each wfs
  double *diamPup;  // Nssp
  double *XPup;     // 0
  double *YPup;     // 0
  double *thetaML;  // 0
  double *alphaX;   // wfs positions
  double *alphaY;   // wfs positions
  double *GsAlt;    // LGS altitude (0)

  double lgs_cst;     // 0
  double spot_width;  // 1
  double lgs_depth;   // 10000
  double lgs_alt;     // 90000.
  int32_t nlgs;           // 0

  int64_t *ioff_d;  // Cumulative number of nvalid ssp with ioff[0] = 0
  int64_t *Nssp_d;  // Nssp
  int64_t *Nsubap_d;

  double *alphaX_d;
  double *alphaY_d;
  double *GsAlt_d;
  double *diamPup_d;
  double *thetaML_d;
  double *X_d;
  double *Y_d;
  double *XPup_d;
  double *YPup_d;

  int64_t max_Nl0;
  int64_t *indexL0_d;
  double *L0diff_d;
  double *h_d;
  double *cn2_d;
  double *tabDPHI_d;
  double *u_d;
  double *v_d;
  double *sspSizeL_d;

  cudaStream_t matcov_stream;
};

struct cphim_struct {
  double DiamTel;
  int64_t Ndphi;  // Useless ?
  int64_t Nlayer;
  double pasDPHI;
  double pasDu;  // Useless
  int64_t int_npts;
  int64_t Nw;
  int64_t Nx;
  int64_t Ndm;
  int64_t Nactu;
  int64_t *Nactu_tot;
  int64_t *indLayerDm;
  int64_t *NlayerDM;
  int64_t *Nssp;
  int64_t *Nsubap;  // Number of valid subap for each wfs
  double *diamPup;
  double *XPup;
  double *YPup;
  double *thetaML;
  double *alphaX;
  double *alphaY;
  double *GsAlt;
  double x0;
  double y0;

  double lgs_cst;
  double spot_width;
  double lgs_depth;
  double lgs_alt;
  int32_t nlgs;

  int64_t *ioff_d;
  int64_t *Nssp_d;
  int64_t *Nactu_tot_d;
  int64_t *indLayerDm_d;
  int64_t *NlayerDM_d;
  double *alphaX_d;
  double *alphaY_d;
  double *GsAlt_d;
  double *diamPup_d;
  double *thetaML_d;
  double *X_d;
  double *Y_d;
  double *XPup_d;
  double *YPup_d;

  int64_t max_Nl0;
  int64_t *indexL0_d;
  int64_t *Nsubap_d;  // Number of valid subap for each wfs
  double *L0diff_d;
  double *h_d;
  double *hDm_d;
  double *cn2_d;
  double *k2_d;
  double *tabDPHI_d;
  double *tab_int_x;
  double *tab_int_y;
  double *xact_d;
  double *yact_d;
  double *u_d;
  double *v_d;
  double *dx_d;
  double *sspSizeL_d;
  double FoV;

  cudaStream_t cphim_stream;
};

void process_err(cudaError_t e, const char *str);
void matts_gpu_gb(double *data, int32_t nrows, int32_t ncols, int32_t xoffset, int32_t yoffset,
                  int32_t lda, struct tomo_struct tomo,
                  struct gtomo_struct *tomo_gpu);
void init_tomo_gpu_gb(struct gtomo_struct *tomo_gpu, SutraAtmos *atmos,
                      SutraSensors *sensors, double diamTel, double cobs);
void free_tomo_gpu_gb(struct gtomo_struct *tomo_gpu);
void update_tomo_atm_gpu_gb(struct gtomo_struct *tomo_gpu,
                            SutraSensors *sensors, SutraAtmos *atmos,
                            double *L0, double *cn2, double *alphaX,
                            double *alphaY);
void update_tomo_sys_gpu_gb(struct gtomo_struct *tomo_gpu,
                            SutraSensors *sensors, double *alphaX,
                            double *alphaY);
void update_cphim_atm(struct cphim_struct *cphim_struct, SutraSensors *sensors,
                      SutraAtmos *atmos, double *L0, double *cn2,
                      double *alphaX, double *alphaY);
void update_cphim_sys(struct cphim_struct *cphim_struct, SutraSensors *sensors,
                      double *alphaX, double *alphaY, double *xactu,
                      double *yactu, double *X, double *Y, int64_t *NlayerDm,
                      int64_t *indLayerDm, double *alt_dm, double *pitch,
                      double *k2, double FoV);
void matcov_gpu_3(double *data, int32_t nrows, int32_t ncols, int32_t xoffset, int32_t yoffset,
                  int32_t lda, struct tomo_struct tomo,
                  struct gtomo_struct *tomo_gpu);
void matcov_gpu_4(float *data, int32_t nrows, int32_t ncols, int32_t xoffset, int32_t yoffset,
                  int32_t lda, struct gtomo_struct *tomo_gpu, SutraAtmos *atmos,
                  SutraSensors *sensors, double *alphaX, double *alphaY);
void CPHIM(float *data, int32_t nrows, int32_t ncols, int32_t xoffset, int32_t yoffset, int32_t lda,
           struct cphim_struct *cphim_struct, SutraAtmos *atmos,
           SutraSensors *sensors, double *alphaX, double *alphaY,
           CarmaDevice *device);
void generateXY(struct gtomo_struct *tomo_gpu, SutraSensors *sensors);
void tab_dphi_gpu_gb(double *tab_dphi, struct gtomo_struct *tomo_gpu,
                     int64_t Ndphi, double *L0diff_d, int32_t Nl0, double convert);
void sub_pos_gpu_gb(struct gtomo_struct *tomo_gpu, int64_t Nlayer);
void sub_pos_cphim(struct cphim_struct *cphim_struct, int64_t Nlayer, int64_t Nw,
                   int64_t Nsubap);
void tab_u831J0(double *tab_int_x, double *tab_int_y, int64_t npts);
void cuda_zcen(double *idata, double *odata, int32_t N, CarmaDevice *device);
void cumsum(double *odata, double *idata, int32_t N);
void init_cphim_struct(struct cphim_struct *cphim_struct, SutraAtmos *atmos,
                       SutraSensors *sensors, SutraDms *dms, double diamTel);
void free_cphim_struct(struct cphim_struct *cphim_struct);
void test_DPHI_highpass(double R, double x0, int64_t npts, CarmaDevice *device);

#endif  // SUTRA_CONTROLLER_UTILS_H_
