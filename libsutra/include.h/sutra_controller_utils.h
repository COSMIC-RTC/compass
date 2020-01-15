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

//! \file      sutra_controller_utils.h
//! \ingroup   libsutra
//! \brief     this file provides the controller utilities to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef SUTRA_CONTROLLER_UTILS_H_
#define SUTRA_CONTROLLER_UTILS_H_

#include <sutra_atmos.h>
#include <sutra_dm.h>
#include <sutra_sensors.h>
#include <cstdio>

struct gtomo_struct {
  double DiamTel;   // Telescope diameter [m]
  double obs;       // Central obstruction
  double pasDPHI;   // DPHI step
  long Nw;          // Number of WFS
  long Nx;          // Total number of valid ssp
  long *Nssp;       // Number of nxsub for each wfs
  long *Nsubap;     // Number of valid subap for each wfs
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
  int nlgs;           // 0

  long *ioff_d;  // Cumulative number of nvalid ssp with ioff[0] = 0
  long *Nssp_d;  // Nssp
  long *Nsubap_d;

  double *alphaX_d;
  double *alphaY_d;
  double *GsAlt_d;
  double *diamPup_d;
  double *thetaML_d;
  double *X_d;
  double *Y_d;
  double *XPup_d;
  double *YPup_d;

  long max_Nl0;
  long *indexL0_d;
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
  long Ndphi;  // Useless ?
  long Nlayer;
  double pasDPHI;
  double pasDu;  // Useless
  long int_npts;
  long Nw;
  long Nx;
  long Ndm;
  long Nactu;
  long *Nactu_tot;
  long *indLayerDm;
  long *NlayerDM;
  long *Nssp;
  long *Nsubap;  // Number of valid subap for each wfs
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
  int nlgs;

  long *ioff_d;
  long *Nssp_d;
  long *Nactu_tot_d;
  long *indLayerDm_d;
  long *NlayerDM_d;
  double *alphaX_d;
  double *alphaY_d;
  double *GsAlt_d;
  double *diamPup_d;
  double *thetaML_d;
  double *X_d;
  double *Y_d;
  double *XPup_d;
  double *YPup_d;

  long max_Nl0;
  long *indexL0_d;
  long *Nsubap_d;  // Number of valid subap for each wfs
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
void matts_gpu_gb(double *data, int nrows, int ncols, int xoffset, int yoffset,
                  int lda, struct tomo_struct tomo,
                  struct gtomo_struct *tomo_gpu);
void init_tomo_gpu_gb(struct gtomo_struct *tomo_gpu, sutra_atmos *atmos,
                      sutra_sensors *sensors, double diamTel, double cobs);
void free_tomo_gpu_gb(struct gtomo_struct *tomo_gpu);
void update_tomo_atm_gpu_gb(struct gtomo_struct *tomo_gpu,
                            sutra_sensors *sensors, sutra_atmos *atmos,
                            double *L0, double *cn2, double *alphaX,
                            double *alphaY);
void update_tomo_sys_gpu_gb(struct gtomo_struct *tomo_gpu,
                            sutra_sensors *sensors, double *alphaX,
                            double *alphaY);
void update_cphim_atm(struct cphim_struct *cphim_struct, sutra_sensors *sensors,
                      sutra_atmos *atmos, double *L0, double *cn2,
                      double *alphaX, double *alphaY);
void update_cphim_sys(struct cphim_struct *cphim_struct, sutra_sensors *sensors,
                      double *alphaX, double *alphaY, double *xactu,
                      double *yactu, double *X, double *Y, long *NlayerDm,
                      long *indLayerDm, double *alt_dm, double *pitch,
                      double *k2, double FoV);
void matcov_gpu_3(double *data, int nrows, int ncols, int xoffset, int yoffset,
                  int lda, struct tomo_struct tomo,
                  struct gtomo_struct *tomo_gpu);
void matcov_gpu_4(float *data, int nrows, int ncols, int xoffset, int yoffset,
                  int lda, struct gtomo_struct *tomo_gpu, sutra_atmos *atmos,
                  sutra_sensors *sensors, double *alphaX, double *alphaY);
void CPHIM(float *data, int nrows, int ncols, int xoffset, int yoffset, int lda,
           struct cphim_struct *cphim_struct, sutra_atmos *atmos,
           sutra_sensors *sensors, double *alphaX, double *alphaY,
           carma_device *device);
void generateXY(struct gtomo_struct *tomo_gpu, sutra_sensors *sensors);
void tab_dphi_gpu_gb(double *tab_dphi, struct gtomo_struct *tomo_gpu,
                     long Ndphi, double *L0diff_d, int Nl0, double convert);
void sub_pos_gpu_gb(struct gtomo_struct *tomo_gpu, long Nlayer);
void sub_pos_cphim(struct cphim_struct *cphim_struct, long Nlayer, long Nw,
                   long Nsubap);
void tab_u831J0(double *tab_int_x, double *tab_int_y, long npts);
void cuda_zcen(double *idata, double *odata, int N, carma_device *device);
void cumsum(double *odata, double *idata, int N);
void init_cphim_struct(struct cphim_struct *cphim_struct, sutra_atmos *atmos,
                       sutra_sensors *sensors, sutra_dms *dms, double diamTel);
void free_cphim_struct(struct cphim_struct *cphim_struct);
void test_DPHI_highpass(double R, double x0, long npts, carma_device *device);

#endif  // SUTRA_CONTROLLER_UTILS_H_
