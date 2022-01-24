// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      sutra_controller_mv.h
//! \ingroup   libsutra
//! \class     sutra_controller_mv
//! \brief     this class provides the controller_mv features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.1
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

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
class sutra_controller_mv : public SutraController<Tcomp, Tout> {
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
  sutra_controller_mv(CarmaContext *context, long nslope,
                      long nactu, float delay, SutraDms *dms, int *idx_dms,
                      int ndm, int *idx_centro, int ncentro);
  sutra_controller_mv(const sutra_controller_mv &controller);
  ~sutra_controller_mv();

  string get_type();

  int svdec_imat();
  int build_cmat(const char *dmtype, char *method);
  int build_cmat(Tcomp cond);
  int frame_delay();
  int comp_com();
  int set_modal_gains(Tcomp *mgain);
  int set_cmat(Tcomp *cmat);
  int set_imat(Tcomp *imat);
  // Florian features
  int load_noisemat(Tcomp *noise);
  int do_covmat(SutraDm *ydm, char *method, int *indx_pup, long dim,
                Tcomp *xpos, Tcomp *ypos, long Nkl, Tcomp norm, Tcomp ampli);
  int do_geomat(CarmaObj<Tcomp> *d_geocov, CarmaObj<Tcomp> *d_IF, long n_pts,
                Tcomp ampli);
  int piston_filt(CarmaObj<Tcomp> *d_statcov);
  int piston_filt_cphim(CarmaObj<Tcomp> *d_cphim, Tcomp *F);
  int filter_cphim(Tcomp *F, Tcomp *Nact);
  int filter_cmat(Tcomp cond);
  int invgen(CarmaObj<Tcomp> *d_mat, Tcomp cond, int job);
  int invgen(CarmaObj<Tcomp> *d_mat, CarmaHostObj<Tcomp> *h_eigen,
             Tcomp cond);
  int invgen_cpu(CarmaObj<Tcomp> *d_mat, CarmaHostObj<Tcomp> *h_eigen,
                 Tcomp cond);
  // int
  // do_statmat(T *statcov,long dim, T *xpos, T *ypos, T norm,
  // CarmaDevice *device);
  int DDiago(CarmaObj<Tcomp> *d_statcov, CarmaObj<Tcomp> *d_geocov);
  int load_covmat(Tcomp *covmat);
  int load_klbasis(Tcomp *klbasis);
  int compute_Cmm(SutraAtmos *atmos, SutraSensors *sensors, double *L0,
                  double *cn2, double *alphaX, double *alphaY, double diamTel,
                  double cobs);
  int compute_Cphim(SutraAtmos *atmos, SutraSensors *sensors, SutraDms *dms,
                    double *L0, double *cn2, double *alphaX, double *alphaY,
                    double *X, double *Y, double *xactu, double *yactu,
                    double diamTel, double *k2, long *NlayerDm,
                    long *indLayerDm, double FoV, double *pitch,
                    double *alt_dm);
};

#endif  // _SUTRA_CONTROLLER_MV_H_
