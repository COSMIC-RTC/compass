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

//! \file      sutra_rtc.h
//! \ingroup   libsutra
//! \class     SutraRtc
//! \brief     this class provides the rtc features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_RTC_H_
#define _SUTRA_RTC_H_

#include <sutra_centroider_bpcog.h>
#include <sutra_centroider_cog.h>
#include <sutra_centroider_corr.h>
#include <sutra_centroider_maskedPix.h>
#include <sutra_centroider_pyr.h>
#include <sutra_centroider_tcog.h>
#include <sutra_centroider_wcog.h>
#include <sutra_controller_cured.h>
#include <sutra_controller_generic.h>
#include <sutra_controller_geo.h>
//#include <sutra_controller_kalman.h>
#include <sutra_controller_ls.h>
#include <sutra_controller_mv.h>

template <typename Tin, typename T, typename Tout>
class SutraRtc {
 public:
  vector<SutraCentroider<Tin, T> *> d_centro;
  vector<SutraController<T, Tout> *> d_control;

 public:
  SutraRtc();
  ~SutraRtc();
  int add_centroider(CarmaContext *context, long nvalid, float offset,
                     float scale, bool filter_TT, long device,
                     std::string typec);

  int add_centroider(CarmaContext *context, long nvalid, float offset,
                     float scale, bool filter_TT, long device,
                     std::string typec, SutraWfs *wfs);

  int add_controller(CarmaContext *context, int nvalid, int nslope, int nactu,
                     float delay, long device, std::string typec,
                     SutraDms *dms = nullptr, int *idx_dms = nullptr,
                     int ndm = 0, int *idx_centro = nullptr, int ncentro = 0,
                     int Nphi = 0, bool wfs_direction = false, int nstates = 0);

  int remove_centroider(int ncentro);
  int remove_controller(int ncontrol);

  int do_imat(int ncntrl, SutraDms *ydms);

  int do_imat_basis(int ncntrl, SutraDms *ydm, int nModes, T *m2v,
                    T *pushAmpl);

  int do_imat_geom(int ncntrl, SutraDms *ydm, int type);

  int comp_images_imat(SutraDms *ydm);

  int do_calibrate_img();
  int do_calibrate_img(int ncntrl);
  int do_centroids();
  int do_centroids(int ncntrl);
  int do_centroids(int ncntrl, bool noise);
  int do_centroids_geom(int ncntrl);
  int do_centroids_ref(int ncntrl);
  int do_control(int ncntrl);
  int do_clipping(int ncntrl);
  int apply_control(int ncntrl, bool compVoltage = true);
  int comp_voltage(int ncntrl);
  int remove_ref(int ncntrl);
  int set_centroids_ref(float *centroids_ref);

 private:
  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_imat_impl(int ncntrl, SutraDms *ydm, std::true_type);
  int do_imat_impl(int ncntrl, SutraDms *ydm, std::false_type);

  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_imat_basis_impl(int ncntrl, SutraDms *ydm, int nModes, T *m2v,
                     T *pushAmpl, std::true_type);
  int do_imat_basis_impl(int ncntrl, SutraDms *ydm, int nModes, T *m2v,
                         T *pushAmpl, std::false_type);
  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_imat_geom_impl(int ncntrl, SutraDms *ydm, int type, std::true_type);
  int do_imat_geom_impl(int ncntrl, SutraDms *ydm, int type, std::false_type);

  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_centroids_geom_impl(int ncntrl, std::true_type);
  int do_centroids_geom_impl(int ncntrl, std::false_type);

  template <typename Q = T>
  typename std::enable_if<!std::is_same<Q, half>::value, int>::type
  add_centroider_impl(CarmaContext *context, long nvalid, float offset,
                      float scale, bool filter_TT, long device,
                      std::string typec, SutraWfs *wfs, std::false_type);
  int add_centroider_impl(CarmaContext *context, long nvalid, float offset,
                          float scale, bool filter_TT, long device,
                          std::string typec, SutraWfs *wfs, std::true_type);

  template <typename Q = T>
  typename std::enable_if<!std::is_same<Q, half>::value, int>::type
  add_controller_impl(CarmaContext *context,
                      vector<SutraController<T, Tout> *> &d_control,
                      int nvalid, int nslope, int nactu, float delay,
                      long device, std::string typec, SutraDms *dms,
                      int *idx_dms, int ndm, int *idx_centro, int ncentro,
                      int Nphi, bool wfs_direction, int nstates,
                      std::false_type);
  int add_controller_impl(CarmaContext *context,
                          vector<SutraController<T, Tout> *> &d_control,
                          int nvalid, int nslope, int nactu, float delay,
                          long device, std::string typec, SutraDms *dms,
                          int *idx_dms, int ndm, int *idx_centro, int ncentro,
                          int Nphi, bool wfs_direction, int nstates,
                          std::true_type);
};

#endif  // _SUTRA_RTC_H_
