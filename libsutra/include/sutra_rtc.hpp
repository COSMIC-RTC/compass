// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_rtc.hpp
//! \ingroup   libsutra
//! \class     SutraRtc
//! \brief     this class provides the rtc features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_RTC_H_
#define _SUTRA_RTC_H_

#include <sutra_centroider_bpcog.hpp>
#include <sutra_centroider_cog.hpp>
#include <sutra_centroider_corr.hpp>
#include <sutra_centroider_maskedPix.hpp>
#include <sutra_centroider_pyr.hpp>
#include <sutra_centroider_tcog.hpp>
#include <sutra_centroider_wcog.hpp>
#include <sutra_controller_generic.hpp>
#include <sutra_controller_generic_linear.hpp>
#include <sutra_controller_geo.hpp>
#include <sutra_controller_ls.hpp>
#include <sutra_controller_mv.hpp>

template <typename Tin, typename Tcomp, typename Tout>
class SutraRtc {
 public:
  vector<SutraCentroider<Tin, Tcomp> *> d_centro;
  vector<SutraController<Tcomp, Tout> *> d_control;

 public:
  SutraRtc();
  ~SutraRtc();
  int32_t add_centroider(CarmaContext *context, int64_t nvalid, float offset,
                     float scale, bool filter_TT, int64_t device,
                     std::string typec);

  int32_t add_centroider(CarmaContext *context, int64_t nvalid, float offset,
                     float scale, bool filter_TT, int64_t device,
                     std::string typec, SutraWfs *wfs);


/**
 * @brief Add a SutraController object in the RTC
 *
 * @param context       : CarmaContext: carma context
 * @param typec         : string      : Controller type
 * @param device        : int64_t        : GPU device index
 * @param delay         : float       : Loop delay [frames]
 * @param nslope        : int32_t         : Number of slopes
 * @param nactu         : int32_t         : Number of actuators to command
 * @param nslope_buffers: int32_t         : (optional) Number of historic slopes vectors to use
 * @param nstates       : int32_t         : (optional) Number of states in state vector
 * @param nstate_buffers: int32_t         : (optional) Number of historic state vectors to use
 * @param nmodes        : int32_t         : (optional) Number of modes in mode vector
 * @param niir_in       : int32_t         : (optional) Number of input mode vectors for iir filter
 * @param niir_out      : int32_t         : (optional) Number of output mode vectors for iir filter
 * @param polc          : bool        : (optional) Activate the Pseudo Open Loop Control if available
 * @param is_modal      : bool        : (optional) Activate projection from modes to actu if available
 * @param dms           : SutraDms*   : (optional) SutraDms object
 * @param idx_dms       : int32_t*        : (optional) index of DM in SutraDms to command
 * @param ndm           : int32_t         : (optional) Number of DM to command
 * @param idx_centro    : int32_t*        : (optional) Index of centoiders in sutra_rtc.d_centro to handle
 * @param ncentro       : int32_t         : (optional) Number of centroiders handled
 * @param Nphi          : int32_t         : (optional) Number of pixels in the pupil
 * @param wfs_direction : bool        : (optional) Flag for ROKET
 * @return int32_t          : exit status
 */
  int32_t add_controller(CarmaContext *context, std::string typec,int64_t device, float delay, int32_t nslope,
                    int32_t nactu, int32_t nslope_buffers = 0, int32_t nstates = 0, int32_t nstate_buffers = 0,
                    int32_t nmodes = 0, int32_t niir_in = 0, int32_t niir_out = 0, bool polc = false,
                    bool is_modal = false, SutraDms *dms = nullptr, int32_t *idx_dms = nullptr,
                    int32_t ndm = 0, int32_t *idx_centro = nullptr, int32_t ncentro = 0, int32_t Nphi = 0,
                    bool wfs_direction = false);



  int32_t remove_centroider(int32_t ncentro);
  int32_t remove_controller(int32_t ncontrol);

  int32_t do_imat(int32_t ncntrl, SutraDms *ydms, int32_t kernconv);

  int32_t do_imat_basis(int32_t ncntrl, SutraDms *ydm, int32_t nModes, Tcomp *m2v,
                    Tcomp *pushAmpl, int32_t kernconv);

  int32_t do_imat_geom(int32_t ncntrl, SutraDms *ydm, int32_t type);

  int32_t comp_images_imat(SutraDms *ydm, int32_t kernconv);

  int32_t do_calibrate_img();
  int32_t do_calibrate_img(int32_t ncntrl);
  int32_t do_centroids();
  int32_t do_centroids(int32_t ncntrl);
  int32_t do_centroids(int32_t ncntrl, bool noise);
  int32_t do_centroids_geom(int32_t ncntrl, int32_t type = 0);
  int32_t do_centroids_ref(int32_t ncntrl);
  int32_t do_control(int32_t ncntrl);
  int32_t do_clipping(int32_t ncntrl);
  int32_t apply_control(int32_t ncntrl, bool compVoltage = true);
  int32_t comp_voltage(int32_t ncntrl);
  int32_t remove_ref(int32_t ncntrl);
  int32_t set_centroids_ref(float *centroids_ref);

 private:
  template <typename Q = Tcomp>
  typename std::enable_if<std::is_same<Q, float>::value, int32_t>::type
  do_imat_impl(int32_t ncntrl, SutraDms *ydm, int32_t kernconv, std::true_type);
  int32_t do_imat_impl(int32_t ncntrl, SutraDms *ydm, int32_t kernconv, std::false_type);

  template <typename Q = Tcomp>
  typename std::enable_if<std::is_same<Q, float>::value, int32_t>::type
  do_imat_basis_impl(int32_t ncntrl, SutraDms *ydm, int32_t nModes, Tcomp *m2v,
                     Tcomp *pushAmpl, int32_t kernconv, std::true_type);
  int32_t do_imat_basis_impl(int32_t ncntrl, SutraDms *ydm, int32_t nModes, Tcomp *m2v,
                         Tcomp *pushAmpl, int32_t kernconv, std::false_type);
  template <typename Q = Tcomp>
  typename std::enable_if<std::is_same<Q, float>::value, int32_t>::type
  do_imat_geom_impl(int32_t ncntrl, SutraDms *ydm, int32_t type, std::true_type);
  int32_t do_imat_geom_impl(int32_t ncntrl, SutraDms *ydm, int32_t type, std::false_type);

  template <typename Q = Tcomp>
  typename std::enable_if<std::is_same<Q, float>::value, int32_t>::type
  do_centroids_geom_impl(int32_t ncntrl, int32_t type, std::true_type);
  int32_t do_centroids_geom_impl(int32_t ncntrl, int32_t type, std::false_type);

  template <typename Q = Tcomp>
  typename std::enable_if<!std::is_same<Q, half>::value, int32_t>::type
  add_centroider_impl(CarmaContext *context, int64_t nvalid, float offset,
                      float scale, bool filter_TT, int64_t device,
                      std::string typec, SutraWfs *wfs, std::false_type);
  int32_t add_centroider_impl(CarmaContext *context, int64_t nvalid, float offset,
                          float scale, bool filter_TT, int64_t device,
                          std::string typec, SutraWfs *wfs, std::true_type);

  template <typename Q = Tcomp>
  typename std::enable_if<!std::is_same<Q, half>::value, int32_t>::type
  add_controller_impl(CarmaContext *context,
                      vector<SutraController<Tcomp, Tout> *> &d_control,
                      std::string typec,int64_t device, float delay,  int32_t nslope, int32_t nactu,
                      int32_t nslope_buffers, int32_t nstates, int32_t nstate_buffers, int32_t nmodes,
                      int32_t niir_in, int32_t niir_out, bool polc,bool is_modal,
                      SutraDms *dms, int32_t *idx_dms, int32_t ndm, int32_t *idx_centro, int32_t ncentro,
                      int32_t Nphi, bool wfs_direction,
                      std::false_type);



  int32_t add_controller_impl(CarmaContext *context,
                          vector<SutraController<Tcomp, Tout> *> &d_control,
                          std::string typec,int64_t device, float delay,  int32_t nslope, int32_t nactu,
                          int32_t nslope_buffers, int32_t nstates, int32_t nstate_buffers, int32_t nmodes,
                          int32_t niir_in, int32_t niir_out, bool polc,bool is_modal,
                          SutraDms *dms, int32_t *idx_dms, int32_t ndm, int32_t *idx_centro, int32_t ncentro,
                          int32_t Nphi, bool wfs_direction,
                          std::true_type);
};

#endif  // _SUTRA_RTC_H_
