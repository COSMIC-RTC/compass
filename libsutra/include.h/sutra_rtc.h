// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_rtc.h
//! \ingroup   libsutra
//! \class     SutraRtc
//! \brief     this class provides the rtc features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

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
#include <sutra_controller_generic_linear.h>
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


/**
 * @brief Add a SutraController object in the RTC
 *
 * @param context       : CarmaContext: carma context
 * @param typec         : string      : Controller type
 * @param device        : long        : GPU device index
 * @param delay         : float       : Loop delay [frames]
 * @param nslope        : int         : Number of slopes
 * @param nactu         : int         : Number of actuators to command
 * @param nslope_buffers: int         : (optional) Number of historic slopes vectors to use
 * @param nstates       : int         : (optional) Number of states in state vector
 * @param nstate_buffers: int         : (optional) Number of historic state vectors to use
 * @param nmodes        : int         : (optional) Number of modes in mode vector
 * @param niir_in       : int         : (optional) Number of input mode vectors for iir filter
 * @param niir_out      : int         : (optional) Number of output mode vectors for iir filter
 * @param polc          : bool        : (optional) Activate the Pseudo Open Loop Control if available
 * @param is_modal      : bool        : (optional) Activate projection from modes to actu if available
 * @param dms           : SutraDms*   : (optional) SutraDms object
 * @param idx_dms       : int*        : (optional) index of DM in SutraDms to command
 * @param ndm           : int         : (optional) Number of DM to command
 * @param idx_centro    : int*        : (optional) Index of centoiders in sutra_rtc.d_centro to handle
 * @param ncentro       : int         : (optional) Number of centroiders handled
 * @param Nphi          : int         : (optional) Number of pixels in the pupil
 * @param wfs_direction : bool        : (optional) Flag for ROKET
 * @return int          : exit status
 */
  int add_controller(CarmaContext *context, std::string typec,long device, float delay, int nslope,
                    int nactu, int nslope_buffers = 0, int nstates = 0, int nstate_buffers = 0,
                    int nmodes = 0, int niir_in = 0, int niir_out = 0, bool polc = false,
                    bool is_modal = false, SutraDms *dms = nullptr, int *idx_dms = nullptr,
                    int ndm = 0, int *idx_centro = nullptr, int ncentro = 0, int Nphi = 0,
                    bool wfs_direction = false);



  int remove_centroider(int ncentro);
  int remove_controller(int ncontrol);

  int do_imat(int ncntrl, SutraDms *ydms, int kernconv);

  int do_imat_basis(int ncntrl, SutraDms *ydm, int nModes, T *m2v,
                    T *pushAmpl, int kernconv);

  int do_imat_geom(int ncntrl, SutraDms *ydm, int type);

  int comp_images_imat(SutraDms *ydm, int kernconv);

  int do_calibrate_img();
  int do_calibrate_img(int ncntrl);
  int do_centroids();
  int do_centroids(int ncntrl);
  int do_centroids(int ncntrl, bool noise);
  int do_centroids_geom(int ncntrl, int type = 0);
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
  do_imat_impl(int ncntrl, SutraDms *ydm, int kernconv, std::true_type);
  int do_imat_impl(int ncntrl, SutraDms *ydm, int kernconv, std::false_type);

  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_imat_basis_impl(int ncntrl, SutraDms *ydm, int nModes, T *m2v,
                     T *pushAmpl, int kernconv, std::true_type);
  int do_imat_basis_impl(int ncntrl, SutraDms *ydm, int nModes, T *m2v,
                         T *pushAmpl, int kernconv, std::false_type);
  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_imat_geom_impl(int ncntrl, SutraDms *ydm, int type, std::true_type);
  int do_imat_geom_impl(int ncntrl, SutraDms *ydm, int type, std::false_type);

  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_centroids_geom_impl(int ncntrl, int type, std::true_type);
  int do_centroids_geom_impl(int ncntrl, int type, std::false_type);

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
                      std::string typec,long device, float delay,  int nslope, int nactu,
                      int nslope_buffers, int nstates, int nstate_buffers, int nmodes,
                      int niir_in, int niir_out, bool polc,bool is_modal,
                      SutraDms *dms, int *idx_dms, int ndm, int *idx_centro, int ncentro,
                      int Nphi, bool wfs_direction,
                      std::false_type);



  int add_controller_impl(CarmaContext *context,
                          vector<SutraController<T, Tout> *> &d_control,
                          std::string typec,long device, float delay,  int nslope, int nactu,
                          int nslope_buffers, int nstates, int nstate_buffers, int nmodes,
                          int niir_in, int niir_out, bool polc,bool is_modal,
                          SutraDms *dms, int *idx_dms, int ndm, int *idx_centro, int ncentro,
                          int Nphi, bool wfs_direction,
                          std::true_type);
};

#endif  // _SUTRA_RTC_H_
