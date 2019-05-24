/**
 * \file sutra_rtc.h
 *
 * \class sutra_rtc
 *
 * \ingroup libsutra
 *
 * \brief this class provides the rtc features to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 1.0
 *
 * \date 2011/01/28
 *
 */
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
class sutra_rtc {
 public:
  vector<sutra_centroider<Tin, T> *> d_centro;
  vector<sutra_controller<T, Tout> *> d_control;

 public:
  sutra_rtc();
  ~sutra_rtc();
  int add_centroider(carma_context *context, long nvalid, float offset,
                     float scale, bool filter_TT, long device,
                     std::string typec);

  int add_centroider(carma_context *context, long nvalid, float offset,
                     float scale, bool filter_TT, long device,
                     std::string typec, sutra_wfs *wfs);

  int add_controller(carma_context *context, int nvalid, int nslope, int nactu,
                     float delay, long device, std::string typec,
                     sutra_dms *dms = nullptr, int *idx_dms = nullptr,
                     int ndm = 0, int Nphi = 0, bool wfs_direction = false);

  int remove_centroider(int ncentro);
  int remove_controller(int ncontrol);

  int do_imat(int ncntrl, sutra_dms *ydms);

  int do_imat_basis(int ncntrl, sutra_dms *ydm, int nModes, T *m2v,
                    T *pushAmpl);

  int do_imat_geom(int ncntrl, sutra_dms *ydm, int type);

  int comp_images_imat(sutra_dms *ydm);

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
  do_imat_impl(int ncntrl, sutra_dms *ydm, std::true_type);
  int do_imat_impl(int ncntrl, sutra_dms *ydm, std::false_type);

  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_imat_basis_impl(int ncntrl, sutra_dms *ydm, int nModes, T *m2v,
                     T *pushAmpl, std::true_type);
  int do_imat_basis_impl(int ncntrl, sutra_dms *ydm, int nModes, T *m2v,
                         T *pushAmpl, std::false_type);
  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_imat_geom_impl(int ncntrl, sutra_dms *ydm, int type, std::true_type);
  int do_imat_geom_impl(int ncntrl, sutra_dms *ydm, int type, std::false_type);

  template <typename Q = T>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  do_centroids_geom_impl(int ncntrl, std::true_type);
  int do_centroids_geom_impl(int ncntrl, std::false_type);

  template <typename Q = T>
  typename std::enable_if<!std::is_same<Q, half>::value, int>::type
  add_centroider_impl(carma_context *context, long nvalid, float offset,
                      float scale, bool filter_TT, long device,
                      std::string typec, sutra_wfs *wfs, std::false_type);
  int add_centroider_impl(carma_context *context, long nvalid, float offset,
                          float scale, bool filter_TT, long device,
                          std::string typec, sutra_wfs *wfs, std::true_type);

  template <typename Q = T>
  typename std::enable_if<!std::is_same<Q, half>::value, int>::type
  add_controller_impl(carma_context *context,
                      vector<sutra_controller<T, Tout> *> &d_control,
                      int nvalid, int nslope, int nactu, float delay,
                      long device, std::string typec, sutra_dms *dms,
                      int *idx_dms, int ndm, int Nphi, bool wfs_direction,
                      std::false_type);
  int add_controller_impl(carma_context *context,
                          vector<sutra_controller<T, Tout> *> &d_control,
                          int nvalid, int nslope, int nactu, float delay,
                          long device, std::string typec, sutra_dms *dms,
                          int *idx_dms, int ndm, int Nphi, bool wfs_direction,
                          std::true_type);
};

#endif  // _SUTRA_RTC_H_
