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

template <typename T>
class sutra_rtc {
 public:
  vector<sutra_centroider<T> *> d_centro;
  vector<sutra_controller<T> *> d_control;

 public:
  sutra_rtc();
  ~sutra_rtc();
  int add_centroider(carma_context *context, long nvalid, T offset, T scale,
                     long device, char *typec);
  int add_centroider(carma_context *context, long nvalid, T offset, T scale,
                     long device, char *typec, sutra_wfs *wfs);
  int rm_centroider();
  int add_controller_geo(carma_context *context, int nactu, int Nphi, T delay,
                         long device, sutra_dms *dms, int *idx_dms, int ndm,
                         bool wfs_direction);
  int add_controller(carma_context *context, int nvalid, int nslope, int nactu,
                     T delay, long device, char *typec,
                     sutra_dms *dms = nullptr, int *idx_dms = nullptr,
                     int ndm = 0, int Nphi = 0, bool wfs_direction = false);

  int rm_controller();

  int do_imat(int ncntrl, sutra_dms *ydms);
  int do_imat_basis(int ncntrl, sutra_dms *ydm, int nModes, T *m2v,
                    T *pushAmpl);
  int do_imat_geom(int ncntrl, sutra_dms *ydm, int type);
  int comp_images_imat(sutra_dms *ydm);

  int do_centroids();
  int do_centroids(int ncntrl);
  int do_centroids(int ncntrl, bool noise);
  int do_centroids_geom(int ncntrl);
  int do_centroids_ref(int ncntrl);
  int do_control(int ncntrl);
  int do_clipping(int ncntrl, T min, T max);
  int apply_control(int ncntrl, sutra_dms *ydm, bool compVoltage = true);
  int comp_voltage(int ncntrl);
  int remove_ref(int ncntrl);
  int set_centroids_ref(float *centroids_ref);
};

#endif  // _SUTRA_RTC_H_
