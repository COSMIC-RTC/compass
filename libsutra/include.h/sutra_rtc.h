#ifndef _SUTRA_RTC_H_
#define _SUTRA_RTC_H_

#include <sutra_centroider_bpcog.h>
#include <sutra_centroider_cog.h>
#include <sutra_centroider_corr.h>
#include <sutra_centroider_pyr.h>
#include <sutra_centroider_tcog.h>
#include <sutra_centroider_wcog.h>
#include <sutra_controller_cured.h>
#include <sutra_controller_generic.h>
#include <sutra_controller_geo.h>
#include <sutra_controller_kalman.h>
#include <sutra_controller_ls.h>
#include <sutra_controller_mv.h>

class sutra_rtc {
 public:
  vector<sutra_centroider *> d_centro;
  vector<sutra_controller *> d_control;

 public:
  sutra_rtc();
  ~sutra_rtc();

  int add_centroider(carma_context *context, long nvalid, float offset,
                     float scale, long device, char *typec,
                     sutra_wfs *wfs = nullptr);
  int rm_centroider();
  int add_controller_geo(carma_context *context, int nactu, int Nphi,
                         float delay, long device, sutra_dms *dms, int *idx_dms,
                         int ndm, bool wfs_direction);
  int add_controller(carma_context *context, int nslope, int nactu, float delay,
                     long device, char *typec, sutra_dms *dms = nullptr,
                     int *idx_dms = nullptr, int ndm = 0, int Nphi = 0,
                     bool wfs_direction = false);

  int rm_controller();

  int do_imat(int ncntrl, sutra_dms *ydms);
  int do_imatkl(int ncntrl, sutra_dms *ydms);
  int do_imatkl4pzt(int ncntrl, sutra_dms *ydms);
  int do_imat_geom(int ncntrl, sutra_dms *ydm, int type);

  int do_centroids();
  int do_centroids(int ncntrl);
  int do_centroids(int ncntrl, bool noise);
  int do_centroids(int nctrl, float *bincube, int npix, int ntot);
  int do_centroids_geom(int ncntrl);
  int do_centroids_ref(int ncntrl);
  int do_control(int ncntrl);
  int do_clipping(int ncntrl, float min, float max);
  int apply_control(int ncntrl, sutra_dms *ydm);
  int comp_voltage(int ncntrl);
  int remove_ref(int ncntrl);
  int set_centroids_ref(int ncntrl, float *centroids_ref);
  int get_centroids_ref(int ncntrl, float *centroids_ref);
};

#endif  // _SUTRA_RTC_H_
