#ifndef _SUTRA_RTC_H_
#define _SUTRA_RTC_H_

#include <sutra_controller_geo.h>
#include <sutra_controller_ls.h>
#include <sutra_controller_mv.h>
#include <sutra_controller_cured.h>
#include <sutra_controller_kalman.h>
#include <sutra_centroider_bpcog.h>
#include <sutra_centroider_cog.h>
#include <sutra_centroider_corr.h>
#include <sutra_centroider_pyr.h>
#include <sutra_centroider_roof.h>
#include <sutra_centroider_tcog.h>
#include <sutra_centroider_wcog.h>

using namespace std;

class sutra_rtc {
public:
  int device;
  int delay;

  vector<sutra_centroider *> d_centro;
  vector<sutra_controller *> d_control;

  carma_context *current_context;

public:
  sutra_rtc(carma_context *context);
  sutra_rtc(const sutra_rtc& yrtc);
  ~sutra_rtc();

  int
    add_centroider(sutra_sensors *sensors, int nwfs, long nvalid, float offset, float scale, long device,
		 char *typec);
  int
    rm_centroider();
  int
  add_controller_geo(int nactu, int Nphi, long delay, long device);
  int
    add_controller(int nactu, long delay, long device, const char *typec);
  int
    rm_controller();

  int
    do_imat(int ncntrl, sutra_dms *ydm);
  int
    do_imatkl(int ncntrl, sutra_dms *ydms);
  int
    do_imatkl4pzt(int ncntrl, sutra_dms *ydms);
  int
    do_imat_geom(int ncntrl, sutra_dms *ydm, int type);

  int
    do_centroids();
  int
    do_centroids(int ncntrl);
  int
    do_centroids(int ncntrl, bool imat);
  int
    do_centroids_geom(int ncntrl);
  int
    do_control(int ncntrl);
  int
  	apply_control(int ncntrl, sutra_dms *ydm);

};

#endif // _SUTRA_RTC_H_
