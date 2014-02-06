#ifndef _SUTRA_RTC_H_
#define _SUTRA_RTC_H_

#include <sutra_controler.h>
#include <sutra_centroider.h>
#include <sutra_centroider_bpcog.h>
#include <sutra_centroider_cog.h>
#include <sutra_centroider_corr.h>
#include <sutra_centroider_pyr.h>
#include <sutra_centroider_tcog.h>
#include <sutra_centroider_wcog.h>

using namespace std;

class sutra_rtc {
 public:
  int                        device;
  int                        delay;

  vector<sutra_centroider *>  d_centro;
  vector<sutra_controler *>   d_control;

  carma_context *current_context;

 public:
  sutra_rtc(carma_context *context);
  sutra_rtc(const sutra_rtc& yrtc);
  ~sutra_rtc();

  int add_centroider(long nwfs, long nvalid,float offset, float scale, long device, char *typec);
  int rm_centroider();

  int add_controler(long nactu, long delay, long device,const char *typec);
  int rm_controler();

  int do_imat(int ncntrl, sutra_sensors *sensors, sutra_dms *ydm);
  int do_imat_geom(int ncntrl, sutra_sensors *sensors, sutra_dms *ydm,int type);

  int do_centroids(sutra_sensors *sensors);
  int do_centroids(int ncntrl, sutra_sensors *sensors);
  int do_centroids(int ncntrl, sutra_sensors *sensors, bool imat);

  int do_control(int ncntrl, sutra_dms *ydm);

};

#endif // _SUTRA_RTC_H_

