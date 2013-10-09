#ifndef _YOGA_RTC_H_
#define _YOGA_RTC_H_

#include <yoga_centroider.h>
#include <yoga_controler.h>

using namespace std;

class yoga_rtc {
 public:
  int                        device;
  int                        delay;

  vector<yoga_centroider *>  d_centro;
  vector<yoga_controler *>   d_control;

  yoga_context *current_context;

 public:
  yoga_rtc(yoga_context *context);
  yoga_rtc(const yoga_rtc& yrtc);
  ~yoga_rtc();

  int add_centroider(long nwfs, long nvalid,float offset, float scale, long device, char *typec);
  int rm_centroider();

  int add_controler(long nactu, long delay, long device,const char *typec);
  int rm_controler();

  int do_imat(int ncntrl, yoga_sensors *sensors, yoga_dms *ydm);
  int do_imat_geom(int ncntrl, yoga_sensors *sensors, yoga_dms *ydm,int type);

  int do_centroids(yoga_sensors *sensors);
  int do_centroids(int ncntrl, yoga_sensors *sensors);
  int do_centroids(int ncntrl, yoga_sensors *sensors, bool imat);

  int do_control(int ncntrl, yoga_dms *ydm);

};

#endif // _YOGA_RTC_H_

