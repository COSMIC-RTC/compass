//kp_timer.h


#ifndef __SEGER__KP_TIMER_H__
#define __SEGER__KP_TIMER_H__

#include <sys/time.h>

class kp_timer
{
 public:
   kp_timer();
   
   void start();
   void pause();
   double rez();
 private:
   bool is_start;
   struct timeval tv;
   long double a;
};

#endif
