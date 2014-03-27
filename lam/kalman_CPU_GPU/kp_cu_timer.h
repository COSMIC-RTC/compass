//kp_cu_timer.h

#ifndef __SEGER__KP_CU_TIMER_H__
#define __SEGER__KP_CU_TIMER_H__
#include <cuda_runtime.h>
#include "helper_timer.h"

class kp_cu_timer
{
 public:
   kp_cu_timer();
   ~kp_cu_timer();
   
   void start();
   void pause();
   void reset();
   float rez();
 private:
   bool is_start;
   StopWatchInterface* timer;
};

#endif



