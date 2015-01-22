//kp_cu_timer.cpp

#include "kp_cu_timer.h"
#include <iostream>
#include <cstdlib>
using namespace std;

kp_cu_timer::kp_cu_timer()
{
   is_start = false;
   sdkCreateTimer(&timer);
   sdkResetTimer(&timer);
}
kp_cu_timer::~kp_cu_timer()
{
   sdkDeleteTimer(&timer);
}

void kp_cu_timer::start()
{
   if (is_start)
     {
	cerr<<"Error | kp_cu_timer::start | has already been started"<<endl;
	exit(EXIT_FAILURE);
     }
   cudaDeviceSynchronize();
   sdkStartTimer(&timer);   
   is_start = true;
}


void kp_cu_timer::pause()
{
   if (!is_start)
     {
	cerr<<"Error | kp_cu_timer::pause | has already been stoped"<<endl;
	exit(EXIT_FAILURE);
     }
   cudaDeviceSynchronize();
   sdkStopTimer(&timer);
   is_start = false;
}

void kp_cu_timer::reset()
{
   if (is_start)
     {
	cerr<<"Error | kp_cu_timer::reset | timer must be stoped"<<endl;
	exit(EXIT_FAILURE);
     }
   sdkResetTimer(&timer); 
}

float kp_cu_timer::rez()
{
   if (is_start)
     {
	cerr<<"Error | kp_cu_timer::rez | timer must be stopped"<<endl;
	exit(EXIT_FAILURE);
     }

   return sdkGetTimerValue(&timer)/1000;
   
}
