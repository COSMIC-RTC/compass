// kp_timer.cpp

#include "kp_timer.h"
#include <cstdlib>
#include <iostream>
using namespace std;

kp_timer::kp_timer() {
  is_start = 0;
  a = 0;
}
//
void kp_timer::start() {
  if (is_start) {
    cerr << "Error | kp_timer::start | is already have been started" << endl;
    exit(EXIT_FAILURE);
  }
  gettimeofday(&tv, NULL);
  is_start = true;
}
//
void kp_timer::pause() {
  if (!is_start) {
    cerr << "Error | kp_timer::pause | is already have been stoped" << endl;
    exit(EXIT_FAILURE);
  }
  struct timeval tv2;
  gettimeofday(&tv2, NULL);
  a += (double)(tv2.tv_sec - tv.tv_sec) + 1e-6 * (tv2.tv_usec - tv.tv_usec);
  is_start = false;
}
//
double kp_timer::rez() {
  if (is_start) {
    cerr << "Error | kp_timer::rez | it should be stopped" << endl;
    exit(EXIT_FAILURE);
  }
  return a;
}
//
