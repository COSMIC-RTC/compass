/**
 * \file sutra_atmos.h
 *
 * \class sutra_atmos
 *
 * \ingroup libsutra
 *
 * \brief this class provides the atmos features to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 4.2.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _SUTRA_ATMOS_H_
#define _SUTRA_ATMOS_H_

#include <sutra_tscreen.h>

using std::pair;
using std::vector;

class sutra_atmos {
 public:
  int nscreens;
  vector<sutra_tscreen *> d_screens;
  float r0;
  carma_context *current_context;

 public:
  sutra_atmos(carma_context *context, int nscreens, float global_r0, float *r0,
              long *size, long *size2, float *altitude, float *windspeed,
              float *winddir, float *deltax, float *deltay, int device);
  ~sutra_atmos();

  int init_screen(int idx, float *h_A, float *h_B, unsigned int *h_istencilx,
                  unsigned int *h_istencily, int seed);

  int add_screen(float altitude, long size, long stencilSize, float amplitude,
                 float windspeed, float winddir, float deltax, float deltay,
                 int device);
  int del_screen(const int idx);
  int refresh_screen(int idx);

  int move_atmos();
  int set_global_r0(float r0);
  int set_frac(float *frac);
  int set_seed(int idx, float seed);
};

#endif  // _SUTRA_ATMOS_H_
