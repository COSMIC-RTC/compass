#ifndef _SUTRA_ATMOS_H_
#define _SUTRA_ATMOS_H_

#include <sutra_tscreen.h>
#include <map>

using std::map;
using std::pair;

class sutra_atmos {
 public:
  int nscreens;
  map<float, sutra_tscreen *> d_screens;
  float r0;
  carma_context *current_context;

 public:
  sutra_atmos(carma_context *context, int nscreens, float global_r0, float *r0,
              long *size, long *size2, float *altitude, float *windspeed,
              float *winddir, float *deltax, float *deltay, int device);
  ~sutra_atmos();

  int init_screen(float alt, float *h_A, float *h_B, unsigned int *h_istencilx,
                  unsigned int *h_istencily, int seed);

  int add_screen(float alt, long size, long stencilSize, float amplitude,
                 float windspeed, float winddir, float deltax, float deltay,
                 int device);
  int del_screen(const float alt);
  int refresh_screen(float altitude);

  int move_atmos();
  int set_global_r0(float r0);
  int set_seed(float alt, float seed);
};

#endif  // _SUTRA_ATMOS_H_
