#ifndef _SUTRA_TURBU_H_
#define _SUTRA_TURBU_H_

#include <vector>
#include <map>
#include <sutra_ao_utils.h>
#include <sutra_phase.h>

using namespace std;

class sutra_tscreen {
public:

  int device; // The device # 
  sutra_phase *d_tscreen; // The phase screen   
  carma_obj<float> *d_tscreen_o; // Additional space of the same size as the phase screen
  carma_obj<float> *d_A; // A matrix for extrusion
  carma_obj<float> *d_B; // B matrix for extrusion
  carma_obj<unsigned int> *d_istencilx; // stencil for column extrusion
  carma_obj<unsigned int> *d_istencily; // stencil for row extrusion
  carma_obj<float> *d_z; // tmp array for extrusion process
  carma_obj<float> *d_noise; // tmp array containing random numbers
  carma_obj<float> *d_ytmp; // contains the extrude update (row or column)
  long screen_size; // size of phase screens
  float amplitude; // amplitude for extrusion (r0^-5/6)
  float altitude;
  float windspeed;
  float winddir;
  float deltax; // number of rows to extrude per iteration
  float deltay; // number of lines to extrude per iteration
  //internal
  float accumx;
  float accumy;
  cudaChannelFormatDesc channelDesc; // Channel descriptor for texture memory access

  carma_obj<cuFloatComplex> *d_tscreen_c; // Additional space for von karman screen generation
  float norm_vk;
  bool vk_on;
  carma_context *current_context;

public:
  sutra_tscreen(carma_context *context, long size, long size2, float amplitude,
      float altitude, float windspeed, float winddir, float deltax,
      float deltay, int device);
  sutra_tscreen(const sutra_tscreen& tscreen);
  ~sutra_tscreen();

  int
  init_screen(float *h_A, float *h_B, unsigned int *h_istencilx,
      unsigned int *h_istencily, int seed);
  int
  extrude(int dir);
  int
  init_vk(int seed, int pupd);
  int
  generate_vk(float l0, int nalias);
};

class sutra_atmos {
public:
  int nscreens;
  map<float, sutra_tscreen *> d_screens;
  float r0;
  carma_obj<float> *d_pupil;
  carma_context *current_context;

public:
  sutra_atmos(carma_context *context, int nscreens, float *r0, long *size,
      long *size2, float *altitude, float *windspeed, float *winddir,
      float *deltax, float *deltay, float *pupil, int device);
  ~sutra_atmos();

  int
  init_screen(float alt, float *h_A, float *h_B, unsigned int *h_istencilx,
      unsigned int *h_istencily, int seed);
  int
  move();
};

int
gene_vonkarman(cuFloatComplex *d_odata, float *d_idata, float k0, int nalias,
    int nx, int ny, int block_size);
int
norm_pscreen(float *d_odata, float *d_idata, int nx, int ny, float norm_fact,
    int device);

extern "C" {
}
;

#endif // _SUTRA_TURBU_H_
