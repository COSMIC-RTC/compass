#ifndef _SUTRA_TARGET_H_
#define _SUTRA_TARGET_H_

#include <map>
#include <string>
#include <vector>
#include <carma.h>
#include <sutra_turbu.h>
#include <sutra_lgs.h>
#include <sutra_dm.h>
#include <sutra_telemetry.h>

using namespace std;
using std::string;

class sutra_source {
public:

  int device; // device # 
  float tposx; // x position of target on the sky  
  float tposy; // y position of target on the sky
  long npos; // number of points in the pupil
  float mag; // brightness of target
  float lambda; // imaging lambda
  float zp; // imaging zero point
  float scale; // phase scale
  bool lgs; // flag for lgs
  string type; // type of source : target / wfs
  int block_size; // optimum block size of device

  float strehl_se; // short exposure strehl
  float strehl_le; // long exposure strehl
  float ref_strehl; // reference for strehl computation
  int strehl_counter; // counter for le strehl computation
  float phase_var; // current phase variance in the pupil
  float phase_var_avg; // average phase variance in the pupil
  int phase_var_count; // counter for average phase variance in the pupil

  sutra_phase *d_phase; // phase for this target
  carma_host_obj<float> *phase_telemetry; // 
  sutra_lgs *d_lgs; // the lgs object
  carma_obj<float> *object; // the object intensity map
  carma_obj<float> *d_pupil; // the pupil mask
  carma_obj<float> *d_image; // the resulting image for this target
  carma_obj<float> *d_leimage; // the long exposure image for this target
  carma_obj<cuFloatComplex> *d_amplipup; // the complex amplitude in the pupil plane
  carma_obj<float> *d_phasepts; // the valid phase points in the pupil (target only)
  carma_obj<int> *d_wherephase; // positions of valid phase points in the pupil (target only)
  map<type_screen, float> xoff; // x reference for raytracing
  map<type_screen, float> yoff; // y reference for raytracing
  carma_context *current_context;

public:
  sutra_source(carma_context *context, float xpos, float ypos, float lambda,
      float mag, long size, string type, float *pupil, int device);
  sutra_source(carma_context *context, float xpos, float ypos, float lambda,
      float mag, long size, string type, int device);
  sutra_source(const sutra_source& source);
  ~sutra_source();
  inline int
  init_source(carma_context *context, float xpos, float ypos, float lambda,
      float mag, long size, string type, int device);
  int
  add_layer(string type, float alt, float xoff, float yoff);
  int
  remove_layer(string type, float alt);
  int
  get_phase(float *dest);
  int
  raytrace(sutra_atmos *atmos, bool async);
  int
  raytrace(sutra_atmos *atmos);
  int
  raytrace(sutra_dms *ydms, int rst, bool async);
  int
  raytrace(sutra_dms *ydms, int rst);
  int
  raytrace_shm(sutra_atmos *atmos);
  int
  comp_image(int puponly);
  int
  init_strehlmeter();
  int
  comp_strehl();
};

class sutra_target {
public:
  int ntargets;
  vector<sutra_source *> d_targets;

public:
  sutra_target(carma_context *context, int ntargets, float *xpos, float *ypos,
      float *lambda, float *mag, long *sizes, float *pupil, int device);
  ~sutra_target();

  int
  get_phase(int ntarget, float *dest);

};

int
target_texraytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
    int Ny, float xoff, float yoff, int Ntot, cudaChannelFormatDesc channelDesc,
    int device);
int
target_raytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
    float xoff, float yoff, int block_size);
int
target_lgs_raytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
    float xoff, float yoff, float delta, int block_size);
int
target_raytrace_async(carma_streams streams, float *d_odata, float *d_idata,
    int nx, int ny, int Nx, float xoff, float yoff, int block_size);
int
target_raytrace_async(carma_host_obj<float> *phase_telemetry, float *d_odata,
    float *d_idata, int nx, int ny, int Nx, float xoff, float yoff,
    int block_size);
int
fft_goodsize(long size);
int
fill_amplipup(cuFloatComplex *amplipup, float *phase, float *mask, float scale,
    int puponly, int nx, int ny, int Nx, int device);

#endif // _SUTRA_TARGET_H_
