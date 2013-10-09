#ifndef _YOGA_TARGET_H_
#define _YOGA_TARGET_H_

#include <map>
#include <string>
#include <vector>
#include <yoga.h>
#include <yoga_turbu.h>
#include <yoga_lgs.h>
#include <yoga_dm.h>
#include <yoga_telemetry.h>

using namespace std;
using std::string;

class yoga_source {
 public:

  int                       device;      // device # 
  float                     tposx;       // x position of target on the sky  
  float                     tposy;       // y position of target on the sky
  long                      npos;        // number of points in the pupil
  float                     mag;         // brightness of target
  float                     lambda;      // imaging lambda
  float                     zp;          // imaging zero point
  float                     scale;       // phase scale
  bool                      lgs;         // flag for lgs
  string                    type;        // type of source : target / wfs
  int                       block_size;  // optimum block size of device

  float                     strehl_se;       // short exposure strehl
  float                     strehl_le;       // long exposure strehl
  float                     ref_strehl;      // reference for strehl computation
  int                       strehl_counter;  // counter for le strehl computation
  float                     phase_var;       // current phase variance in the pupil
  float                     phase_var_avg;   // average phase variance in the pupil
  int                       phase_var_count; // counter for average phase variance in the pupil

  yoga_phase                *d_phase;    // phase for this target
  yoga_host_obj<float>      *phase_telemetry;  // 
  yoga_lgs                  *d_lgs;      // the lgs object
  yoga_obj<float>           *object;     // the object intensity map
  yoga_obj<float>           *d_pupil;    // the pupil mask
  yoga_obj<float>           *d_image;    // the resulting image for this target
  yoga_obj<float>           *d_leimage;  // the long exposure image for this target
  yoga_obj<cuFloatComplex>  *d_amplipup; // the complex amplitude in the pupil plane
  yoga_obj<float>           *d_phasepts; // the valid phase points in the pupil (target only)
  yoga_obj<int>             *d_wherephase;  // positions of valid phase points in the pupil (target only)
  map<type_screen,float>    xoff;        // x reference for raytracing
  map<type_screen,float>    yoff;        // y reference for raytracing
  yoga_context *current_context;

 public:
  yoga_source(yoga_context *context, float xpos,float ypos,float lambda,float mag,long size,string type, float *pupil, int device);
  yoga_source(yoga_context *context, float xpos,float ypos,float lambda,float mag,long size,string type, int device);
  yoga_source(const yoga_source& source);
  ~yoga_source();
  int init_source(yoga_context *context, float xpos,float ypos,float lambda,float mag,long size,string type, int device);
  int add_layer(string type,float alt,float xoff, float yoff);
  int remove_layer(string type,float alt);
  int get_phase(float *dest);
  int raytrace(yoga_atmos *atmos, bool async);
  int raytrace(yoga_atmos *atmos);
  int raytrace(yoga_dms *ydms, int rst, bool async);
  int raytrace(yoga_dms *ydms, int rst);
  int raytrace_shm(yoga_atmos *atmos);
  int comp_image(int puponly);
  int init_strehlmeter();
  int comp_strehl();
};

class yoga_target {
 public:
  int                    ntargets;
  vector<yoga_source *>  d_targets;

 public:
  yoga_target(yoga_context *context, int ntargets,float *xpos,float *ypos,float *lambda,float *mag,long *sizes, float *pupil, int device);
  ~yoga_target();

  int get_phase(int ntarget,float *dest);

};

int target_texraytrace(float *d_odata,float *d_idata,int nx, int ny,int Nx, int Ny, float xoff, 
		       float yoff, int Ntot, cudaChannelFormatDesc channelDesc, int device);
int target_raytrace(float *d_odata,float *d_idata,int nx, int ny,int Nx,float xoff, float yoff,  
		    int block_size);
int target_raytrace_async(yoga_streams streams, float *d_odata,float *d_idata,int nx, int ny,int Nx,float xoff, float yoff,
		    int block_size);
int target_raytrace_async(yoga_host_obj<float> *phase_telemetry,float *d_odata,float *d_idata,
			  int nx, int ny,int Nx,float xoff, float yoff,int block_size);
int fft_goodsize(long size);
int fill_amplipup(cuFloatComplex *amplipup, float *phase, float *mask, float scale, int puponly, int nx, int ny, int Nx, int device);

#endif // _YOGA_TARGET_H_
