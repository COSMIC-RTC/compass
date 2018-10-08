#ifndef _SUTRA_SOURCE_H_
#define _SUTRA_SOURCE_H_

#include <carma.h>
#include <sutra_atmos.h>
#include <sutra_dm.h>
#include <sutra_lgs.h>
#include <sutra_source.h>
#include <sutra_telemetry.h>
#include <sutra_telescope.h>
#include <sutra_utils.h>
#include <map>
#include <string>
#include <vector>

using std::map;
using std::string;

typedef std::pair<std::string, int> type_screen;

class sutra_source {
 public:
  carma_context *current_context;
  int device;      // device #
  float posx;      // x position of target on the sky
  float posy;      // y position of target on the sky
  long npts;       // number of points in the pupil
  float mag;       // brightness of target
  float lambda;    // imaging lambda
  float zp;        // imaging zero point
  float scale;     // phase scale
  bool lgs;        // flag for lgs
  float G;         // Magnifying factor for WFS misalignment
  float thetaML;   // Pupil rotation angle
  float dx;        // WFS misalignment
  float dy;        // WFS misalignment
  string type;     // type of source : target / wfs
  int block_size;  // optimum block size of device

  float strehl_se;      // short exposure strehl
  float strehl_le;      // long exposure strehl
  float ref_strehl;     // reference for strehl computation
  int strehl_counter;   // counter for le strehl computation
  float phase_var;      // current phase variance in the pupil
  float phase_var_avg;  // average phase variance in the pupil
  int phase_var_count;  // counter for average phase variance in the pupil

  sutra_phase *d_phase;  // phase for this target
  // INTRO PHASE INSTRU
  // sutra_phase *d_phase_instru;
  //
  carma_host_obj<float> *phase_telemetry;  //
  sutra_lgs *d_lgs;                        // the lgs object
  carma_obj<float> *d_pupil;               // the pupil mask
  carma_obj<float> *d_image_se;  // the resulting image for this target
  carma_obj<float> *d_image_le;  // the long exposure image for this target
  carma_obj<cuFloatComplex>
      *d_amplipup;  // the complex amplitude in the pupil plane
  carma_obj<float>
      *d_phasepts;  // the valid phase points in the pupil (target only)
  carma_obj<int> *d_wherephase;  // positions of valid phase points in the pupil
  // (target only)
  map<type_screen, float> xoff;    // x reference for raytracing
  map<type_screen, float> yoff;    // y reference for raytracing
  carma_obj<float> *d_ncpa_phase;  // ncpa phase

 public:
  sutra_source(carma_context *context, float xpos, float ypos, float lambda,
               float mag, float zerop, long size, string type,
               carma_obj<float> *pupil, int Npts, int device);
  sutra_source(carma_context *context, float xpos, float ypos, float lambda,
               float mag, float zerop, long size, string type, int device);
  ~sutra_source();
  inline int init_source(carma_context *context, float xpos, float ypos,
                         float lambda, float mag, float zerop, long size,
                         string type, int device);
  int add_layer(string type, int idx, float xoff, float yoff);
  int remove_layer(string type, int idx);
  // int get_phase(float *h_dest);
  int raytrace(sutra_telescope *tel, bool rst = false);
  int raytrace(sutra_atmos *atmos, bool async = false);
  int raytrace(sutra_dms *ydms, bool rst = false, bool do_phase_var = true,
               bool async = false);
  int raytrace(sutra_telescope *tel, sutra_atmos *atmos, sutra_dms *ydms,
               bool do_phase_var = true, bool async = false);
  int raytrace(bool rst = false);
  int raytrace_shm(sutra_atmos *atmos);
  int comp_image(int puponly = 0, bool comp_le = true);
  int init_strehlmeter();
  int reset_strehlmeter();
  int comp_strehl();
  int reset_phase();
};

int target_texraytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
                       int Ny, float xoff, float yoff, int Ntot,
                       cudaChannelFormatDesc channelDesc, carma_device *device);
int target_raytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
                    float xoff, float yoff, float G, float thetaML, float dx,
                    float dy, int block_size);
int target_lgs_raytrace(float *d_odata, float *d_idata, int nx, int ny, int Nx,
                        float xoff, float yoff, float delta, int block_size);
int target_raytrace_async(carma_streams streams, float *d_odata, float *d_idata,
                          int nx, int ny, int Nx, float xoff, float yoff,
                          int block_size);
int target_raytrace_async(carma_host_obj<float> *phase_telemetry,
                          float *d_odata, float *d_idata, int nx, int ny,
                          int Nx, float xoff, float yoff, int block_size);
int fft_goodsize(long size);

// ATTEMPT AT ADDING PHASE_INSTRU
/*int
fill_amplipup(cuFloatComplex *amplipup, float *phase, float *phase_instru, float
*mask, float scale,
    int puponly, int nx, int ny, int Nx, carma_device *device);*/

int fill_amplipup(cuFloatComplex *amplipup, float *phase, float *mask,
                  float scale, int puponly, int nx, int ny, int Nx,
                  carma_device *device);

#endif  // _SUTRA_SOURCE_H_
