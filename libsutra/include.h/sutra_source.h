/**
 * \file sutra_source.h
 *
 * \class sutra_source
 *
 * \ingroup libsutra
 *
 * \brief this class provides the source features to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 4.2.0
 *
 * \date 2011/01/28
 *
 */
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
  /// device #
  int device;
  /// x position of target on the sky
  float posx;
  /// y position of target on the sky
  float posy;
  /// number of points in the pupil
  long npts;
  /// brightness of target
  float mag;
  /// imaging lambda
  float lambda;
  /// imaging zero point
  float zp;
  /// phase scale
  float scale;
  /// flag for lgs
  bool lgs;
  /// Magnifying factor for WFS misalignment
  float G;
  /// Pupil rotation angle
  float thetaML;
  /// WFS misalignment
  float dx;
  /// WFS misalignment
  float dy;
  /// type of source : target / wfs
  string type;
  /// optimum block size of device
  int block_size;

  /// short exposure strehl
  float strehl_se;
  /// long exposure strehl
  float strehl_le;
  /// reference for strehl computation
  float ref_strehl;
  /// counter for le strehl computation
  int strehl_counter;
  /// current phase variance in the pupil
  float phase_var;
  /// average phase variance in the pupil
  float phase_var_avg;
  /// counter for average phase variance in the pupil
  int phase_var_count;

  /// phase for this target
  sutra_phase *d_phase;
  // INTRO PHASE INSTRU
  // sutra_phase *d_phase_instru;
  carma_host_obj<float> *phase_telemetry;
  /// the lgs object
  sutra_lgs *d_lgs;
  /// the pupil mask
  carma_obj<float> *d_pupil;
  /// the resulting image for this target
  carma_obj<float> *d_image_se;
  /// the long exposure image for this target
  carma_obj<float> *d_image_le;
  /// the complex amplitude in the pupil plane
  carma_obj<cuFloatComplex> *d_amplipup;
  /// the valid phase points in the pupil (target only)
  carma_obj<float> *d_phasepts;
  /// positions of valid phase points in the pupil (target only)
  carma_obj<int> *d_wherephase;
  /// x reference for raytracing
  map<type_screen, float> xoff;
  /// y reference for raytracing
  map<type_screen, float> yoff;
  /// ncpa phase
  carma_obj<float> *d_ncpa_phase;
  /// temporary array for accurate strehl computation
  carma_obj<float> *d_smallimg;
  const int d_smallimg_size = 3;

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
  int comp_strehl(bool do_fit);
  int reset_phase();

 private:
  /**
   * @brief fit the strehl with a sinc
   *
   * Utilise la “croix” de 3 pixels centraux qui encadrent le max
   * pour fitter des paraboles qui determinent la position du maximum,
   * puis calcule l’interpolation exacte en ce point via la formule
   * des sinus cardinaux qui s’applique a un signal bien echantillonne.
   *
   * @param d_img full image of size img_size*img_size
   * @param ind_max position of the maximum in d_img
   * @param img_size size of the d_img leading dimension
   * @return float Strehl fitted
   */
  float fitmax2x1dSinc(float *d_img, int ind_max, int img_size);
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
