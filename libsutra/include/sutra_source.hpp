// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_source.hpp
//! \ingroup   libsutra
//! \class     SutraSource
//! \brief     this class provides the source features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_SOURCE_H_
#define _SUTRA_SOURCE_H_

#include <carma.hpp>
#include <sutra_atmos.hpp>
#include <sutra_dm.hpp>
#include <sutra_lgs.hpp>
#include <sutra_source.hpp>
#include <sutra_telemetry.hpp>
#include <sutra_telescope.hpp>
#include <sutra_utils.hpp>
#include <map>
#include <string>
#include <vector>

using std::map;
using std::string;

typedef std::pair<std::string, int32_t> type_screen;

class SutraSource {
 public:
  CarmaContext *current_context;
  /// device #
  int32_t device;
  /// x position of target on the sky
  float posx;
  /// y position of target on the sky
  float posy;
  /// number of points in the pupil
  int64_t npts;
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
  int32_t block_size;

  /// short exposure strehl
  float strehl_se;
  /// int64_t exposure strehl
  float strehl_le;
  /// reference for strehl computation
  float ref_strehl;
  /// counter for le strehl computation
  int32_t strehl_counter;
  /// current phase variance in the pupil
  float phase_var;
  /// average phase variance in the pupil
  float phase_var_avg;
  /// counter for average phase variance in the pupil
  int32_t phase_var_count;

  /// phase for this target
  SutraPhase *d_phase;
  // INTRO PHASE INSTRU
  // SutraPhase *d_phase_instru;
  CarmaHostObj<float> *phase_telemetry;
  /// the lgs object
  SutraLGS *d_lgs;
  /// the pupil mask
  CarmaObj<float> *d_pupil;
  /// the resulting image for this target
  CarmaObj<float> *d_image_se;
  /// the int64_t exposure image for this target
  CarmaObj<float> *d_image_le;
  /// the complex amplitude in the pupil plane
  CarmaObj<cuFloatComplex> *d_amplipup;
  /// the valid phase points in the pupil (target only)
  CarmaObj<float> *d_phasepts;
  /// positions of valid phase points in the pupil (target only)
  CarmaObj<int32_t> *d_wherephase;
  /// x reference for raytracing
  map<type_screen, float> xoff;
  /// y reference for raytracing
  map<type_screen, float> yoff;
  /// ncpa phase
  CarmaObj<float> *d_ncpa_phase;
  /// temporary array for accurate strehl computation
  CarmaObj<float> *d_smallimg;
  const int32_t d_smallimg_size = 3;

 public:
  SutraSource(CarmaContext *context, float xpos, float ypos, float lambda,
               float mag, float zerop, int64_t size, string type,
               CarmaObj<float> *pupil, int32_t Npts, int32_t device);
  SutraSource(CarmaContext *context, float xpos, float ypos, float lambda,
               float mag, float zerop, int64_t size, string type, int32_t device);
  ~SutraSource();
  inline int32_t init_source(CarmaContext *context, float xpos, float ypos,
                         float lambda, float mag, float zerop, int64_t size,
                         string type, int32_t device);
  int32_t add_layer(string type, int32_t idx, float xoff, float yoff);
  int32_t remove_layer(string type, int32_t idx);
  // int32_t get_phase(float *h_dest);
  int32_t raytrace(SutraTelescope *tel, bool rst = false);
  int32_t raytrace(SutraAtmos *atmos, bool async = false);
  int32_t raytrace(SutraDms *ydms, bool rst = false, bool do_phase_var = true,
               bool async = false);
  int32_t raytrace(SutraTelescope *tel, SutraAtmos *atmos, SutraDms *ydms,
               bool do_phase_var = true, bool async = false);
  int32_t raytrace(bool rst = false);
  int32_t raytrace_shm(SutraAtmos *atmos);
  int32_t comp_image(int32_t puponly = 0, bool comp_le = true);
  int32_t init_strehlmeter();
  int32_t reset_strehlmeter();
  int32_t comp_strehl(bool do_fit);
  int32_t reset_phase();

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
  float fitmax2x1dSinc(float *d_img, int32_t ind_max, int32_t img_size);
};

int32_t target_raytrace(float *d_odata, float *d_idata, int32_t nx, int32_t ny, int32_t Nx,
                    float xoff, float yoff, float G, float thetaML, float dx,
                    float dy, int32_t block_size, float delta);
int32_t target_raytrace_async(CarmaStreams streams, float *d_odata, float *d_idata,
                          int32_t nx, int32_t ny, int32_t Nx, float xoff, float yoff,
                          int32_t block_size);
int32_t target_raytrace_async(CarmaHostObj<float> *phase_telemetry,
                          float *d_odata, float *d_idata, int32_t nx, int32_t ny,
                          int32_t Nx, float xoff, float yoff, int32_t block_size);
int32_t fft_goodsize(int64_t size);

// ATTEMPT AT ADDING PHASE_INSTRU
/*int32_t
fill_amplipup(cuFloatComplex *amplipup, float *phase, float *phase_instru, float
*mask, float scale,
    int32_t puponly, int32_t nx, int32_t ny, int32_t Nx, CarmaDevice *device);*/

int32_t fill_amplipup(cuFloatComplex *amplipup, float *phase, float *mask,
                  float scale, int32_t puponly, int32_t nx, int32_t ny, int32_t Nx,
                  CarmaDevice *device);

#endif  // _SUTRA_SOURCE_H_
