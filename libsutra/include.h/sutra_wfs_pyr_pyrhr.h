// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_wfs_pyr_pyrhr.h
//! \ingroup   libsutra
//! \class     SutraWfs_PyrHR
//! \brief     this class provides the wfs_pyr_pyrhr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_WFS_PYR_PYRHR_H_
#define _SUTRA_WFS_PYR_PYRHR_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_telescope.h>
#include <sutra_wfs.h>
#include <map>
#include <vector>

using std::string;
class SutraWfs_PyrHR : public SutraWfs {
 public:
  int64_t npupils;
  bool compute_pyrfocalplane;
  CarmaObj<float> *d_hrimg;
  CarmaObj<float> *d_psum;
  CarmaObj<float> *d_pyrfocalplane;
  CarmaObj<cuFloatComplex> *d_phalfxy;
  CarmaObj<cuFloatComplex> *d_poffsets;
  CarmaObj<float> *d_modu_gather;
  CarmaHostObj<float> *pyr_cx;
  CarmaHostObj<float> *pyr_cy;
  CarmaHostObj<float> *pyr_mod_weights;


 public:
  SutraWfs_PyrHR(CarmaContext *context, SutraTelescope *d_tel,
                      CarmaObj<cuFloatComplex> *d_camplipup,
                      CarmaObj<cuFloatComplex> *d_camplifoc,
                      CarmaObj<cuFloatComplex> *d_fttotim, int64_t nxsub,
                      int64_t nvalid, int64_t npupils, int64_t npix, int64_t nphase,
                      int64_t nrebin, int64_t nfft, int64_t ntot, int64_t npup, float pdiam,
                      float nphotons, float nphot4imat, int32_t lgs, bool fakecam,
                      int32_t max_flux_per_pix, int32_t max_pix_value, bool roket,
                      int32_t device);
  SutraWfs_PyrHR(CarmaContext *context, SutraTelescope *d_tel,
                      CarmaObj<cuFloatComplex> *d_camplipup,
                      CarmaObj<cuFloatComplex> *d_camplifoc,
                      CarmaObj<cuFloatComplex> *d_fttotim, int64_t nxsub,
                      int64_t nvalid, int64_t npupils, int64_t npix, int64_t nphase,
                      int64_t nrebin, int64_t nfft, int64_t ntot, int64_t npup, float pdiam,
                      float nphotons, float nphot4imat, int32_t lgs, bool fakecam,
                      int32_t max_flux_per_pix, int32_t max_pix_value, bool roket,
                      int32_t nbdevices, int32_t *devices);
  ~SutraWfs_PyrHR();

  int32_t load_arrays(cuFloatComplex *halfxy, float *cx, float *cy, float *weights,
                 float *sincar, float *submask, int32_t *validsubsx,
                 int32_t *validsubsy, int32_t *phasemap, float *fluxPerSub,
                 float *ttprojmat);
  int32_t set_submask(float *submask);
  int32_t set_phalfxy(cuFloatComplex *phalfxy);

  int32_t fill_binimage(int32_t async = 0);
  int32_t comp_image(bool noise = true);
  void comp_modulation(int32_t cpt);

  int32_t copy_valid_pix(float *img, int32_t *validx, int32_t *validy, int32_t im_dim);
  int32_t set_pyr_modulation_points(float *cx, float *cy, int32_t npts);
  int32_t set_pyr_modulation_points(float *cx, float *cy, float *weights, int32_t npts);
  int32_t set_pyr_mod_weights(float *weights, int32_t npts);

  int32_t define_mpi_rank(int32_t rank, int32_t size) { return EXIT_SUCCESS; }
  int32_t allocate_buffers(map<vector<int32_t>, cufftHandle *> campli_plans,
                       map<vector<int32_t>, cufftHandle *> fttotim_plans) {
    return EXIT_SUCCESS;
  }
  int32_t comp_nphot(float ittime, float optthroughput, float diam, float cobs,
                 float zerop, float gsmag);

 private:
  int32_t comp_generic();
  std::vector<CarmaObj<cuFloatComplex> *> d_camplipup_ngpu;
  std::vector<CarmaObj<cuFloatComplex> *> d_camplifoc_ngpu;
  std::vector<CarmaObj<cuFloatComplex> *> d_phalfxy_ngpu;
  std::vector<CarmaObj<cuFloatComplex> *> d_fttotim_ngpu;
  std::vector<CarmaObj<float> *> d_pyrfocalplane_ngpu;
  std::vector<CarmaObj<float> *> d_screen_ngpu;
  std::vector<CarmaObj<float> *> d_hrimg_ngpu;
  std::vector<CarmaObj<float> *> d_submask_ngpu;
};

#endif  // _SUTRA_WFS_PYR_PYRHR_H_
