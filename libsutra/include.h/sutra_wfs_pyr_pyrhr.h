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
//! \version   5.4.2
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
  long npupils;
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
                      CarmaObj<cuFloatComplex> *d_fttotim, long nxsub,
                      long nvalid, long npupils, long npix, long nphase,
                      long nrebin, long nfft, long ntot, long npup, float pdiam,
                      float nphotons, float nphot4imat, int lgs, bool fakecam,
                      int max_flux_per_pix, int max_pix_value, bool roket,
                      int device);
  SutraWfs_PyrHR(CarmaContext *context, SutraTelescope *d_tel,
                      CarmaObj<cuFloatComplex> *d_camplipup,
                      CarmaObj<cuFloatComplex> *d_camplifoc,
                      CarmaObj<cuFloatComplex> *d_fttotim, long nxsub,
                      long nvalid, long npupils, long npix, long nphase,
                      long nrebin, long nfft, long ntot, long npup, float pdiam,
                      float nphotons, float nphot4imat, int lgs, bool fakecam,
                      int max_flux_per_pix, int max_pix_value, bool roket,
                      int nbdevices, int *devices);
  ~SutraWfs_PyrHR();

  int load_arrays(cuFloatComplex *halfxy, float *cx, float *cy, float *weights,
                 float *sincar, float *submask, int *validsubsx,
                 int *validsubsy, int *phasemap, float *fluxPerSub,
                 float *ttprojmat);
  int set_submask(float *submask);
  int set_phalfxy(cuFloatComplex *phalfxy);

  int fill_binimage(int async = 0);
  int comp_image(bool noise = true);
  void comp_modulation(int cpt);

  int copy_valid_pix(float *img, int *validx, int *validy, int im_dim);
  int set_pyr_modulation_points(float *cx, float *cy, int npts);
  int set_pyr_modulation_points(float *cx, float *cy, float *weights, int npts);
  int set_pyr_mod_weights(float *weights, int npts);

  int define_mpi_rank(int rank, int size) { return EXIT_SUCCESS; }
  int allocate_buffers(map<vector<int>, cufftHandle *> campli_plans,
                       map<vector<int>, cufftHandle *> fttotim_plans) {
    return EXIT_SUCCESS;
  }
  int comp_nphot(float ittime, float optthroughput, float diam, float cobs,
                 float zerop, float gsmag);

 private:
  int comp_generic();
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
