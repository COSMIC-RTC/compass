// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_wfs_pyr_pyrhr.h
//! \ingroup   libsutra
//! \class     sutra_wfs_pyr_pyrhr
//! \brief     this class provides the wfs_pyr_pyrhr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
class sutra_wfs_pyr_pyrhr : public sutra_wfs {
 public:
  long npupils;
  bool compute_pyrfocalplane;
  carma_obj<float> *d_hrimg;
  carma_obj<float> *d_submask;
  carma_obj<float> *d_psum;
  carma_obj<float> *d_pyrfocalplane;
  carma_obj<cuFloatComplex> *d_phalfxy;
  carma_obj<cuFloatComplex> *d_poffsets;
  carma_obj<float> *d_modu_gather;
  carma_host_obj<float> *pyr_cx;
  carma_host_obj<float> *pyr_cy;
  carma_host_obj<float> *pyr_mod_weights;


 public:
  sutra_wfs_pyr_pyrhr(carma_context *context, sutra_telescope *d_tel,
                      carma_obj<cuFloatComplex> *d_camplipup,
                      carma_obj<cuFloatComplex> *d_camplifoc,
                      carma_obj<cuFloatComplex> *d_fttotim, long nxsub,
                      long nvalid, long npupils, long npix, long nphase,
                      long nrebin, long nfft, long ntot, long npup, float pdiam,
                      float nphotons, float nphot4imat, int lgs, bool fakecam,
                      int maxFluxPerPix, int maxPixValue, bool roket,
                      int device);
  sutra_wfs_pyr_pyrhr(carma_context *context, sutra_telescope *d_tel,
                      carma_obj<cuFloatComplex> *d_camplipup,
                      carma_obj<cuFloatComplex> *d_camplifoc,
                      carma_obj<cuFloatComplex> *d_fttotim, long nxsub,
                      long nvalid, long npupils, long npix, long nphase,
                      long nrebin, long nfft, long ntot, long npup, float pdiam,
                      float nphotons, float nphot4imat, int lgs, bool fakecam,
                      int maxFluxPerPix, int maxPixValue, bool roket,
                      int nbdevices, int *devices);
  ~sutra_wfs_pyr_pyrhr();

  int loadarrays(cuFloatComplex *halfxy, float *cx, float *cy, float *weights,
                 float *sincar, float *submask, int *validsubsx,
                 int *validsubsy, int *phasemap, float *fluxPerSub);
  int set_submask(float *submask);
  int set_phalfxy(cuFloatComplex *phalfxy);

  int fill_binimage(int async = 0);
  int comp_image(bool noise = true);
  void comp_modulation(int cpt);

  int copyValidPix(float *img, int *validx, int *validy, int im_dim);
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
  std::vector<carma_obj<cuFloatComplex> *> d_camplipup_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_camplifoc_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_phalfxy_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_fttotim_ngpu;
  std::vector<carma_obj<float> *> d_pyrfocalplane_ngpu;
  std::vector<carma_obj<float> *> d_screen_ngpu;
  std::vector<carma_obj<float> *> d_hrimg_ngpu;
  std::vector<carma_obj<float> *> d_submask_ngpu;
};

#endif  // _SUTRA_WFS_PYR_PYRHR_H_
