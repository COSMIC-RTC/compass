// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_wfs_sh.h
//! \ingroup   libsutra
//! \class     SutraWfsSH
//! \brief     this class provides the wfs_sh features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24

#ifndef _SUTRA_WFS_SH_H_
#define _SUTRA_WFS_SH_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>
#include <map>
#include <vector>

class SutraWfsSH : public SutraWfs {
 public:
  // sh only
  CarmaObj<int> *d_binmap;
  CarmaObj<int> *d_validpuppixx;  // nxsub
  CarmaObj<int> *d_validpuppixy;  // nxsub
  CarmaObj<cuFloatComplex> *d_fsamplipup; // Field stop computation arrays
  CarmaObj<cuFloatComplex> *d_fsamplifoc; // Field stop computation arrays
  cufftHandle *fsampli_plan;

 public:
  SutraWfsSH(CarmaContext *context, SutraTelescope *d_tel,
               CarmaObj<cuFloatComplex> *d_camplipup,
               CarmaObj<cuFloatComplex> *d_camplifoc,
               CarmaObj<cuFloatComplex> *d_fttotim, long nxsub, long nvalid,
               long npix, long nphase, long nrebin, long nfft, long ntot,
               long npup, float pdiam, float nphotons, float nphot4imat,
               int lgs, bool fakecam, int max_flux_per_pix, int max_pix_value,
               bool is_low_order, bool roket, int device);
  SutraWfsSH(const SutraWfsSH &wfs);
  ~SutraWfsSH();

  int define_mpi_rank(int rank, int size);
  int allocate_buffers(map<vector<int>, cufftHandle *> campli_plans,
                       map<vector<int>, cufftHandle *> fttotim_plans);

  int load_arrays(int *phasemap, int *hrmap, int *binmap, float *offsets,
                 float *fluxPerSub, int *validsubsx, int *validsubsy,
                 int *istart, int *jstart, float *ttprojmat,
                cuFloatComplex *kernel);

  int fill_binimage(int async);
  int comp_image(bool noise = true);
  int comp_nphot(float ittime, float optthroughput, float diam, int nxsub,
                 float zerop = 0, float gsmag = 0, float lgsreturnperwatt = 0,
                 float laserpower = 0);
  int set_bincube(float *bincube, int nElem);
  int set_field_stop(map<vector<int>, cufftHandle *> campli_plans, float* field_stop, int N);

 private:
  int comp_generic();
};

#endif  // _SUTRA_WFS_SH_H_
