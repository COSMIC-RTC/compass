// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_pbcog.h
//! \ingroup   libsutra
//! \class     sutra_centroider_pbcog
//! \brief     this class provides the centroider_pbcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CENTROIDER_BPCOG_H_
#define _SUTRA_CENTROIDER_BPCOG_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class SutraCentroiderBpcog : public SutraCentroider<Tin, T> {
 public:
  int nmax;
  CarmaObj<T> *d_bpix;
  CarmaObj<uint> *d_bpind;

 public:
  SutraCentroiderBpcog(CarmaContext *context, SutraWfs *wfs, long nvalid,
                         float offset, float scale, bool filter_TT, int device,
                         int nmax);
  SutraCentroiderBpcog(const SutraCentroiderBpcog &centroider);
  ~SutraCentroiderBpcog();

  string get_type();

  int init_nmax(int nmax);
  int set_nmax(int nmax);

  int get_cog(float *cube, float *intensities, T *centroids, int nvalid,
              int npix, int ntot, cudaStream_t stream=0);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
};
template <class T>
void get_centroids(int size, int threads, int blocks, int n, float *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy,
                   float *intensities, int nbpix, float scale, float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device, cudaStream_t stream=0);

#endif  // _SUTRA_CENTROIDER_H_
