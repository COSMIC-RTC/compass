// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_tcog.h
//! \ingroup   libsutra
//! \class     SutraCentroiderTcog
//! \brief     this class provides the centroider_tcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_TCOG_H_
#define _SUTRA_CENTROIDER_TCOG_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class SutraCentroiderTcog : public SutraCentroider<Tin, T> {
 public:
  float threshold;

 public:
  SutraCentroiderTcog(CarmaContext *context, SutraWfs *wfs, long nvalid,
                        float offset, float scale, bool filter_TT, int device);
  SutraCentroiderTcog(const SutraCentroiderTcog &centroider);
  ~SutraCentroiderTcog();

  string get_type();

  int set_threshold(float threshold);

  int get_cog(float *cube, float *intensities, T *centroids, int nvalid,
              int npix, int ntot, cudaStream_t stream=0);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, float *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy,
                   float *intensities, float threshold, float scale,
                   float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device);

#endif  // _SUTRA_CENTROIDER_TCOG_H_
