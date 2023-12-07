// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_tcog.h
//! \ingroup   libsutra
//! \class     SutraCentroiderTcog
//! \brief     this class provides the centroider_tcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_TCOG_H_
#define _SUTRA_CENTROIDER_TCOG_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class SutraCentroiderTcog : public SutraCentroider<Tin, T> {
 public:
  float threshold;

 public:
  SutraCentroiderTcog(CarmaContext *context, SutraWfs *wfs, int64_t nvalid,
                        float offset, float scale, bool filter_TT, int32_t device);
  SutraCentroiderTcog(const SutraCentroiderTcog &centroider);
  ~SutraCentroiderTcog();

  string get_type();

  int32_t set_threshold(float threshold);

  int32_t get_cog(float *cube, float *intensities, T *centroids, int32_t nvalid,
              int32_t npix, int32_t ntot, cudaStream_t stream=0);
  int32_t get_cog(float *intensities, T *slopes, bool noise);
  int32_t get_cog();
};

template <class T>
void get_centroids(int32_t size, int32_t threads, int32_t blocks, int32_t n, float *d_idata,
                   T *d_odata, T *ref, int32_t *validx, int32_t *validy,
                   float *intensities, float threshold, float scale,
                   float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device, cudaStream_t stream=0);

#endif  // _SUTRA_CENTROIDER_TCOG_H_
