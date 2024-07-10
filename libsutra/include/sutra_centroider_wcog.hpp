// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_wcog.hpp
//! \ingroup   libsutra
//! \class     SutraCentroiderWcog
//! \brief     this class provides the centroider_wcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_WCOG_H_
#define _SUTRA_CENTROIDER_WCOG_H_

#include <sutra_centroider.hpp>

template <class Tin, class T>
class SutraCentroiderWcog : public SutraCentroider<Tin, T> {
 public:
  float threshold;
  CarmaObj<float> *d_weights;

 public:
  SutraCentroiderWcog(CarmaContext *context, SutraWfs *wfs, int64_t nvalid,
                        float offset, float scale, bool filter_TT, int32_t device);
  SutraCentroiderWcog(const SutraCentroiderWcog &centroider);
  ~SutraCentroiderWcog();

  string get_type();

  int32_t set_threshold(float threshold);

  int32_t init_weights();
  int32_t load_weights(float *weights, int32_t ndim);

  int32_t get_cog(float *cube, float *intensities, T *centroids, int32_t nvalid,
              int32_t npix, int32_t ntot, cudaStream_t stream=0);
  int32_t get_cog(float *intensities, T *slopes, bool noise);
  int32_t get_cog();
};

template <class T>
void get_centroids(int32_t size, int32_t threads, int32_t blocks, int32_t n, float *d_idata,
                   T *d_odata, T *ref, int32_t *validx, int32_t *validy,
                   float *intensities, float *d_weights, float threshold, float scale,
                   float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device, cudaStream_t stream=0);

template <class T>
int32_t fill_weights(T *d_out, T *d_in, int32_t npix, int32_t N, CarmaDevice *device);
#endif  // _SUTRA_CENTROIDER_WCOG_H_
