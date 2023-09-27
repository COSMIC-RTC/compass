// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_wcog.h
//! \ingroup   libsutra
//! \class     SutraCentroiderWcog
//! \brief     this class provides the centroider_wcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_WCOG_H_
#define _SUTRA_CENTROIDER_WCOG_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class SutraCentroiderWcog : public SutraCentroider<Tin, T> {
 public:
  float threshold;
  CarmaObj<float> *d_weights;

 public:
  SutraCentroiderWcog(CarmaContext *context, SutraWfs *wfs, long nvalid,
                        float offset, float scale, bool filter_TT, int device);
  SutraCentroiderWcog(const SutraCentroiderWcog &centroider);
  ~SutraCentroiderWcog();

  string get_type();

  int set_threshold(float threshold);

  int init_weights();
  int load_weights(float *weights, int ndim);

  int get_cog(float *cube, float *intensities, T *centroids, int nvalid,
              int npix, int ntot, cudaStream_t stream=0);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, float *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy,
                   float *intensities, float *d_weights, float threshold, float scale,
                   float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device, cudaStream_t stream=0);

template <class T>
int fill_weights(T *d_out, T *d_in, int npix, int N, CarmaDevice *device);
#endif  // _SUTRA_CENTROIDER_WCOG_H_
