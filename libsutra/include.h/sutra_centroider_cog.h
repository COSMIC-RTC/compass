// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_cog.h
//! \ingroup   libsutra
//! \class     SutraCentroiderCog
//! \brief     this class provides the centroider_cog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_COG_H_
#define _SUTRA_CENTROIDER_COG_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class SutraCentroiderCog : public SutraCentroider<Tin, T> {
 public:
  SutraCentroiderCog(CarmaContext *context, SutraWfs *wfs, long nvalid,
                       float offset, float scale, bool filter_TT, int device);
  SutraCentroiderCog(const SutraCentroiderCog &centroider);
  ~SutraCentroiderCog();

  string get_type();

  int get_cog(float *cube, float *intensities, T *centroids, int nvalid,
              int npix, int ntot, cudaStream_t stream=0);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, float *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy,
                   float *intensities, float scale, float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device, cudaStream_t stream=0);

template <class T>
void get_centroids_async(int threads, int blocks, int n, CarmaStreams *streams,
                         T *d_idata, T *d_odata, T *alpha, float scale,
                         float offset);

#endif  // _SUTRA_CENTROIDER_COG_H_
