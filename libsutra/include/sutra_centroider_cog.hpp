// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_cog.hpp
//! \ingroup   libsutra
//! \class     SutraCentroiderCog
//! \brief     this class provides the centroider_cog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_COG_H_
#define _SUTRA_CENTROIDER_COG_H_

#include <sutra_centroider.hpp>

template <class Tin, class T>
class SutraCentroiderCog : public SutraCentroider<Tin, T> {
 public:
  SutraCentroiderCog(CarmaContext *context, SutraWfs *wfs, int64_t nvalid,
                       float offset, float scale, bool filter_TT, int32_t device);
  SutraCentroiderCog(const SutraCentroiderCog &centroider);
  ~SutraCentroiderCog();

  string get_type();

  int32_t get_cog(float *cube, float *intensities, T *centroids, int32_t nvalid,
              int32_t npix, int32_t ntot, cudaStream_t stream=0);
  int32_t get_cog(float *intensities, T *slopes, bool noise);
  int32_t get_cog();
};

template <class T>
void get_centroids(int32_t size, int32_t threads, int32_t blocks, int32_t n, float *d_idata,
                   T *d_odata, T *ref, int32_t *validx, int32_t *validy,
                   float *intensities, float scale, float offset,
                   SlopeOrder slope_order,
                   CarmaDevice *device, cudaStream_t stream=0);

template <class T>
void get_centroids_async(int32_t threads, int32_t blocks, int32_t n, CarmaStreams *streams,
                         T *d_idata, T *d_odata, T *alpha, float scale,
                         float offset);

#endif  // _SUTRA_CENTROIDER_COG_H_
