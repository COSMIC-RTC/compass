// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_centroider_tcog.hpp
//! \ingroup   libsutra
//! \class     SutraCentroiderTcog
//! \brief     this class provides the centroider_tcog features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_TCOG_H_
#define _SUTRA_CENTROIDER_TCOG_H_

#include <sutra_centroider.hpp>

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
