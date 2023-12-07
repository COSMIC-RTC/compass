// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
// -----------------------------------------------------------------------------

//! \file      sutra_centroider_maskedPix.h
//! \ingroup   libsutra
//! \class     SutraCentroiderMaskedPix
//! \brief     this class provides the centroider_maskedPix features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CENTROIDER_MASKEDPIX_H_
#define _SUTRA_CENTROIDER_MASKEDPIX_H_

#include <sutra_centroider.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <string>

template <class Tin, class T>
class SutraCentroiderMaskedPix : public SutraCentroider<Tin, T> {
 public:
  CarmaObj<T> *d_selected_pix;
  CarmaObj<T> *d_mask;

 public:
  SutraCentroiderMaskedPix(CarmaContext *context, SutraWfs *wfs,
                             int64_t nvalid, float offset,
                             float scale, bool filter_TT, int32_t device);

  ~SutraCentroiderMaskedPix();

  string get_type();

  int32_t get_maskedPix(float *img, float *intensities, T *centroids, int32_t *subindx,
                    int32_t *subindy, int32_t ns, cudaStream_t stream=0);
  int32_t get_cog(float *img, float *intensities, T *centroids, int32_t nvalid,
              int32_t npix, int32_t ntot, cudaStream_t stream=0);
  int32_t get_cog(float *intensities, T *slopes, bool noise);
  int32_t get_cog();
  int32_t fill_selected_pix(T *pix);
  int32_t fill_mask();
};

void fill_intensities(float *intensities, float *img, int32_t *subindx,
                      int32_t *subindy, int32_t ns, int32_t nslopes, CarmaDevice *device, cudaStream_t stream=0);
template <class T>
void get_masked_pix(T *centroids, T *ref, float *img, int32_t *subindx, int32_t *subindy,
                  float *psum, int32_t ns, int32_t nslopes, CarmaDevice *device, cudaStream_t stream=0);
template <class T>
void pyr_fill_selected_pix(T *img, int32_t img_sizex, T *pix, int32_t *subindx, int32_t *subindy,
                      int32_t nvalid, CarmaDevice *device);
template <class T>
void pyr_fill_mask(T *mask, int32_t img_sizex, int32_t *subindx, int32_t *subindy,
                      int32_t nvalid, CarmaDevice *device);

#endif  // _SUTRA_CENTROIDER_MASKEDPIX_H_
