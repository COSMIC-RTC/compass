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
                             long nvalid, float offset,
                             float scale, bool filter_TT, int device);

  ~SutraCentroiderMaskedPix();

  string get_type();

  int get_maskedPix(float *img, float *intensities, T *centroids, int *subindx,
                    int *subindy, int ns, cudaStream_t stream=0);
  int get_cog(float *img, float *intensities, T *centroids, int nvalid,
              int npix, int ntot, cudaStream_t stream=0);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
  int fill_selected_pix(T *pix);
  int fill_mask();
};

void fill_intensities(float *intensities, float *img, int *subindx,
                      int *subindy, int ns, int nslopes, CarmaDevice *device, cudaStream_t stream=0);
template <class T>
void get_masked_pix(T *centroids, T *ref, float *img, int *subindx, int *subindy,
                  float *psum, int ns, int nslopes, CarmaDevice *device, cudaStream_t stream=0);
template <class T>
void pyr_fill_selected_pix(T *img, int img_sizex, T *pix, int *subindx, int *subindy,
                      int nvalid, CarmaDevice *device);
template <class T>
void pyr_fill_mask(T *mask, int img_sizex, int *subindx, int *subindy,
                      int nvalid, CarmaDevice *device);

#endif  // _SUTRA_CENTROIDER_MASKEDPIX_H_
