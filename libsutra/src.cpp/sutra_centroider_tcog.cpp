// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_tcog.cpp
//! \ingroup   libsutra
//! \class     SutraCentroiderTcog
//! \brief     this class provides the centroider_tcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_centroider_tcog.h>
#include <string>

template <class Tin, class T>
SutraCentroiderTcog<Tin, T>::SutraCentroiderTcog(CarmaContext *context,
                                                     SutraWfs *wfs,
                                                     long nvalid, float offset,
                                                     float scale,
                                                     bool filter_TT, int device)
    : SutraCentroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT,
                               device) {
  context->set_active_device(device, 1);

  this->nslopes = 2 * nvalid;
  this->threshold = 0.f;
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new CarmaObj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();

  if (filter_TT) {
    this->init_TT_filter();
  }
}

template <class Tin, class T>
SutraCentroiderTcog<Tin, T>::~SutraCentroiderTcog() {}

template <class Tin, class T>
string SutraCentroiderTcog<Tin, T>::get_type() {
  return "tcog";
}

template <class Tin, class T>
int SutraCentroiderTcog<Tin, T>::set_threshold(float threshold) {
  this->threshold = threshold;

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int SutraCentroiderTcog<Tin, T>::get_cog(float *img, float *intensities,
                                           T *centroids, int nvalid, int npix,
                                           int ntot, cudaStream_t stream) {
  this->current_context->set_active_device(this->device, 1);

  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->get_data(), this->d_validx->get_data(),
                this->d_validy->get_data(), intensities, this->threshold,
                this->scale, this->offset,
                this->slope_order,
                this->current_context->get_device(this->device), stream);

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int SutraCentroiderTcog<Tin, T>::get_cog(float *intensities, T *slopes,
                                           bool noise) {
  if (this->wfs != nullptr) {
    if (noise || this->wfs->roket == false)
      return this->get_cog(*(this->wfs->d_binimg), intensities, slopes,
                           this->wfs->nvalid, this->wfs->npix,
                           this->wfs->d_binimg->get_dims()[1]);
    else
      return this->get_cog(*(this->wfs->d_binimg_notnoisy), intensities, slopes,
                           this->wfs->nvalid, this->wfs->npix,
                           this->wfs->d_binimg->get_dims()[1]);
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template <class Tin, class T>
int SutraCentroiderTcog<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class SutraCentroiderTcog<float, float>;
template class SutraCentroiderTcog<uint16_t, float>;

#ifdef CAN_DO_HALF
template <>
int SutraCentroiderTcog<float, half>::get_cog(float *intensities,
                                                half *slopes, bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int SutraCentroiderTcog<uint16_t, half>::get_cog(float *intensities,
                                                   half *slopes, bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int SutraCentroiderTcog<float, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int SutraCentroiderTcog<uint16_t, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template class SutraCentroiderTcog<float, half>;
template class SutraCentroiderTcog<uint16_t, half>;
#endif
