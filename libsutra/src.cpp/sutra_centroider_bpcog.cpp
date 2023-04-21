// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_pbcog.h
//! \ingroup   libsutra
//! \class     sutra_centroider_pbcog
//! \brief     this class provides the centroider_pbcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24

#include <sutra_centroider_bpcog.h>

template <class Tin, class T>
SutraCentroiderBpcog<Tin, T>::SutraCentroiderBpcog(
    CarmaContext *context, SutraWfs *wfs, long nvalid, float offset,
    float scale, bool filter_TT, int device, int nmax)
    : SutraCentroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT,
                               device) {
  this->nslopes = 2 * nvalid;
  this->nmax = nmax;
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new CarmaObj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();

  if (filter_TT) {
    this->init_TT_filter();
  }

  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = nmax;
  dims_data[2] = nvalid;

  this->d_bpix = new CarmaObj<T>(this->current_context, dims_data);
  this->d_bpind = new CarmaObj<uint>(this->current_context, dims_data);
}

template <class Tin, class T>
SutraCentroiderBpcog<Tin, T>::~SutraCentroiderBpcog() {
  delete this->d_bpix;
  delete this->d_bpind;
}

template <class Tin, class T>
string SutraCentroiderBpcog<Tin, T>::get_type() {
  return "bpcog";
}

template <class Tin, class T>
int SutraCentroiderBpcog<Tin, T>::set_nmax(int nmax) {
  this->current_context->set_active_device(this->device, 1);
  if(nmax != this->nmax) {
    this->nmax = nmax;
    delete this->d_bpix;
    delete this->d_bpind;

    long dims_data[3];
    dims_data[0] = 2;
    dims_data[1] = nmax;
    dims_data[2] = this->nvalid;

    this->d_bpix = new CarmaObj<T>(this->current_context, dims_data);
    this->d_bpind = new CarmaObj<uint>(this->current_context, dims_data);
  }
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int SutraCentroiderBpcog<Tin, T>::get_cog(float *img, float *intensities,
                                            T *centroids, int nvalid, int npix,
                                            int ntot, cudaStream_t stream) {
  this->current_context->set_active_device(this->device, 1);

  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->get_data(), this->d_validx->get_data(),
                this->d_validy->get_data(), intensities, this->nmax, this->scale,
                this->offset, this->slope_order, this->current_context->get_device(this->device));

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }
  // brightest pixels cog
  // subap_sortmax<T>(npix * npix, nvalid, cube, this->d_bpix->get_data(),
  //                      this->d_bpind->get_data(), this->nmax,
  //                      current_context->get_device(device));

  // subap_bpcentro<T>(this->nmax, nvalid, npix, this->d_bpix->get_data(),
  //                       this->d_bpind->get_data(), centroids, this->scale,
  //                       this->offset);

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int SutraCentroiderBpcog<Tin, T>::get_cog(float *intensities, T *slopes,
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
int SutraCentroiderBpcog<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class SutraCentroiderBpcog<float, float>;
template class SutraCentroiderBpcog<uint16_t, float>;

#ifdef CAN_DO_HALF
template <>
int SutraCentroiderBpcog<float, half>::get_cog(float *intensities,
                                                 half *slopes, bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int SutraCentroiderBpcog<uint16_t, half>::get_cog(float *intensities,
                                                    half *slopes, bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int SutraCentroiderBpcog<float, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int SutraCentroiderBpcog<uint16_t, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template class SutraCentroiderBpcog<float, half>;
template class SutraCentroiderBpcog<uint16_t, half>;
#endif
