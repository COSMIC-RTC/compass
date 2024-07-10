// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_pyr.cpp
//! \ingroup   libsutra
//! \class     SutraCentroiderPyr
//! \brief     this class provides the centroider_pyr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_centroider_pyr.hpp>
#include <iostream>
#include <string>

template <class Tin, class T>
SutraCentroiderPyr<Tin, T>::SutraCentroiderPyr(CarmaContext *context,
                                                   SutraWfs *wfs, int64_t nvalid,
                                                   float offset, float scale,
                                                   bool filter_TT, int32_t device)
    : SutraCentroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT,
                               device) {
  context->set_active_device(device, 1);

  this->nslopes = 2 * nvalid;

  this->pyr_type = "pyrhr";

  this->valid_thresh = 1e-4;

  // centroider method by default nosin_global
  this->method = Method_CoG(false, false);

  this->d_intensities->init_reduceCub();
  int64_t dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new CarmaObj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();

  if (filter_TT) {
    this->init_TT_filter();
  }
}

template <class Tin, class T>
SutraCentroiderPyr<Tin, T>::~SutraCentroiderPyr() {}

template <class Tin, class T>
string SutraCentroiderPyr<Tin, T>::get_type() {
  return this->pyr_type;
}

template <class Tin, class T>
int32_t SutraCentroiderPyr<Tin, T>::set_valid_thresh(T valid_thresh) {
  this->valid_thresh = valid_thresh;
  return EXIT_SUCCESS;
}
template <class Tin, class T>
T SutraCentroiderPyr<Tin, T>::get_valid_thresh() {
  return this->valid_thresh;
}

template <class Tin, class T>
int32_t SutraCentroiderPyr<Tin, T>::set_method(Method_CoG method) {
  this->method = method;
  return EXIT_SUCCESS;
}

template <class Tin, class T>
Method_CoG SutraCentroiderPyr<Tin, T>::get_method() {
  return this->method;
}

template <class Tin, class T>
string SutraCentroiderPyr<Tin, T>::get_method_str() {
  return Method_CoG::str(this->method);
}

template <class Tin, class T>
int32_t SutraCentroiderPyr<Tin, T>::get_cog(float *cube, float *intensities,
                                          T *centroids, int32_t nvalid, int32_t npix,
                                          int32_t ntot, cudaStream_t stream) {
  // TODO(Implement SutraCentroiderPyr<Tin, T>::get_cog)

  return get_pyr(cube, intensities, centroids, this->d_validx->get_data(),
                 this->d_validy->get_data(), this->nvalid,
                 this->d_img->get_dims(1), 4, stream);
}

template <class Tin, class T>
int32_t SutraCentroiderPyr<Tin, T>::get_pyr(float *cube, float *intensities,
                                          T *centroids, int32_t *subindx,
                                          int32_t *subindy, int32_t nvalid, int32_t ns,
                                          int32_t nim, cudaStream_t stream) {
  this->current_context->set_active_device(this->device, 1);

  pyr_intensities(this->d_intensities->get_data(), cube, subindx, subindy, ns,
                  nvalid, nim, this->current_context->get_device(this->device), stream);

  if (!(this->method.is_local)) {
    // float p_sum = reduce<float>(this->d_intensities->get_data(), nvalid);
    this->d_intensities->reduceCub(stream);
    fillvalues<float>(this->d_intensities->get_data(),
                      this->d_intensities->get_o_data(), nvalid,
                      this->current_context->get_device(this->device), stream);
  }

  pyr_slopes(centroids, this->d_centroids_ref->get_data(), cube, subindx,
              subindy, this->d_intensities->get_data(), ns, nvalid, this->scale,
              this->valid_thresh,
              this->method.is_sinus,  // if we are using a sin method
              this->slope_order,
              this->current_context->get_device(this->device), stream);

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderPyr<Tin, T>::get_cog(float *intensities, T *slopes,
                                          bool noise) {
  if (this->wfs != nullptr) {
    if (this->pyr_type == "pyr" || this->pyr_type == "roof")
      return this->get_pyr(*(this->wfs->d_bincube), intensities, slopes,
                           *(this->wfs->d_validsubsx),
                           *(this->wfs->d_validsubsy), this->wfs->nvalid,
                           this->wfs->nfft / this->wfs->nrebin, 4);
    else if (this->pyr_type == "pyrhr") {
      if (noise || this->wfs->roket == false) {
        return this->get_pyr(*(this->wfs->d_binimg), intensities, slopes,
                             *(this->wfs->d_validsubsx),
                             *(this->wfs->d_validsubsy), this->wfs->nvalid,
                             this->wfs->nfft / this->wfs->nrebin, 4);
      } else
        return this->get_pyr(*(this->wfs->d_binimg_notnoisy), intensities,
                             slopes, *(this->wfs->d_validsubsx),
                             *(this->wfs->d_validsubsy), this->wfs->nvalid,
                             this->wfs->nfft / this->wfs->nrebin, 4);
    }
  }

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template <class Tin, class T>
int32_t SutraCentroiderPyr<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class SutraCentroiderPyr<float, float>;
template class SutraCentroiderPyr<uint16_t, float>;
