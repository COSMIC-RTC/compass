// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_wcog.cpp
//! \ingroup   libsutra
//! \class     SutraCentroiderWcog
//! \brief     this class provides the centroider_wcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_centroider_wcog.h>
#include <string>

template <class Tin, class T>
SutraCentroiderWcog<Tin, T>::SutraCentroiderWcog(CarmaContext *context,
                                                     SutraWfs *wfs,
                                                     int64_t nvalid, float offset,
                                                     float scale,
                                                     bool filter_TT, int32_t device)
    : SutraCentroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT,
                               device) {
  context->set_active_device(device, 1);

  this->nslopes = 2 * nvalid;
  this->d_weights = 0L;
  int64_t dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new CarmaObj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();
  this->threshold = 0;

  if (filter_TT) {
    this->init_TT_filter();
  }
}

template <class Tin, class T>
SutraCentroiderWcog<Tin, T>::~SutraCentroiderWcog() {}

template <class Tin, class T>
string SutraCentroiderWcog<Tin, T>::get_type() {
  return "wcog";
}

template <class Tin, class T>
int32_t SutraCentroiderWcog<Tin, T>::init_weights() {
  this->current_context->set_active_device(this->device, 1);
  if (this->d_weights != 0L) delete this->d_weights;

  int64_t *dims_data3 = new int64_t[4];
  dims_data3[0] = 3;
  dims_data3[1] = this->npix;
  dims_data3[2] = this->npix;
  dims_data3[3] = this->nvalid;

  this->current_context->set_active_device(this->device, 1);
  this->d_weights = new CarmaObj<float>(this->current_context, dims_data3);

  delete[] dims_data3;

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderWcog<Tin, T>::load_weights(float *weights, int32_t ndim) {
  if (ndim == 3)
    this->d_weights->host2device(weights);
  else {
    // weights is a 2d array
    // same weight for each subap
    float *tmp;  ///< Input data
    carma_safe_call(
        cudaMalloc((void **)&tmp, sizeof(float) * this->npix * this->npix));
    carma_safe_call(cudaMemcpy(tmp, weights,
                             sizeof(float) * this->npix * this->npix,
                             cudaMemcpyHostToDevice));
    fill_weights<float>(*(this->d_weights), tmp, this->npix,
                       this->d_weights->get_nb_elements(),
                       this->current_context->get_device(this->device));
    carma_safe_call(cudaFree(tmp));
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderWcog<Tin, T>::set_threshold(float threshold) {
  this->threshold = threshold;

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderWcog<Tin, T>::get_cog(float *img, float *intensities,
                                           T *centroids, int32_t nvalid, int32_t npix,
                                           int32_t ntot, cudaStream_t stream) {
  // wcog
  // TODO: Implement SutraCentroiderWcog<Tin, T>::get_cog_async
  // subap_reduce<T>(ntot, npix * npix, nvalid, cube, intensities,
  //                     *(this->d_weights),
  //                     this->current_context->get_device(device));

  // get_centroids<T>(ntot, npix * npix, nvalid, npix, cube, centroids,
  //                      intensities, *(this->d_weights), this->scale,
  //                      this->offset,
  //                      this->current_context->get_device(device));
  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->get_data(), this->d_validx->get_data(),
                this->d_validy->get_data(), intensities,
                this->d_weights->get_data(), this->threshold, this->scale, this->offset,
                this->slope_order,
                this->current_context->get_device(this->device), stream);

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderWcog<Tin, T>::get_cog(float *intensities, T *slopes,
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
int32_t SutraCentroiderWcog<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class SutraCentroiderWcog<float, float>;
template class SutraCentroiderWcog<uint16_t, float>;
