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

//! \file      sutra_centroider_maskedPix.cpp
//! \ingroup   libsutra
//! \class     SutraCentroiderMaskedPix
//! \brief     this class provides the centroider_maskedPix features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_centroider_maskedPix.hpp>

#include <string>

template <class Tin, class T>
SutraCentroiderMaskedPix<Tin, T>::SutraCentroiderMaskedPix(
    CarmaContext *context, SutraWfs *wfs, int64_t nvalid, float offset,
    float scale, bool filter_TT, int32_t device)
    : SutraCentroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT,
                              device) {
  context->set_active_device(device, 1);

  if (wfs != nullptr)
    this->nslopes = wfs->d_validsubsx->get_nb_elements();
  else
    this->nslopes = nvalid;
  int64_t dims_data[2] = {1, this->nslopes};
  if (this->d_intensities != nullptr) delete this->d_intensities;
  this->d_intensities = new CarmaObj<float>(context, dims_data);
  this->d_intensities->init_reduceCub();
  int64_t dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new CarmaObj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->memset(1);
  this->d_selected_pix = nullptr;
  this->d_mask = nullptr;

  if (filter_TT) {
    this->init_TT_filter();
  }
}

template <class Tin, class T>
SutraCentroiderMaskedPix<Tin, T>::~SutraCentroiderMaskedPix() {}

template <class Tin, class T>
string SutraCentroiderMaskedPix<Tin, T>::get_type() {
  return "maskedpix";  // TODO: Fix an use constants !!!
}

template <class Tin, class T>
int32_t SutraCentroiderMaskedPix<Tin, T>::fill_selected_pix(T *pix) {
  this->current_context->set_active_device(this->device, 1);

  const int64_t *dims_data2;
  if (this->d_img != nullptr) {
    dims_data2 = this->d_img->get_dims();
  } else if (this->wfs != nullptr) {
    dims_data2 = this->wfs->d_binimg->get_dims();
  } else {
    std::cerr << "Image not initialized" << std::endl;
    return EXIT_FAILURE;
  }
  if (this->d_selected_pix == nullptr) {
    this->d_selected_pix = new CarmaObj<T>(this->current_context, dims_data2);
  }
  this->d_selected_pix->reset();
  pyr_fill_selected_pix(this->d_selected_pix->get_data(), dims_data2[1], pix,
                        this->d_validx->get_data(), this->d_validy->get_data(),
                        this->d_validx->get_nb_elements(),
                        this->current_context->get_device(this->device));
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderMaskedPix<Tin, T>::fill_mask() {
  this->current_context->set_active_device(this->device, 1);

  const int64_t *dims_data2;
  if (this->d_img != nullptr) {
    dims_data2 = this->d_img->get_dims();
  } else if (this->wfs != nullptr) {
    dims_data2 = this->wfs->d_binimg->get_dims();
  } else {
    std::cerr << "Image not initialized" << std::endl;
    return EXIT_FAILURE;
  }
  if (this->d_mask == nullptr) {
    this->d_mask = new CarmaObj<T>(this->current_context, dims_data2);
  }
  this->d_mask->reset();
  pyr_fill_mask(this->d_mask->get_data(), dims_data2[1],
                this->d_validx->get_data(), this->d_validy->get_data(),
                this->d_validx->get_nb_elements(),
                this->current_context->get_device(this->device));
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderMaskedPix<Tin, T>::get_cog(float *img, float *intensities,
                                              T *centroids, int32_t nvalid,
                                              int32_t npix, int32_t ntot,
                                              cudaStream_t stream) {
  // TODO(Implement SutraCentroiderMaskedPix<Tin, T>::get_cog)

  return get_maskedPix(img, intensities, centroids, this->d_validx->get_data(),
                       this->d_validy->get_data(), ntot, stream);
}

template <class Tin, class T>
int32_t SutraCentroiderMaskedPix<Tin, T>::get_maskedPix(float *img,
                                                    float *intensities,
                                                    T *centroids, int32_t *subindx,
                                                    int32_t *subindy, int32_t ns,
                                                    cudaStream_t stream) {
  this->current_context->set_active_device(this->device, 1);

  fill_intensities(this->d_intensities->get_data(), img, subindx, subindy, ns,
                   this->nslopes,
                   this->current_context->get_device(this->device), stream);

  // T p_sum = reduce<T>(this->d_intensities->get_data(), this->nslopes);
  this->d_intensities->reduceCub(stream);

  get_masked_pix<T>(centroids, this->d_centroids_ref->get_data(), img, subindx,
                    subindy, this->d_intensities->get_o_data(), ns,
                    this->nslopes,
                    this->current_context->get_device(this->device), stream);

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderMaskedPix<Tin, T>::get_cog(float *intensities, T *slopes,
                                              bool noise) {
  if (this->wfs != nullptr) {
    if (this->wfs->type == "pyrhr") {
      if (noise || this->wfs->roket == false) {
        return this->get_maskedPix(*(this->wfs->d_binimg), intensities, slopes,
                                   *(this->wfs->d_validsubsx),
                                   *(this->wfs->d_validsubsy),
                                   this->wfs->nfft / this->wfs->nrebin);
      } else
        return this->get_maskedPix(*(this->wfs->d_binimg_notnoisy), intensities,
                                   slopes, *(this->wfs->d_validsubsx),
                                   *(this->wfs->d_validsubsy),
                                   this->wfs->nfft / this->wfs->nrebin);
    } else
      DEBUG_TRACE("WFS must be pyrhr");
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template <class Tin, class T>
int32_t SutraCentroiderMaskedPix<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class SutraCentroiderMaskedPix<float, float>;
template class SutraCentroiderMaskedPix<uint16_t, float>;