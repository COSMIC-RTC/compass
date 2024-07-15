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

//! \file      sutra_centroider_pbcog.hpp
//! \ingroup   libsutra
//! \class     sutra_centroider_pbcog
//! \brief     this class provides the centroider_pbcog features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_centroider_bpcog.hpp>

template <class Tin, class T>
SutraCentroiderBpcog<Tin, T>::SutraCentroiderBpcog(
    CarmaContext *context, SutraWfs *wfs, int64_t nvalid, float offset,
    float scale, bool filter_TT, int32_t device, int32_t nmax)
    : SutraCentroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT,
                               device) {
  this->nslopes = 2 * nvalid;
  this->nmax = nmax;
  int64_t dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new CarmaObj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();

  if (filter_TT) {
    this->init_TT_filter();
  }

  int64_t dims_data[3];
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
int32_t SutraCentroiderBpcog<Tin, T>::set_nmax(int32_t nmax) {
  this->current_context->set_active_device(this->device, 1);
  if(nmax != this->nmax) {
    this->nmax = nmax;
    delete this->d_bpix;
    delete this->d_bpind;

    int64_t dims_data[3];
    dims_data[0] = 2;
    dims_data[1] = nmax;
    dims_data[2] = this->nvalid;

    this->d_bpix = new CarmaObj<T>(this->current_context, dims_data);
    this->d_bpind = new CarmaObj<uint>(this->current_context, dims_data);
  }
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderBpcog<Tin, T>::get_cog(float *img, float *intensities,
                                            T *centroids, int32_t nvalid, int32_t npix,
                                            int32_t ntot, cudaStream_t stream) {
  this->current_context->set_active_device(this->device, 1);

  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->get_data(), this->d_validx->get_data(),
                this->d_validy->get_data(), intensities, this->nmax, this->scale,
                this->offset, this->slope_order, this->current_context->get_device(this->device), stream);

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int32_t SutraCentroiderBpcog<Tin, T>::get_cog(float *intensities, T *slopes,
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
int32_t SutraCentroiderBpcog<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class SutraCentroiderBpcog<float, float>;
template class SutraCentroiderBpcog<uint16_t, float>;
