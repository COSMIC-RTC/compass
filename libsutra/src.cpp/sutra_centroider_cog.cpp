// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_centroider_cog.cpp
//! \ingroup   libsutra
//! \class     sutra_centroider_cog
//! \brief     this class provides the centroider_cog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#include <sutra_centroider_cog.h>
#include <string>

template <class Tin, class T>
sutra_centroider_cog<Tin, T>::sutra_centroider_cog(carma_context *context,
                                                   sutra_wfs *wfs, long nvalid,
                                                   float offset, float scale,
                                                   bool filter_TT, int device)
    : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT,
                               device) {
  this->nslopes = 2 * nvalid;
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();
  if (filter_TT) {
    this->init_TT_filter();
  }
}

template <class Tin, class T>
sutra_centroider_cog<Tin, T>::~sutra_centroider_cog() {}

template <class Tin, class T>
string sutra_centroider_cog<Tin, T>::get_type() {
  return "cog";
}

template <class Tin, class T>
int sutra_centroider_cog<Tin, T>::get_cog(float *img, float *intensities,
                                          T *centroids, int nvalid, int npix,
                                          int ntot) {
  this->current_context->set_activeDevice(this->device, 1);

  // subap_reduce(ntot, (npix * npix), nvalid, img, ref,
  //              this->current_context->get_device(device));
  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->getData(), this->d_validx->getData(),
                this->d_validy->getData(), intensities, this->scale,
                this->offset, this->current_context->get_device(this->device));

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_cog<Tin, T>::get_cog(float *intensities, T *slopes,
                                          bool noise) {
  if (this->wfs != nullptr) {
    if (noise || this->wfs->roket == false) {
      return this->get_cog(*(this->wfs->d_binimg), intensities, slopes,
                           this->wfs->nvalid_tot, this->wfs->npix,
                           this->wfs->d_binimg->getDims()[1]);
    } else {
      return this->get_cog(*(this->wfs->d_binimg_notnoisy), intensities, slopes,
                           this->wfs->nvalid_tot, this->wfs->npix,
                           this->wfs->d_binimg->getDims()[1]);
    }
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_cog<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_SUCCESS;
}

template class sutra_centroider_cog<float, float>;
template class sutra_centroider_cog<uint16_t, float>;

#ifdef CAN_DO_HALF
template <>
int sutra_centroider_cog<float, half>::get_cog(float *intensities, half *slopes,
                                               bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_cog<uint16_t, half>::get_cog(float *intensities,
                                                  half *slopes, bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_cog<float, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_cog<uint16_t, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template class sutra_centroider_cog<float, half>;
template class sutra_centroider_cog<uint16_t, half>;
#endif
