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

//! \file      sutra_centroider_tcog.cpp
//! \ingroup   libsutra
//! \class     SutraCentroiderTcog
//! \brief     this class provides the centroider_tcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
                this->current_context->get_device(this->device));

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }

  // TODO: Implement get_cog_async
  // #ifndef USE_OLD
  //   // Modif Nono !!!
  //   // Now subap_reduce modify cube to apply the threshold on the cube so
  //   // get_centroids should apply threshold a second time
  //   subap_reduce_new(ntot, npix * npix, nvalid, cube, intensities,
  //                    this->threshold,
  //                    this->current_context->get_device(device));

  //   get_centroids(ntot, npix * npix, nvalid, npix, cube, centroids,
  //   intensities,
  //                 this->threshold, this->scale, this->offset,
  //                 this->current_context->get_device(device));
  // #else
  //   subap_reduce(ntot, npix * npix, nvalid, cube, intensities,
  //   this->threshold,
  //                this->current_context->get_device(device));

  //   get_centroids(ntot, npix * npix, nvalid, npix, cube, centroids,
  //   intensities,
  //                 this->threshold, this->scale, this->offset,
  //                 this->current_context->get_device(device));
  // #endif
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
