// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the
//  terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for
//  the simulation of AO systems.
//
//  The final product includes a software package for simulating all the
//  critical subcomponents of AO, particularly in the context of the ELT and a
//  real-time core based on several control approaches, with performances
//  consistent with its integration into an instrument. Taking advantage of the
//  specific hardware architecture of the GPU, the COMPASS tool allows to
//  achieve adequate execution speeds to conduct large simulation campaigns
//  called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to
//  both testspecific components of AO of the E-ELT (such as wavefront analysis
//  device with a pyramid or elongated Laser star), and various systems
//  configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//  details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with COMPASS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_centroider_maskedPix.cpp
//! \ingroup   libsutra
//! \class     sutra_centroider_maskedPix
//! \brief     this class provides the centroider_maskedPix features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_centroider_maskedPix.h>
#include <string>

template <class Tin, class T>
sutra_centroider_maskedPix<Tin, T>::sutra_centroider_maskedPix(
    carma_context *context, sutra_wfs *wfs, long nvalid, long npupils,
    float offset, float scale, bool filter_TT, int device)
    : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT,
                               device) {
  context->set_activeDevice(device, 1);

  if (wfs != nullptr)
    this->nslopes = wfs->d_validsubsx->getNbElem();
  else
    this->nslopes = nvalid * npupils;
  long dims_data[2] = {1, this->nslopes};
  if (this->d_intensities != nullptr) delete this->d_intensities;
  this->d_intensities = new carma_obj<float>(context, dims_data);
  this->d_intensities->init_reduceCub();
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();
  this->d_selected_pix = nullptr;
  this->d_mask = nullptr;

  if (filter_TT) {
    this->init_TT_filter();
  }
}

template <class Tin, class T>
sutra_centroider_maskedPix<Tin, T>::~sutra_centroider_maskedPix() {}

template <class Tin, class T>
string sutra_centroider_maskedPix<Tin, T>::get_type() {
  return "maskedpix";  // TODO: Fix an use constants !!!
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::fill_selected_pix(carma_obj<T> *pix) {
  this->current_context->set_activeDevice(this->device, 1);

  const long *dims_data2;
  if (this->d_img != nullptr) {
    dims_data2 = this->d_img->getDims();
  } else if (this->wfs != nullptr) {
    dims_data2 = this->wfs->d_binimg->getDims();
  } else {
    std::cerr << "Image not initialized" << std::endl;
    return EXIT_FAILURE;
  }
  if (this->d_selected_pix == nullptr) {
    this->d_selected_pix = new carma_obj<T>(this->current_context, dims_data2);
  }
  this->d_selected_pix->reset();
  pyr_fill_selected_pix(this->d_selected_pix->getData(), dims_data2[1],
                        pix->getData(), this->d_validx->getData(),
                        this->d_validy->getData(), this->d_validx->getNbElem(),
                        this->current_context->get_device(this->device));
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::fill_mask() {
  this->current_context->set_activeDevice(this->device, 1);

  const long *dims_data2;
  if (this->d_img != nullptr) {
    dims_data2 = this->d_img->getDims();
  } else if (this->wfs != nullptr) {
    dims_data2 = this->wfs->d_binimg->getDims();
  } else {
    std::cerr << "Image not initialized" << std::endl;
    return EXIT_FAILURE;
  }
  if (this->d_mask == nullptr) {
    this->d_mask = new carma_obj<T>(this->current_context, dims_data2);
  }
  this->d_mask->reset();
  pyr_fill_mask(this->d_mask->getData(), dims_data2[1],
                        this->d_validx->getData(),
                        this->d_validy->getData(), this->d_validx->getNbElem(),
                        this->current_context->get_device(this->device));
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::get_cog(float *img, float *intensities,
                                                T *centroids, int nvalid,
                                                int npix, int ntot, cudaStream_t stream) {
  // TODO(Implement sutra_centroider_maskedPix<Tin, T>::get_cog)

  return get_maskedPix(img, intensities, centroids, this->d_validx->getData(),
                       this->d_validy->getData(), this->nvalid, ntot);
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::get_maskedPix(
    float *img, float *intensities, T *centroids, int *subindx, int *subindy,
    int nvalid, int ns) {
  this->current_context->set_activeDevice(this->device, 1);

  fill_intensities(this->d_intensities->getData(), img, subindx, subindy, ns,
                   this->nslopes,
                   this->current_context->get_device(this->device));

  // T p_sum = reduce<T>(this->d_intensities->getData(), this->nslopes);
  this->d_intensities->reduceCub();

  getMaskedPix<T>(centroids, this->d_centroids_ref->getData(), img, subindx,
                  subindy, this->d_intensities->getOData(), ns, this->nslopes,
                  this->current_context->get_device(this->device));

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::get_cog(float *intensities, T *slopes,
                                                bool noise) {
  if (this->wfs != nullptr) {
    if (this->wfs->type == "pyrhr") {
      if (noise || this->wfs->roket == false) {
        return this->get_maskedPix(
            *(this->wfs->d_binimg), intensities, slopes,
            *(this->wfs->d_validsubsx), *(this->wfs->d_validsubsy),
            this->wfs->nvalid, this->wfs->nfft / this->wfs->nrebin);
      } else
        return this->get_maskedPix(
            *(this->wfs->d_binimg_notnoisy), intensities, slopes,
            *(this->wfs->d_validsubsx), *(this->wfs->d_validsubsy),
            this->wfs->nvalid, this->wfs->nfft / this->wfs->nrebin);
    } else
      DEBUG_TRACE("WFS must be pyrhr");
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class sutra_centroider_maskedPix<float, float>;
template class sutra_centroider_maskedPix<uint16_t, float>;
#ifdef CAN_DO_HALF
template <>
int sutra_centroider_maskedPix<float, half>::get_cog(float *intensities,
                                                     half *slopes, bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_maskedPix<uint16_t, half>::get_cog(float *intensities,
                                                        half *slopes,
                                                        bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_maskedPix<float, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_maskedPix<uint16_t, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}
template class sutra_centroider_maskedPix<float, half>;
template class sutra_centroider_maskedPix<uint16_t, half>;
#endif
