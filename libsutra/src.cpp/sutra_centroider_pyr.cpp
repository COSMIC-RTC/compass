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

//! \file      sutra_centroider_pyr.cpp
//! \ingroup   libsutra
//! \class     sutra_centroider_pyr
//! \brief     this class provides the centroider_pyr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_centroider_pyr.h>
#include <iostream>
#include <string>

template <class Tin, class T>
sutra_centroider_pyr<Tin, T>::sutra_centroider_pyr(carma_context *context,
                                                   sutra_wfs *wfs, long nvalid,
                                                   float offset, float scale,
                                                   bool filter_TT, int device)
    : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT,
                               device) {
  context->set_activeDevice(device, 1);

  this->nslopes = 2 * nvalid;

  this->pyr_type = "pyrhr";

  this->valid_thresh = 1e-4;

  // centroider method by default nosin_global
  this->method = Method_CoG(false, false);

  this->d_intensities->init_reduceCub();
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();

  if (filter_TT) {
    this->init_TT_filter();
  }
}

template <class Tin, class T>
sutra_centroider_pyr<Tin, T>::~sutra_centroider_pyr() {}

template <class Tin, class T>
string sutra_centroider_pyr<Tin, T>::get_type() {
  return this->pyr_type;
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::set_valid_thresh(T valid_thresh) {
  this->valid_thresh = valid_thresh;
  return EXIT_SUCCESS;
}
template <class Tin, class T>
T sutra_centroider_pyr<Tin, T>::get_valid_thresh() {
  return this->valid_thresh;
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::set_method(Method_CoG method) {
  this->method = method;
  return EXIT_SUCCESS;
}

template <class Tin, class T>
Method_CoG sutra_centroider_pyr<Tin, T>::get_method() {
  return this->method;
}

template <class Tin, class T>
string sutra_centroider_pyr<Tin, T>::get_method_str() {
  return Method_CoG::str(this->method);
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::get_cog(float *cube, float *intensities,
                                          T *centroids, int nvalid, int npix,
                                          int ntot) {
  // TODO(Implement sutra_centroider_pyr<Tin, T>::get_cog)

  return get_pyr(cube, intensities, centroids, this->d_validx->getData(),
                 this->d_validy->getData(), this->nvalid,
                 this->d_img->getDims(1), 4);
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::get_pyr(float *cube, float *intensities,
                                          T *centroids, int *subindx,
                                          int *subindy, int nvalid, int ns,
                                          int nim) {
  this->current_context->set_activeDevice(this->device, 1);

  pyr_intensities(this->d_intensities->getData(), cube, subindx, subindy, ns,
                  nvalid, nim, this->current_context->get_device(this->device));

  if (!(this->method.isLocal)) {
    // float p_sum = reduce<float>(this->d_intensities->getData(), nvalid);
    this->d_intensities->reduceCub();
    fillvalues<float>(this->d_intensities->getData(),
                      this->d_intensities->getOData(), nvalid,
                      this->current_context->get_device(this->device));
  }

  pyr2_slopes(centroids, this->d_centroids_ref->getData(), cube, subindx,
              subindy, this->d_intensities->getData(), ns, nvalid, this->scale,
              this->valid_thresh,
              this->method.isSinus,  // if we are using a sin method
              this->current_context->get_device(this->device));

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::get_cog(float *intensities, T *slopes,
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
int sutra_centroider_pyr<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class sutra_centroider_pyr<float, float>;
template class sutra_centroider_pyr<uint16_t, float>;
