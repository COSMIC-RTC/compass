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

//! \file      sutra_centroider.cpp
//! \ingroup   libsutra
//! \class     SutraCentroider
//! \brief     this class provides the centroider features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_centroider.h>

template <class Tin, class Tout>
SutraCentroider<Tin, Tout>::SutraCentroider(CarmaContext *context,
                                            SutraWfs *wfs, long nvalid,
                                            float offset, float scale,
                                            bool filter_TT, int device) {
  this->current_context = context;
  this->device = device;
  context->set_active_device(device, 1);
  this->wfs = wfs;

  this->nvalid = nvalid;
  this->offset = offset;
  this->scale = scale;
  this->nslopes = 0;
  this->npix = 0;
  this->nxsub = 0;
  this->filter_TT = filter_TT;

  long dims_data[2] = {1, this->nvalid};
  this->d_intensities = new CarmaObj<float>(current_context, dims_data);
  this->d_intensities->reset();

  this->d_centroids_ref = nullptr;
  if (this->wfs != nullptr) {
    this->d_img = this->wfs->d_binimg;
  } else {
    this->d_img = nullptr;
  }
  this->d_img_raw = nullptr;
  this->d_validx = nullptr;
  this->d_validy = nullptr;
  this->d_dark = nullptr;
  this->d_flat = nullptr;
  this->d_lutPix = nullptr;
  this->d_bincube = nullptr;
  this->d_validMask = nullptr;
  this->d_centro_filtered = nullptr;
  this->d_ref_Tip = nullptr;
  this->d_ref_Tilt = nullptr;
  this->d_TT_slopes = nullptr;
}

template <class Tin, class Tout>
SutraCentroider<Tin, Tout>::~SutraCentroider() {
  if (this->d_intensities != nullptr) delete this->d_intensities;
  if (this->d_centroids_ref != nullptr) delete this->d_centroids_ref;
  if (this->wfs == nullptr && this->d_img != nullptr) delete this->d_img;
  if (this->d_img_raw != nullptr) delete this->d_img_raw;
  if (this->d_validx != nullptr) delete this->d_validx;
  if (this->d_validy != nullptr) delete this->d_validy;
  if (this->d_dark != nullptr) delete this->d_dark;
  if (this->d_flat != nullptr) delete this->d_flat;
  if (this->d_lutPix != nullptr) delete this->d_lutPix;
  if (this->d_bincube != nullptr) delete this->d_bincube;
  if (this->d_validMask != nullptr) delete this->d_validMask;
  if (this->d_TT_slopes != nullptr) delete this->d_TT_slopes;
  if (this->d_centro_filtered != nullptr) delete this->d_centro_filtered;
  if (this->d_ref_Tip != nullptr) delete this->d_ref_Tip;
  if (this->d_ref_Tilt != nullptr) delete this->d_ref_Tilt;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::set_scale(float scale) {
  this->scale = scale;
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::set_offset(float offset) {
  this->offset = offset;
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::set_nxsub(int nxsub) {
  this->nxsub = nxsub;
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::init_calib(int n, int m) {
  current_context->set_active_device(device, 1);
  long dims_data2[3] = {2, n, m};
  if (this->d_img == nullptr) {
    this->d_img = new CarmaObj<float>(current_context, dims_data2);
  }
  if (this->d_dark == nullptr) {
    this->d_dark = new CarmaObj<float>(current_context, dims_data2);
    this->d_dark->reset();
  }
  if (this->d_flat == nullptr) {
    this->d_flat = new CarmaObj<float>(current_context, dims_data2);
    this->d_flat->memset(1.f);
  }
  if (this->d_lutPix == nullptr) {
    long dims_data1[3] = {1, n * m};
    this->d_lutPix = new CarmaObj<int>(current_context, dims_data1);
    std::vector<int> h_lutPix(n * m);
    for (int i = 0; i < n * m; ++i) h_lutPix[i] = i;
    this->d_lutPix->host2device(h_lutPix.data());
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::init_roi(int N) {
  current_context->set_active_device(device, 1);
  if (this->d_validx == nullptr) {
    long dims_data[2] = {1, N};
    this->d_validx = new CarmaObj<int>(current_context, dims_data);
    this->d_validy = new CarmaObj<int>(current_context, dims_data);
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::set_dark(float *dark, int n) {
  current_context->set_active_device(device, 1);
  if (this->d_dark == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_dark = new CarmaObj<float>(current_context, dims_data2);
  }
  this->d_dark->host2device(dark);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::set_flat(float *flat, int n) {
  current_context->set_active_device(device, 1);
  if (this->d_flat == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_flat = new CarmaObj<float>(current_context, dims_data2);
  }
  this->d_flat->host2device(flat);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::set_lutPix(int *lutPix, int n) {
  current_context->set_active_device(device, 1);
  if (this->d_lutPix == nullptr) {
    long dims_data1[2] = {1, n};
    this->d_lutPix = new CarmaObj<int>(current_context, dims_data1);
  }
  this->d_lutPix->host2device(lutPix);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::get_validMask() {
  this->current_context->set_active_device(this->device, 1);
  if (this->d_validMask == nullptr) {
    if (this->d_img == nullptr) {
      std::cout << "RTC image has not been initialized" << std::endl;
      return EXIT_FAILURE;
    }
    this->d_validMask = new CarmaObj<int>(current_context, d_img->get_dims());
    this->d_validMask->reset();
  }

  fill_validMask(this->d_validMask->get_dims(1), this->npix, this->nvalid,
                 this->d_validMask->get_data(), this->d_validx->get_data(),
                 this->d_validy->get_data(),
                 current_context->get_device(device));

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::calibrate_img_validPix(cudaStream_t stream) {
  current_context->set_active_device(device, 1);

  if (this->d_img_raw == nullptr) {
    std::cout << "Image not initialized\n" << std::endl;
    return EXIT_FAILURE;
  }

  const long *dims = this->d_img_raw->get_dims();
  init_calib(dims[1], dims[2]);

  this->d_img->reset();
  string type = get_type();

  if ((type == "maskedpix") || (type == "pyrhr")) {
    return calibration_validPix_pyr<Tin>(
        this->d_img_raw->get_data(), this->d_img->get_data(),
        this->d_dark->get_data(), this->d_flat->get_data(),
        this->d_lutPix->get_data(), this->d_validx->get_data(),
        this->d_validy->get_data(), this->d_validy->get_nb_elements(), dims[1],
        this->current_context->get_device(this->device), stream);
  } else {
    return calibration_validPix_sh<Tin>(
        this->npix, this->d_img->get_dims(1), this->nvalid,
        this->d_img_raw->get_data(), this->d_img->get_data(),
        this->d_dark->get_data(), this->d_flat->get_data(),
        this->d_lutPix->get_data(), this->d_validx->get_data(),
        this->d_validy->get_data(),
        this->current_context->get_device(this->device), stream);
  }
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::calibrate_img(cudaStream_t stream) {
  current_context->set_active_device(device, 1);

  if (this->d_img_raw == nullptr) {
    std::cout << "Image not initialized\n" << std::endl;
    return EXIT_FAILURE;
  }

  const long *dims = this->d_img_raw->get_dims();
  init_calib(dims[1], dims[2]);

  calibration<Tin>(this->d_img_raw->get_data(), this->d_img->get_data(),
                   this->d_dark->get_data(), this->d_flat->get_data(),
                   this->d_lutPix->get_data(), this->d_img->get_nb_elements(),
                   this->current_context->get_device(this->device), stream);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::load_img(CarmaObj<Tin> *img) {
  return this->load_img(img->get_data(), img->get_dims(1), img->get_device());
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::load_img(Tin *img, int n) {
  return this->load_img(img, n, -1);
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::load_img(Tin *img, int n, int location) {
  return this->load_img(img, n, n, location);
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::load_img(Tin *img, int m, int n, int location) {
  init_img_raw(m, n);
  if (location < 0) {  // img data on host
    this->d_img_raw->host2device(img);
  } else {  // img data on device
    this->d_img_raw->copy_from(img, this->d_img_raw->get_nb_elements());
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::init_img_raw(int m, int n) {
  current_context->set_active_device(device, 1);
  if (this->d_img_raw == nullptr) {
    long dims_data2[3] = {2, m, n};
    this->d_img_raw = new CarmaObj<Tin>(current_context, dims_data2);
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::set_npix(int npix) {
  this->npix = npix;

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::load_validpos(int *ivalid, int *jvalid, int N) {
  current_context->set_active_device(device, 1);
  if (this->d_validx == nullptr) {
    this->init_roi(N);
  }

  this->d_validx->host2device(ivalid);
  this->d_validy->host2device(jvalid);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::set_centroids_ref(float *centroids_ref) {
  this->d_centroids_ref->host2device(centroids_ref);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::init_TT_filter() {
  this->current_context->set_active_device(device, 1);
  long dims_data[2] = {1, 2};
  this->d_TT_slopes = new CarmaObj<float>(this->current_context, dims_data);
  dims_data[1] = this->nslopes;
  this->d_centro_filtered =
      new CarmaObj<float>(this->current_context, dims_data);
  this->d_ref_Tip = new CarmaObj<float>(this->current_context, dims_data);
  this->d_ref_Tilt = new CarmaObj<float>(this->current_context, dims_data);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::apply_TT_filter(Tout *centroids) {
  return this->apply_TT_filter_impl(centroids, std::is_same<Tout, float>());
}

template <class Tin, class Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
SutraCentroider<Tin, Tout>::apply_TT_filter_impl(Tout *centroids,
                                                 std::true_type) {
  this->d_centro_filtered->copy_from(centroids, this->nslopes);

  float tip = this->d_centro_filtered->dot(this->d_ref_Tip, 1, 1);
  float tilt = this->d_centro_filtered->dot(this->d_ref_Tilt, 1, 1);

  this->d_centro_filtered->axpy(-1.f * tip, this->d_ref_Tip, 1, 1);
  this->d_centro_filtered->axpy(-1.f * tilt, this->d_ref_Tilt, 1, 1);

  this->d_centro_filtered->copy_into(centroids, this->nslopes);

  float TT_data[2] = {tip, tilt};
  this->d_TT_slopes->host2device(TT_data);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int SutraCentroider<Tin, Tout>::apply_TT_filter_impl(Tout *centroids,
                                                     std::false_type) {
  DEBUG_TRACE("Tip/tilt filtering is only implemented in single precision");
  return EXIT_SUCCESS;
}

template class SutraCentroider<float, float>;
template class SutraCentroider<uint16_t, float>;
#ifdef CAN_DO_HALF
template class SutraCentroider<float, half>;
template class SutraCentroider<uint16_t, half>;
#endif
