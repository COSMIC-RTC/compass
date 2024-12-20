// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_centroider.cpp
//! \ingroup   libsutra
//! \class     SutraCentroider
//! \brief     this class provides the centroider features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_centroider.hpp>

template <class Tin, class Tout>
SutraCentroider<Tin, Tout>::SutraCentroider(CarmaContext *context,
                                            SutraWfs *wfs, int64_t nvalid,
                                            float offset, float scale,
                                            bool filter_TT, int32_t device) {
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

  int64_t dims_data[2] = {1, this->nvalid};
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
int32_t SutraCentroider<Tin, Tout>::set_scale(float scale) {
  this->scale = scale;
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::set_offset(float offset) {
  this->offset = offset;
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::set_nxsub(int32_t nxsub) {
  this->nxsub = nxsub;
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::init_calib(int32_t n, int32_t m) {
  current_context->set_active_device(device, 1);
  int64_t dims_data2[3] = {2, n, m};
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
    int64_t dims_data1[3] = {1, n * m};
    this->d_lutPix = new CarmaObj<int32_t>(current_context, dims_data1);
    std::vector<int32_t> h_lutPix(n * m);
    for (int32_t i = 0; i < n * m; ++i) h_lutPix[i] = i;
    this->d_lutPix->host2device(h_lutPix.data());
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::init_roi(int32_t N) {
  current_context->set_active_device(device, 1);
  if (this->d_validx == nullptr) {
    int64_t dims_data[2] = {1, N};
    this->d_validx = new CarmaObj<int32_t>(current_context, dims_data);
    this->d_validy = new CarmaObj<int32_t>(current_context, dims_data);
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::set_dark(float *dark, int32_t n) {
  current_context->set_active_device(device, 1);
  if (this->d_dark == nullptr) {
    int64_t dims_data2[3] = {2, n, n};
    this->d_dark = new CarmaObj<float>(current_context, dims_data2);
  }
  this->d_dark->host2device(dark);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::set_flat(float *flat, int32_t n) {
  current_context->set_active_device(device, 1);
  if (this->d_flat == nullptr) {
    int64_t dims_data2[3] = {2, n, n};
    this->d_flat = new CarmaObj<float>(current_context, dims_data2);
  }
  this->d_flat->host2device(flat);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::set_lutPix(int32_t *lutPix, int32_t n) {
  current_context->set_active_device(device, 1);
  if (this->d_lutPix == nullptr) {
    int64_t dims_data1[2] = {1, n};
    this->d_lutPix = new CarmaObj<int32_t>(current_context, dims_data1);
  }
  this->d_lutPix->host2device(lutPix);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::get_validMask() {
  this->current_context->set_active_device(this->device, 1);
  if (this->d_validMask == nullptr) {
    if (this->d_img == nullptr) {
      std::cout << "RTC image has not been initialized" << std::endl;
      return EXIT_FAILURE;
    }
    this->d_validMask = new CarmaObj<int32_t>(current_context, d_img->get_dims());
    this->d_validMask->reset();
  }

  fill_validMask(this->d_validMask->get_dims(1), this->npix, this->nvalid,
                 this->d_validMask->get_data(), this->d_validx->get_data(),
                 this->d_validy->get_data(),
                 current_context->get_device(device));

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::calibrate_img_validPix(cudaStream_t stream) {
  current_context->set_active_device(device, 1);

  if (this->d_img_raw == nullptr) {
    std::cout << "Image not initialized\n" << std::endl;
    return EXIT_FAILURE;
  }

  const int64_t *dims = this->d_img_raw->get_dims();
  init_calib(dims[1], dims[2]);

  this->d_img->reset(stream);
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
int32_t SutraCentroider<Tin, Tout>::calibrate_img(cudaStream_t stream) {
  current_context->set_active_device(device, 1);

  if (this->d_img_raw == nullptr) {
    std::cout << "Image not initialized\n" << std::endl;
    return EXIT_FAILURE;
  }

  const int64_t *dims = this->d_img_raw->get_dims();
  init_calib(dims[1], dims[2]);

  calibration<Tin>(this->d_img_raw->get_data(), this->d_img->get_data(),
                   this->d_dark->get_data(), this->d_flat->get_data(),
                   this->d_lutPix->get_data(), this->d_img->get_nb_elements(),
                   this->current_context->get_device(this->device), stream);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::load_img(CarmaObj<Tin> *img) {
  return this->load_img(img->get_data(), img->get_dims(1), img->get_device());
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::load_img(Tin *img, int32_t n) {
  return this->load_img(img, n, -1);
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::load_img(Tin *img, int32_t n, int32_t location) {
  return this->load_img(img, n, n, location);
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::load_img(Tin *img, int32_t m, int32_t n, int32_t location) {
  init_img_raw(m, n);
  if (location < 0) {  // img data on host
    this->d_img_raw->host2device(img);
  } else {  // img data on device
    this->d_img_raw->copy_from(img, this->d_img_raw->get_nb_elements());
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::init_img_raw(int32_t m, int32_t n) {
  current_context->set_active_device(device, 1);
  if (this->d_img_raw == nullptr) {
    int64_t dims_data2[3] = {2, m, n};
    this->d_img_raw = new CarmaObj<Tin>(current_context, dims_data2);
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::set_npix(int32_t npix) {
  this->npix = npix;

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::load_validpos(int32_t *ivalid, int32_t *jvalid, int32_t N) {
  current_context->set_active_device(device, 1);
  if (this->d_validx == nullptr) {
    this->init_roi(N);
  }

  this->d_validx->host2device(ivalid);
  this->d_validy->host2device(jvalid);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::set_centroids_ref(float *centroids_ref) {
  this->d_centroids_ref->host2device(centroids_ref);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::init_TT_filter() {
  this->current_context->set_active_device(device, 1);
  int64_t dims_data[2] = {1, 2};
  this->d_TT_slopes = new CarmaObj<float>(this->current_context, dims_data);
  dims_data[1] = this->nslopes;
  this->d_centro_filtered =
      new CarmaObj<float>(this->current_context, dims_data);
  this->d_ref_Tip = new CarmaObj<float>(this->current_context, dims_data);
  this->d_ref_Tilt = new CarmaObj<float>(this->current_context, dims_data);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int32_t SutraCentroider<Tin, Tout>::apply_TT_filter(Tout *centroids) {
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

template class SutraCentroider<float, float>;
template class SutraCentroider<uint16_t, float>;
