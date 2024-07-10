// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_cured.cpp
//! \ingroup   libsutra
//! \class     SutraControllerCured
//! \brief     this class provides the controller_cured features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <cured.h>
#include <sutra_controller_cured.hpp>
#include <string>

template <typename T, typename Tout>
SutraControllerCured<T, Tout>::SutraControllerCured(
    CarmaContext *context, int64_t nslopes, int64_t nactu, float delay,
    SutraDms *dms, int32_t *idx_dms, int32_t ndm, int32_t *idx_centro, int32_t ncentro)
    : SutraController<T, Tout>(context, nslopes, nactu, delay, dms,
                                idx_dms, ndm, idx_centro, ncentro),
      ndivs(0),
      tt_flag(false),
      h_syscure(nullptr),
      h_parcure(nullptr) {
  int64_t dims_data1[2] = {1, 0};
  int64_t dims_data2[3] = {2, 0, 0};

  dims_data1[1] = nslopes;
  this->d_centroids = new CarmaObj<T>(context, dims_data1);

  dims_data1[1] = nslopes;
  this->h_centroids = new CarmaHostObj<T>(dims_data1, MA_PAGELOCK);

  dims_data1[1] = nactu;
  this->h_err = new CarmaHostObj<T>(dims_data1, MA_PAGELOCK);
  this->d_err = new CarmaObj<T>(context, dims_data1);

  dims_data2[1] = nslopes;
  dims_data2[2] = nactu;

  this->d_imat = new CarmaObj<T>(context, dims_data2);
  if ((int32_t)delay > 0) {
    dims_data2[1] = nslopes;
    dims_data2[2] = (int32_t)delay + 1;
    this->d_cenbuff = new CarmaObj<T>(context, dims_data2);
  } else
    this->d_cenbuff = 0L;
}

template <typename T, typename Tout>
SutraControllerCured<T, Tout>::~SutraControllerCured() {
  this->current_context->set_active_device(this->device, 1);

  if (this->h_centroids != nullptr) delete this->h_centroids;
  if (this->h_err != nullptr) delete this->h_err;
  if (this->d_imat != nullptr) delete this->d_imat;
  if (this->d_err != nullptr) delete this->d_err;
  if (this->d_cenbuff) delete this->d_cenbuff;
  curefree((sysCure *)this->h_syscure, (parCure *)this->h_parcure);
}

template <typename T, typename Tout>
int32_t SutraControllerCured<T, Tout>::comp_com() {
  this->current_context->set_active_device(this->device, 1);

  // this->frame_delay();
  h_centroids->cpy_obj(this->d_centroids, cudaMemcpyDeviceToHost);

  if (this->tt_flag) {
    cured((sysCure *)this->h_syscure, (parCure *)this->h_parcure,
          this->h_centroids->get_data(), this->h_err->get_data(),
          this->h_err->get_data_at(this->h_err->get_nb_elements() - 2),
          this->h_err->get_data_at(this->h_err->get_nb_elements() - 1));

    //*this->h_err->get_data_at(this->h_err->get_nb_elements()-2) *= -1.0f;
    //*this->h_err->get_data_at(this->h_err->get_nb_elements()-1) *= -1.0f;
  } else {
    cured((sysCure *)this->h_syscure, (parCure *)this->h_parcure,
          this->h_centroids->get_data(), this->h_err->get_data());
  }

  h_err->cpy_obj(this->d_err, cudaMemcpyHostToDevice);

  mult_int(this->d_com->get_data(), this->d_err->get_data(), -1.0f * this->gain,
           this->nactu(), this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int32_t SutraControllerCured<T, Tout>::init_cured(int32_t nxsubs, int32_t *isvalid,
                                                int32_t ndivs, int32_t tt) {
  if (tt > 0)
    this->tt_flag = true;
  else
    this->tt_flag = false;
  this->ndivs = (ndivs > 0) ? ndivs : 1;
  this->h_syscure = (void *)cureSystem(
      nxsubs, this->nslope() / 2.,
      this->tt_flag ? this->nactu() - 2 : this->nactu(), isvalid, this->ndivs);
  this->h_parcure = (void *)cureInit((sysCure *)this->h_syscure);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int32_t SutraControllerCured<T, Tout>::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  if ((int32_t)this->delay > 0) {
    for (int32_t cc = 0; cc < this->delay; cc++)
      shift_buf(&((this->d_cenbuff->get_data())[cc * this->nslope()]), 1,
                this->nslope(),
                this->current_context->get_device(this->device));

    carma_safe_call(cudaMemcpy(
        this->d_cenbuff->get_data_at((int32_t)this->delay * this->nslope()),
        this->d_centroids->get_data(), sizeof(T) * this->nslope(),
        cudaMemcpyDeviceToDevice));

    carma_safe_call(
        cudaMemcpy(this->d_centroids->get_data(), this->d_cenbuff->get_data(),
                   sizeof(T) * this->nslope(), cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

template class SutraControllerCured<float, float>;
template class SutraControllerCured<float, uint16_t>;
