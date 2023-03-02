// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_rtc_cacao.cpp
//! \ingroup   libsutra
//! \class     SutraRtcCacao
//! \brief     this class provides the rtc_cacao features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#include <sutra_rtc_cacao.h>

template <typename Tin, typename Tcomp, typename Tout>
SutraRtcCacao<Tin, Tcomp, Tout>::SutraRtcCacao(std::string interface_cal_frame_name,
                                                   std::string interface_loop_frame_name)
    : SutraRtc<Tin, Tcomp, Tout>(),
      interface_cal_frame_name_(interface_cal_frame_name),
      interface_loop_frame_name_(interface_loop_frame_name),
      framecounter_(0),
      nslp_(0),
      ncmd_(0),
      nvalid_(0),
      is_initialised_(false) {
  interface_cal_frame_.reset();
  interface_loop_frame_.reset();
}

template <typename Tin, typename Tcomp, typename Tout>
SutraRtcCacao<Tin, Tcomp, Tout>::~SutraRtcCacao() {
  if (!is_initialised_) {
    return;
  }
  interface_cal_frame_.reset();
  interface_loop_frame_.reset();
}

template <typename Tin, typename Tcomp, typename Tout>
void SutraRtcCacao<Tin, Tcomp, Tout>::allocate_buffers() {
  if (is_initialised_) {
    return;
  }

  uint32_t size = this->d_centro[0]->d_img->get_dims(1);
  interface_cal_frame_ = std::make_shared<ipc::Cacao<float>>(ipc::Cacao<float>(
      interface_cal_frame_name_, std::vector<uint32_t>{size, size, 10}, -1, 1, 8, 0));

  nslp_ = 0;
  ncmd_ = 0;
  nvalid_ = 0;
  for (unsigned int i = 0; i < this->d_centro.size(); i++) {
    nvalid_ += this->d_centro[i]->d_intensities->get_nb_elements();
  }
  for (unsigned int i = 0; i < this->d_control.size(); i++) {
    nslp_ += this->d_control[i]->nslope();
    ncmd_ += this->d_control[i]->nactu();
  }

  uint32_t size_tot = nvalid_ + nslp_ + ncmd_;
  interface_loop_frame_ = std::make_shared<ipc::Cacao<Tcomp>>(ipc::Cacao<Tcomp>(
      interface_loop_frame_name_, std::vector<uint32_t>{size_tot, 1, 10}, -1, 1, 8, 3));

  strcpy(interface_loop_frame_->image()->kw[0].name, "nvalid");
  interface_loop_frame_->image()->kw[0].type = 'L';
  interface_loop_frame_->image()->kw[0].value.numl = nvalid_;
  strcpy(interface_loop_frame_->image()->kw[0].comment, "number of measurements");

  strcpy(interface_loop_frame_->image()->kw[1].name, "nslp");
  interface_loop_frame_->image()->kw[1].type = 'L';
  interface_loop_frame_->image()->kw[1].value.numl = nslp_;
  strcpy(interface_loop_frame_->image()->kw[1].comment, "number of slopes");

  strcpy(interface_loop_frame_->image()->kw[2].name, "ncmd");
  interface_loop_frame_->image()->kw[2].type = 'L';
  interface_loop_frame_->image()->kw[2].value.numl = ncmd_;
  strcpy(interface_loop_frame_->image()->kw[2].comment, "number of commands");

  is_initialised_ = true;
}

template <typename Tin, typename Tcomp, typename Tout>
void SutraRtcCacao<Tin, Tcomp, Tout>::publish() {
  if (!is_initialised_) {
    allocate_buffers();
  }

  float* iFrame = interface_cal_frame_->outputPtr();

  this->d_centro[0]->d_img->device2host(iFrame);
  interface_cal_frame_->notify();

  Tcomp* zFrame = interface_loop_frame_->outputPtr();

  for (unsigned int i = 0; i < this->d_centro.size(); i++) {
    // this->d_centro[i]->d_intensities->device2host(zFrame);
    zFrame += this->d_centro[i]->d_intensities->get_nb_elements();
  }

  for (unsigned int i = 0; i < this->d_control.size(); i++) {
    this->d_control[i]->d_centroids->device2host(zFrame);
    zFrame += this->d_control[i]->nslope();
  }

  for (unsigned int i = 0; i < this->d_control.size(); i++) {
    this->d_control[i]->d_com_clipped->device2host(zFrame);
    zFrame += this->d_control[i]->nactu();
  }

  interface_loop_frame_->notify();

  framecounter_++;
}
template class SutraRtcCacao<float, float, float>;
template class SutraRtcCacao<uint16_t, float, float>;
template class SutraRtcCacao<float, float, uint16_t>;
template class SutraRtcCacao<uint16_t, float, uint16_t>;
#ifdef CAN_DO_HALF
template class SutraRtcCacao<float, half, float>;
template class SutraRtcCacao<uint16_t, half, float>;
template class SutraRtcCacao<float, half, uint16_t>;
template class SutraRtcCacao<uint16_t, half, uint16_t>;
#endif
