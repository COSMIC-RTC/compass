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

//! \file      sutra_rtc_cacao.cpp
//! \ingroup   libsutra
//! \class     sutra_rtc_cacao
//! \brief     this class provides the rtc_cacao features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_rtc_cacao.h>

template <typename Tin, typename Tcomp, typename Tout>
sutra_rtc_cacao<Tin, Tcomp, Tout>::sutra_rtc_cacao(std::string iCalFrame_name,
                                                   std::string iLoopFrame_name)
    : sutra_rtc<Tin, Tcomp, Tout>(),
      iCalFrame_name_(iCalFrame_name),
      iLoopFrame_name_(iLoopFrame_name),
      framecounter_(0),
      nslp_(0),
      ncmd_(0),
      nvalid_(0),
      is_initialised_(false) {
  iCalFrame_.reset();
  iLoopFrame_.reset();
}

template <typename Tin, typename Tcomp, typename Tout>
sutra_rtc_cacao<Tin, Tcomp, Tout>::~sutra_rtc_cacao() {
  if (!is_initialised_) {
    return;
  }
  iCalFrame_.reset();
  iLoopFrame_.reset();
}

template <typename Tin, typename Tcomp, typename Tout>
void sutra_rtc_cacao<Tin, Tcomp, Tout>::allocateBuffers() {
  if (is_initialised_) {
    return;
  }

  uint32_t size = this->d_centro[0]->d_img->getDims(1);
  iCalFrame_ = std::make_shared<ipc::Cacao<float>>(ipc::Cacao<float>(
      iCalFrame_name_, std::vector<uint32_t>{size, size, 10}, -1, 1, 8, 0));

  nslp_ = 0;
  ncmd_ = 0;
  nvalid_ = 0;
  for (unsigned int i = 0; i < this->d_centro.size(); i++) {
    nvalid_ += this->d_centro[i]->d_intensities->getNbElem();
  }
  for (unsigned int i = 0; i < this->d_control.size(); i++) {
    nslp_ += this->d_control[i]->nslope();
    ncmd_ += this->d_control[i]->nactu();
  }

  uint32_t size_tot = nvalid_ + nslp_ + ncmd_;
  iLoopFrame_ = std::make_shared<ipc::Cacao<Tcomp>>(ipc::Cacao<Tcomp>(
      iLoopFrame_name_, std::vector<uint32_t>{size_tot, 1, 10}, -1, 1, 8, 0));

  is_initialised_ = true;
}

template <typename Tin, typename Tcomp, typename Tout>
void sutra_rtc_cacao<Tin, Tcomp, Tout>::publish() {
  if (!is_initialised_) {
    allocateBuffers();
  }

  float* iFrame = iCalFrame_->outputPtr();

  this->d_centro[0]->d_img->device2host(iFrame);
  iCalFrame_->notify();

  Tcomp* zFrame = iLoopFrame_->outputPtr();

  for (unsigned int i = 0; i < this->d_centro.size(); i++) {
    // this->d_centro[i]->d_intensities->device2host(zFrame);
    zFrame += this->d_centro[i]->d_intensities->getNbElem();
  }

  for (unsigned int i = 0; i < this->d_control.size(); i++) {
    this->d_control[i]->d_centroids->device2host(zFrame);
    zFrame += this->d_control[i]->nslope();
  }

  for (unsigned int i = 0; i < this->d_control.size(); i++) {
    this->d_control[i]->d_comClipped->device2host(zFrame);
    zFrame += this->d_control[i]->nactu();
  }

  iLoopFrame_->notify();

  framecounter_++;
}
template class sutra_rtc_cacao<float, float, float>;
template class sutra_rtc_cacao<uint16_t, float, float>;
template class sutra_rtc_cacao<float, float, uint16_t>;
template class sutra_rtc_cacao<uint16_t, float, uint16_t>;
#ifdef CAN_DO_HALF
template class sutra_rtc_cacao<float, half, float>;
template class sutra_rtc_cacao<uint16_t, half, float>;
template class sutra_rtc_cacao<float, half, uint16_t>;
template class sutra_rtc_cacao<uint16_t, half, uint16_t>;
#endif
