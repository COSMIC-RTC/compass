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

//! \file      sutra_rtc_cacao.h
//! \ingroup   libsutra
//! \class     sutra_rtc_cacao
//! \brief     this class provides the rtc_cacao features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#ifndef SUTRA_RTC_CACAO_H_
#define SUTRA_RTC_CACAO_H_

#include <Cacao.h>

#include <sutra_rtc.h>
#include <sutra_target.h>
#include <sutra_wfs.h>

template <typename Tin, typename Tcomp, typename Tout>
class sutra_rtc_cacao : public sutra_rtc<Tin, Tcomp, Tout> {
 private:
  std::shared_ptr<ipc::Cacao<Tin>> iRawFrame_;
  std::shared_ptr<ipc::Cacao<float>> iCalFrame_;
  std::shared_ptr<ipc::Cacao<Tcomp>> iLoopFrame_;
  std::shared_ptr<ipc::Cacao<Tout>> iCommands_;

  std::string iRawFrame_name_;
  std::string iCalFrame_name_;
  std::string iLoopFrame_name_;
  std::string iCommands_name_;

  long framecounter_;

  int nslp_;
  int ncmd_;
  int nvalid_;

  bool is_initialised_;

 public:
  sutra_rtc_cacao(std::string iCalFrame_name, std::string iLoopFrame_name);
  ~sutra_rtc_cacao();

  void publish();

 private:
  void allocateBuffers();
};

#endif /* SUTRA_RTC_CACAO_H_ */
