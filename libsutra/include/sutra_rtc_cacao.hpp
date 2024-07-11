// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_rtc_cacao.hpp
//! \ingroup   libsutra
//! \class     SutraRtcCacao
//! \brief     this class provides the rtc_cacao features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef SUTRA_RTC_CACAO_H_
#define SUTRA_RTC_CACAO_H_

#include <Cacao.h>

#include <sutra_rtc.hpp>
#include <sutra_target.hpp>
#include <sutra_wfs.hpp>

template <typename Tin, typename Tcomp, typename Tout>
class SutraRtcCacao : public SutraRtc<Tin, Tcomp, Tout> {
 private:
  std::shared_ptr<ipc::Cacao<Tin>> interface_raw_frame_;
  std::shared_ptr<ipc::Cacao<float>> interface_cal_frame_;
  std::shared_ptr<ipc::Cacao<Tcomp>> interface_loop_frame_;
  std::shared_ptr<ipc::Cacao<Tout>> interface_commands_;

  std::string interface_raw_frame_name_;
  std::string interface_cal_frame_name_;
  std::string interface_loop_frame_name_;
  std::string interface_commands_name_;

  int64_t framecounter_;

  int32_t nslp_;
  int32_t ncmd_;
  int32_t nvalid_;

  bool is_initialised_;

 public:
  SutraRtcCacao(std::string interface_cal_frame_name, std::string interface_loop_frame_name);
  ~SutraRtcCacao();

  void publish();

 private:
  void allocate_buffers();
};

#endif /* SUTRA_RTC_CACAO_H_ */
