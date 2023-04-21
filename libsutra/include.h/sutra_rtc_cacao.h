// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_rtc_cacao.h
//! \ingroup   libsutra
//! \class     SutraRtcCacao
//! \brief     this class provides the rtc_cacao features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24


#ifndef SUTRA_RTC_CACAO_H_
#define SUTRA_RTC_CACAO_H_

#include <Cacao.h>

#include <sutra_rtc.h>
#include <sutra_target.h>
#include <sutra_wfs.h>

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

  long framecounter_;

  int nslp_;
  int ncmd_;
  int nvalid_;

  bool is_initialised_;

 public:
  SutraRtcCacao(std::string interface_cal_frame_name, std::string interface_loop_frame_name);
  ~SutraRtcCacao();

  void publish();

 private:
  void allocate_buffers();
};

#endif /* SUTRA_RTC_CACAO_H_ */
