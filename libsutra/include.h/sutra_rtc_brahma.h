// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_rtc_brahma.h
//! \ingroup   libsutra
//! \class     SutraRtcBrahma
//! \brief     this class provides the rtc_brahma features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef SUTRA_RTC_BRAHMA_H_
#define SUTRA_RTC_BRAHMA_H_

#include <BRAHMA_context.h>
#include <sutra_rtc.h>
#include <sutra_rtc_brahmaListenerImpl.h>
#include <sutra_target.h>
#include <sutra_wfs.h>

template <typename T>
class SutraRtcBrahma : public SutraRtc<float, T, float> {
 private:
  DDS::Subscriber_var sub;
  DDS::Publisher_var pub;

  DDS::DataReaderListener_var cmd_listener;
  sutra_rtc_brahmaListenerImpl<T> *cmd_listener_servant;
  DDS::DataReader_var cmd_dr;

  DDS::DataWriter_var superframe_base_dw;
  BRAHMA::SuperFrameDataWriter_var superframe_dw;
  DDS::InstanceHandle_t superframe_handle;

  DDS::DataWriter_var megaframe_base_dw;
  BRAHMA::MegaFrameDataWriter_var megaframe_dw;
  DDS::InstanceHandle_t megaframe_handle;

  CORBA::Octet *buff_wfs;
  CORBA::Octet *buff_wfs_phase;
  CORBA::Octet *buff_intensities;
  CORBA::Octet *buff_slopes;
  CORBA::Octet *buff_commands;
  CORBA::Octet *buff_target;
  CORBA::Octet *buff_target_phase;

  CORBA::ULong *dims_wfs;
  CORBA::ULong *dims_wfs_phase;
  CORBA::ULong *dims_intensities;
  CORBA::ULong *dims_slopes;
  CORBA::ULong *dims_commands;
  CORBA::ULong *dims_target;
  CORBA::ULong *dims_target_phase;

  int64_t framecounter;
  ACE_Mutex lock_;

  int32_t wfs_size;
  int32_t wfs_phase_size;
  SutraSensors *wfs;
  int32_t target_size;
  int32_t target_phase_size;
  SutraTarget *target;

  int32_t nslp;
  int32_t ncmd;
  int32_t nvalid;

  int32_t is_initialised;

 public:
  SutraRtcBrahma(CarmaContext *context, SutraSensors *wfs,
                   SutraTarget *target, ACE_TCHAR *name);
  ~SutraRtcBrahma();

  void publish();

 private:
  void allocate_buffers();
};

#endif /* SUTRA_RTC_BRAHMA_H_ */
