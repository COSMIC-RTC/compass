// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      SutraRtcBrahmaListenerImpl.h
//! \ingroup   libsutra
//! \class     SutraRtcBrahmaListenerImpl
//! \brief     this class provides the rtc_brahmaListenerImpl features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#ifndef SUTRA_RTC_LISTENER_IMPL
#define SUTRA_RTC_LISTENER_IMPL

#include <ace/Synch.h>
#include <dds/DdsDcpsSubscriptionC.h>
#include "libBRAHMATypeSupportImpl.h"

#if !defined(ACE_LACKS_PRAGMA_ONCE)
#pragma once
#endif /* ACE_LACKS_PRAGMA_ONCE */

template <typename T>
class SutraRtcBrahma;

template <typename T>
class libBRAHMACommon_Export SutraRtcBrahmaListenerImpl
    : public virtual OpenDDS::DCPS::LocalObject<DDS::DataReaderListener> {
 private:
  SutraRtcBrahma<T>* rtc;
  ACE_Mutex lock_;

 public:
  // Constructor
  SutraRtcBrahmaListenerImpl();

  // Destructor
  virtual ~SutraRtcBrahmaListenerImpl(void);

  // app-specific
  void attach_rtc(SutraRtcBrahma<T>* rtc);

  // must also override:
  virtual void on_data_available(DDS::DataReader_ptr reader) noexcept(false);
  virtual void on_requested_deadline_missed(
      DDS::DataReader_ptr reader,
      const DDS::RequestedDeadlineMissedStatus& status) noexcept(false);
  virtual void on_requested_incompatible_qos(
      DDS::DataReader_ptr reader,
      const DDS::RequestedIncompatibleQosStatus& status) noexcept(false);
  virtual void on_liveliness_changed(
      DDS::DataReader_ptr reader,
      const DDS::LivelinessChangedStatus& status) noexcept(false);
  virtual void on_subscription_matched(
      DDS::DataReader_ptr reader,
      const DDS::SubscriptionMatchedStatus& status) noexcept(false);
  virtual void on_sample_rejected(
      DDS::DataReader_ptr reader,
      const DDS::SampleRejectedStatus& status) noexcept(false);
  virtual void on_sample_lost(
      DDS::DataReader_ptr reader,
      const DDS::SampleLostStatus& status) noexcept(false);
};
#endif /* SUTRA_RTC_LISTENER_IMPL  */
