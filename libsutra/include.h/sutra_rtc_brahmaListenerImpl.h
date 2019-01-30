#ifndef SUTRA_RTC_LISTENER_IMPL
#define SUTRA_RTC_LISTENER_IMPL

#include <ace/Synch.h>
#include <dds/DdsDcpsSubscriptionC.h>
#include "libBRAHMATypeSupportImpl.h"

#if !defined(ACE_LACKS_PRAGMA_ONCE)
#pragma once
#endif /* ACE_LACKS_PRAGMA_ONCE */

template <typename T>
class sutra_rtc_brahma;

template <typename T>
class libBRAHMACommon_Export sutra_rtc_brahmaListenerImpl
    : public virtual OpenDDS::DCPS::LocalObject<DDS::DataReaderListener> {
 private:
  sutra_rtc_brahma<T>* rtc;
  ACE_Mutex lock_;

 public:
  // Constructor
  sutra_rtc_brahmaListenerImpl();

  // Destructor
  virtual ~sutra_rtc_brahmaListenerImpl(void);

  // app-specific
  void attach_rtc(sutra_rtc_brahma<T>* rtc);

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
