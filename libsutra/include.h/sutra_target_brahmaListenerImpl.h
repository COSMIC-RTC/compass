#ifndef SUTRA_TARGET_LISTENER_IMPL
#define SUTRA_TARGET_LISTENER_IMPL

#include <ace/Synch.h>
#include <dds/DdsDcpsSubscriptionC.h>
#include "libBRAHMATypeSupportImpl.h"

#if !defined(ACE_LACKS_PRAGMA_ONCE)
#pragma once
#endif /* ACE_LACKS_PRAGMA_ONCE */

class sutra_target_brahma;
class libBRAHMACommon_Export sutra_target_brahmaListenerImpl
    : public virtual OpenDDS::DCPS::LocalObject<DDS::DataReaderListener> {
 private:
  sutra_target_brahma* target;
  ACE_Mutex lock_;

 public:
  // Constructor
  sutra_target_brahmaListenerImpl();

  // Destructor
  virtual ~sutra_target_brahmaListenerImpl(void);

  // app-specific
  void attach_target(sutra_target_brahma* target);

  // must also override:
  virtual void on_data_available(DDS::DataReader_ptr reader) throw(
      CORBA::SystemException);
  virtual void on_requested_deadline_missed(
      DDS::DataReader_ptr reader, const DDS::RequestedDeadlineMissedStatus&
                                      status) throw(CORBA::SystemException);
  virtual void on_requested_incompatible_qos(
      DDS::DataReader_ptr reader, const DDS::RequestedIncompatibleQosStatus&
                                      status) throw(CORBA::SystemException);
  virtual void on_liveliness_changed(
      DDS::DataReader_ptr reader,
      const DDS::LivelinessChangedStatus& status) throw(CORBA::SystemException);
  virtual void on_subscription_matched(
      DDS::DataReader_ptr reader, const DDS::SubscriptionMatchedStatus&
                                      status) throw(CORBA::SystemException);
  virtual void on_sample_rejected(
      DDS::DataReader_ptr reader,
      const DDS::SampleRejectedStatus& status) throw(CORBA::SystemException);
  virtual void on_sample_lost(
      DDS::DataReader_ptr reader,
      const DDS::SampleLostStatus& status) throw(CORBA::SystemException);
};
#endif /* SUTRA_TARGET_LISTENER_IMPL  */
