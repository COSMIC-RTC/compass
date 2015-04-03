#ifndef SUTRA_RTC_LISTENER_IMPL
#define SUTRA_RTC_LISTENER_IMPL

#include "libBRAMATypeSupportImpl.h"
#include <dds/DdsDcpsSubscriptionC.h>
#include <ace/Synch.h>

#if !defined (ACE_LACKS_PRAGMA_ONCE)
#pragma once
#endif /* ACE_LACKS_PRAGMA_ONCE */

class sutra_rtc_brama;
class libBRAMACommon_Export sutra_rtc_bramaListenerImpl: public virtual OpenDDS::DCPS::LocalObject<
    DDS::DataReaderListener> {
  private:
    sutra_rtc_brama *rtc;
    ACE_Mutex lock_;
  public:
    //Constructor
    sutra_rtc_bramaListenerImpl();

    //Destructor
    virtual ~sutra_rtc_bramaListenerImpl(void);

    // app-specific
    void attach_rtc(sutra_rtc_brama *rtc);

    // must also override:
    virtual void on_data_available(DDS::DataReader_ptr reader)
        throw (CORBA::SystemException);
    virtual void on_requested_deadline_missed(
        DDS::DataReader_ptr reader,
        const DDS::RequestedDeadlineMissedStatus & status)
            throw (CORBA::SystemException);
    virtual void on_requested_incompatible_qos(
        DDS::DataReader_ptr reader,
        const DDS::RequestedIncompatibleQosStatus & status)
            throw (CORBA::SystemException);
    virtual void on_liveliness_changed(
        DDS::DataReader_ptr reader, const DDS::LivelinessChangedStatus & status)
            throw (CORBA::SystemException);
    virtual void on_subscription_matched(
        DDS::DataReader_ptr reader,
        const DDS::SubscriptionMatchedStatus & status)
            throw (CORBA::SystemException);
    virtual void on_sample_rejected(DDS::DataReader_ptr reader,
                                    const DDS::SampleRejectedStatus& status)
                                        throw (CORBA::SystemException);
    virtual void on_sample_lost(DDS::DataReader_ptr reader,
                                const DDS::SampleLostStatus& status)
                                    throw (CORBA::SystemException);

};
#endif /* SUTRA_RTC_LISTENER_IMPL  */
