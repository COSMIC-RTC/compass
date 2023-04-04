// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_target_brahmaListenerImpl.cpp
//! \ingroup   libsutra
//! \class     SutraTargetBrahmaListenerImpl
//! \brief     this class provides the target_brahmaListenerImpl features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24

#ifdef USE_BRAHMA

#include <sutra_target_brahmaListenerImpl.h>
#include "ace/streams.h"
#include "dds/DCPS/Service_Participant.h"
#include "libBRAHMATypeSupportC.h"
#include "libBRAHMATypeSupportImpl.h"

#include "sutra_target_brahma.h"

// Constructor
SutraTargetBrahmaListenerImpl::SutraTargetBrahmaListenerImpl()
    : target(0L) {}

// Destructor
SutraTargetBrahmaListenerImpl::~SutraTargetBrahmaListenerImpl(void) {}

// app-specific
void SutraTargetBrahmaListenerImpl::attach_target(
    SutraTargetBrahma* target_) {
  target = target_;
}

void SutraTargetBrahmaListenerImpl::on_data_available(
    DDS::DataReader_ptr reader) noexcept(false) {
  try {
    BRAHMA::CommandDataReader_var target_cmd_dr =
        BRAHMA::CommandDataReader::_narrow(reader);

    if (CORBA::is_nil(target_cmd_dr.in())) {
      cerr << "CommandDataReaderListenerImpl:: "
           << "on_data_available:"
           << " _narrow failed." << endl;
      ACE_OS::exit(1);
    }
    int count = 0;
    while (true) {
      DDS::SampleInfo si;
      BRAHMA::Command cmd;

      //      ACE_Time_Value ace_wait(0, 250000);
      //      while (!this->lastCmdApplied)
      //        ACE_OS::sleep(ace_wait);
      {
        // ACE_Guard<ACE_Mutex> guard(this->lock_);
        DDS::ReturnCode_t status = target_cmd_dr->take_next_sample(cmd, si);
        if (status == DDS::RETCODE_OK) {
          ++count;
          {
            ACE_Guard<ACE_Mutex> guard(this->lock_);

            try {
              const string type(cmd.name.in());
              //              DEBUG_TRACE("get command type %s", type.c_str());

              std::vector<string> cmd_splited;
              carma_utils::split(cmd_splited, type, '_');

              if (cmd_splited[0] == "canapassExposure") {
                int ntarget = carma_utils::from_string<int>(cmd_splited[1]);
                int subsample = carma_utils::from_string<int>(cmd_splited[2]);
                DEBUG_TRACE("Updating canapassExposure on target %d to %d",
                            ntarget, subsample);

                target->set_subsample(ntarget, subsample);

                //                if ((cmd.dimensions.length() != 1) &&
                //                (cmd.dimensions[0] == ncmd)) {
                //                  BRAHMA_DEBUG_TRACE("wrong dimensions : %d
                //                  %d",
                //                                    cmd.dimensions.length(),
                //                                    cmd.dimensions[0]);
                //                  BRAHMA_DEBUG_TRACE("it should be : 1 %d",
                //                  ncmd); throw CORBA::BAD_PARAM();
                //                }
              } else {
                throw "Unknown parameter";
              }
            } catch (char const* msg) {
              // BRAHMA_DEBUG_TRACE("%s",msg);
            } catch (CORBA::BAD_PARAM& p) {
              std::ostringstream stm;
              stm << p;
              BRAHMA_DEBUG_TRACE("%s", stm.str().c_str());
            }
          }
        } else if (status == DDS::RETCODE_NO_DATA) {
          //          cerr << "INFO: Event reading complete after " << count <<
          //          " samples." << endl; cerr << "ERROR: reader received " <<
          //          "DDS::RETCODE_NO_DATA!" << endl;
          break;
        } else {
          cerr << "ERROR: read BRAHMACommand: Error: " << status << endl;
          break;
        }
      }
    }
  } catch (CORBA::Exception& e) {
    cerr << "Exception caught in read:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

// must also override:
void SutraTargetBrahmaListenerImpl::on_requested_deadline_missed(
    DDS::DataReader_ptr reader,
    const DDS::RequestedDeadlineMissedStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_requested_deadline_missed" <<
  //  endl;
}
void SutraTargetBrahmaListenerImpl::on_requested_incompatible_qos(
    DDS::DataReader_ptr reader,
    const DDS::RequestedIncompatibleQosStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_requested_incompatible_qos" <<
  //  endl;
}
void SutraTargetBrahmaListenerImpl::on_liveliness_changed(
    DDS::DataReader_ptr reader,
    const DDS::LivelinessChangedStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_liveliness_changed" << endl;
}
void SutraTargetBrahmaListenerImpl::on_subscription_matched(
    DDS::DataReader_ptr reader,
    const DDS::SubscriptionMatchedStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_subscription_matched" << endl;
}
void SutraTargetBrahmaListenerImpl::on_sample_rejected(
    DDS::DataReader_ptr reader,
    const DDS::SampleRejectedStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_sample_rejected" << endl;
}
void SutraTargetBrahmaListenerImpl::on_sample_lost(
    DDS::DataReader_ptr reader,
    const DDS::SampleLostStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_sample_lost" << endl;
}

#endif /* USE_BRAHMA */
