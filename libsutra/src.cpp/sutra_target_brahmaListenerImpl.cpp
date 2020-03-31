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

//! \file      sutra_target_brahmaListenerImpl.cpp
//! \ingroup   libsutra
//! \class     sutra_target_brahmaListenerImpl
//! \brief     this class provides the target_brahmaListenerImpl features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifdef USE_BRAHMA

#include <sutra_target_brahmaListenerImpl.h>
#include "ace/streams.h"
#include "dds/DCPS/Service_Participant.h"
#include "libBRAHMATypeSupportC.h"
#include "libBRAHMATypeSupportImpl.h"

#include "sutra_target_brahma.h"

// Constructor
sutra_target_brahmaListenerImpl::sutra_target_brahmaListenerImpl()
    : target(0L) {}

// Destructor
sutra_target_brahmaListenerImpl::~sutra_target_brahmaListenerImpl(void) {}

// app-specific
void sutra_target_brahmaListenerImpl::attach_target(
    sutra_target_brahma* target_) {
  target = target_;
}

void sutra_target_brahmaListenerImpl::on_data_available(
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
void sutra_target_brahmaListenerImpl::on_requested_deadline_missed(
    DDS::DataReader_ptr reader,
    const DDS::RequestedDeadlineMissedStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_requested_deadline_missed" <<
  //  endl;
}
void sutra_target_brahmaListenerImpl::on_requested_incompatible_qos(
    DDS::DataReader_ptr reader,
    const DDS::RequestedIncompatibleQosStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_requested_incompatible_qos" <<
  //  endl;
}
void sutra_target_brahmaListenerImpl::on_liveliness_changed(
    DDS::DataReader_ptr reader,
    const DDS::LivelinessChangedStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_liveliness_changed" << endl;
}
void sutra_target_brahmaListenerImpl::on_subscription_matched(
    DDS::DataReader_ptr reader,
    const DDS::SubscriptionMatchedStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_subscription_matched" << endl;
}
void sutra_target_brahmaListenerImpl::on_sample_rejected(
    DDS::DataReader_ptr reader,
    const DDS::SampleRejectedStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_sample_rejected" << endl;
}
void sutra_target_brahmaListenerImpl::on_sample_lost(
    DDS::DataReader_ptr reader,
    const DDS::SampleLostStatus& status) noexcept(false) {
  //  cerr << "CommandDataReaderListenerImpl::on_sample_lost" << endl;
}

#endif /* USE_BRAHMA */
