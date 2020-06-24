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

//! \file      sutra_rtc_brahmaListenerImpl.cpp
//! \ingroup   libsutra
//! \class     sutra_rtc_brahmaListenerImpl
//! \brief     this class provides the rtc_brahmaListenerImpl features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifdef USE_BRAHMA

#include <sutra_rtc_brahmaListenerImpl.h>
#include "ace/streams.h"
#include "dds/DCPS/Service_Participant.h"
#include "libBRAHMATypeSupportImpl.h"

#include "sutra_rtc_brahma.h"

// Constructor
template <typename T>
sutra_rtc_brahmaListenerImpl<T>::sutra_rtc_brahmaListenerImpl() : rtc(0L) {}

// Destructor
template <typename T>
sutra_rtc_brahmaListenerImpl<T>::~sutra_rtc_brahmaListenerImpl(void) {}

// app-specific
template <typename T>
void sutra_rtc_brahmaListenerImpl<T>::attach_rtc(sutra_rtc_brahma<T> *rtc_) {
  rtc = rtc_;
}

template <typename T>
void sutra_rtc_brahmaListenerImpl<T>::on_data_available(
    DDS::DataReader_ptr reader) noexcept(false) {
  //  DEBUG_TRACE("Entering in
  //  sutra_rtc_brahmaListenerImpl::on_data_available");

  try {
    BRAHMA::CommandDataReader_var rtc_cmd_dr =
        BRAHMA::CommandDataReader::_narrow(reader);

    if (CORBA::is_nil(rtc_cmd_dr.in())) {
      cerr << "CommandDataReaderListenerImpl:: "
           << "on_data_available:"
           << " _narrow failed." << endl;
      ACE_OS::exit(1);
    }
    DDS::SampleInfo si;
    BRAHMA::Command cmd;
    while (DDS::RETCODE_OK == rtc_cmd_dr->take_next_sample(cmd, si)) {
      if (si.valid_data) {
        //        ACE_Guard < ACE_Mutex > guard(this->lock_);
        try {
          const string type(cmd.name.in());
          //          DEBUG_TRACE("get command type %s", type.c_str());

          std::vector<string> cmd_splited;
          carma_utils::split(cmd_splited, type, '_');

          if (cmd_splited[0] == "gain") {
            int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);
            DEBUG_TRACE("Updating mgain on controller %d", ncontrol);

            unsigned int ncmd = rtc->d_control[ncontrol]->nactu();
            T *data = (T *)cmd.data.get_buffer();

            if ((cmd.dimensions.length() != 1) && (cmd.dimensions[0] == ncmd)) {
              BRAHMA_DEBUG_TRACE("wrong dimensions : %d %d",
                                 cmd.dimensions.length(), cmd.dimensions[0]);
              BRAHMA_DEBUG_TRACE("it should be : 1 %d", ncmd);
              throw CORBA::BAD_PARAM();
            }
            if (rtc->d_control[ncontrol]->get_type() == "ls") {
              sutra_controller_ls *control =
                  dynamic_cast<sutra_controller_ls *>(rtc->d_control[ncontrol]);
              control->set_mgain(data);
              return;
            } else if (rtc->d_control[ncontrol]->get_type() == "mv") {
              sutra_controller_mv *control =
                  dynamic_cast<sutra_controller_mv *>(rtc->d_control[ncontrol]);
              control->set_mgain(data);
              return;
            } else if (rtc->d_control[ncontrol]->get_type() == "generic") {
              sutra_controller_generic *control =
                  dynamic_cast<sutra_controller_generic *>(
                      rtc->d_control[ncontrol]);
              control->set_mgain(data);
              return;
            } else {
              BRAHMA_DEBUG_TRACE("controller %d must be a ls or mv controller",
                                 ncontrol);
              throw CORBA::BAD_PARAM();
            }
          } else if (cmd_splited[0] == "globalGain") {
            int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);
            float gain = carma_utils::from_string<float>(cmd_splited[2]);
            DEBUG_TRACE("Updating gain num %d to %f", ncontrol, gain);
            if (rtc->d_control[ncontrol]->get_type() == "ls") {
              sutra_controller_ls *control =
                  dynamic_cast<sutra_controller_ls *>(rtc->d_control[ncontrol]);
              control->set_gain(gain);
              return;
            } else if (rtc->d_control[ncontrol]->get_type() == "mv") {
              sutra_controller_mv *control =
                  dynamic_cast<sutra_controller_mv *>(rtc->d_control[ncontrol]);
              control->set_gain(gain);
              return;
            } else {
              BRAHMA_DEBUG_TRACE("controller %d must be a ls or mv controller",
                                 ncontrol);
              throw CORBA::BAD_PARAM();
            }
          } else if (cmd_splited[0] == "openLoop") {
            int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);
            int openloop = carma_utils::from_string<int>(cmd_splited[2]);
            DEBUG_TRACE("Updating openloop num %d to %d", ncontrol, openloop);
            rtc->d_control[ncontrol]->set_openloop(openloop);
          } else if (cmd_splited[0] == "CM") {
            int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);

            DEBUG_TRACE("Updating control matrix num %d", ncontrol);

            unsigned int nslope = rtc->d_control[ncontrol]->nslope();
            unsigned int ncmd = rtc->d_control[ncontrol]->nactu();

            if (cmd.dimensions.length() != 2 || cmd.dimensions[0] != ncmd ||
                cmd.dimensions[1] != nslope) {
              std::stringstream ss;
              ss << "wrong dimensions :";
              for (unsigned int i = 0; i <= cmd.dimensions.length(); i++)
                ss << " " << cmd.dimensions[i];
              BRAHMA_DEBUG_TRACE("%s", ss.str().c_str());
              BRAHMA_DEBUG_TRACE("it should be : %d %d", ncmd, nslope);
              throw CORBA::BAD_PARAM();
            }

            T *data = (T *)cmd.data.get_buffer();
            if (rtc->d_control[ncontrol]->get_type() == "ls") {
              sutra_controller_ls *control =
                  dynamic_cast<sutra_controller_ls *>(rtc->d_control[ncontrol]);
              control->set_cmat(data);
              return;
            } else if (rtc->d_control[ncontrol]->get_type() == "mv") {
              sutra_controller_mv *control =
                  dynamic_cast<sutra_controller_mv *>(rtc->d_control[ncontrol]);
              control->set_cmat(data);
              return;
            } else if (rtc->d_control[ncontrol]->get_type() == "generic") {
              sutra_controller_generic *control =
                  dynamic_cast<sutra_controller_generic *>(
                      rtc->d_control[ncontrol]);
              control->set_cmat(data);
              return;
            } else {
              BRAHMA_DEBUG_TRACE(
                  "controller %d must be a ls, mv or generic controller",
                  ncontrol);
              throw CORBA::BAD_PARAM();
            }
          } else if (cmd_splited[0] == "E") {
            int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);

            DEBUG_TRACE("Updating E matrix num %d", ncontrol);

            unsigned int ncmd = rtc->d_control[ncontrol]->nactu();

            if (cmd.dimensions.length() != 2 || cmd.dimensions[0] != ncmd ||
                cmd.dimensions[1] != ncmd) {
              BRAHMA_DEBUG_TRACE("wrong dimensions : %d %d %d",
                                 cmd.dimensions.length(), cmd.dimensions[0],
                                 cmd.dimensions[1]);
              BRAHMA_DEBUG_TRACE("it should be : 2 %d %d", ncmd, ncmd);
              throw CORBA::BAD_PARAM();
            }

            T *data = (T *)cmd.data.get_buffer();
            if (rtc->d_control[ncontrol]->get_type() == "generic") {
              sutra_controller_generic *control =
                  dynamic_cast<sutra_controller_generic *>(
                      rtc->d_control[ncontrol]);
              control->set_matE(data);
              return;
            } else {
              BRAHMA_DEBUG_TRACE("controller %d must be a generic controller",
                                 ncontrol);
              throw CORBA::BAD_PARAM();
            }
          } else if (cmd_splited[0] == "decayFactor") {
            int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);

            DEBUG_TRACE("Updating decayFactor num %d", ncontrol);

            unsigned int ncmd = rtc->d_control[ncontrol]->nactu();

            if (cmd.dimensions.length() != 1 || cmd.dimensions[0] != ncmd) {
              BRAHMA_DEBUG_TRACE("wrong dimensions : %d %d",
                                 cmd.dimensions.length(), cmd.dimensions[0]);
              BRAHMA_DEBUG_TRACE("it should be : 1 %d", ncmd);
              throw CORBA::BAD_PARAM();
            }

            T *data = (T *)cmd.data.get_buffer();
            if (rtc->d_control[ncontrol]->get_type() == "generic") {
              sutra_controller_generic *control =
                  dynamic_cast<sutra_controller_generic *>(
                      rtc->d_control[ncontrol]);
              control->set_decayFactor(data);
              return;
            } else {
              BRAHMA_DEBUG_TRACE("controller %d must be a generic controller",
                                 ncontrol);
              throw CORBA::BAD_PARAM();
            }
          } else if (cmd_splited[0] == "commandLaw") {
            int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);
            string law = cmd_splited[2];

            DEBUG_TRACE("Updating commandLaw num %d to %s", ncontrol,
                        law.c_str());

            if (rtc->d_control[ncontrol]->get_type() == "generic") {
              sutra_controller_generic *control =
                  dynamic_cast<sutra_controller_generic *>(
                      rtc->d_control[ncontrol]);
              control->set_commandlaw(law);
              return;
            } else {
              BRAHMA_DEBUG_TRACE("controller %d must be a generic controller",
                                 ncontrol);
              throw CORBA::BAD_PARAM();
            }
          } else if (cmd_splited[0] == "PertVolt") {
            int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);
            DEBUG_TRACE("Updating perturbation voltages on controller %d",
                        ncontrol);

            unsigned int ncmd = rtc->d_control[ncontrol]->nactu();
            T *data = (T *)cmd.data.get_buffer();

            if (cmd.dimensions.length() == 1) {
              if (cmd.dimensions[0] == ncmd) {
                rtc->d_control[ncontrol]->set_perturbcom(data, 1);
                return;
              } else {
                BRAHMA_DEBUG_TRACE("wrong dimensions : %d %d",
                                   cmd.dimensions.length(), cmd.dimensions[0]);
                BRAHMA_DEBUG_TRACE("it should be : 1 %d", ncmd);
                throw CORBA::BAD_PARAM();
              }
            } else if (cmd.dimensions.length() == 2) {
              if (cmd.dimensions[0] == ncmd) {
                rtc->d_control[ncontrol]->set_perturbcom(data,
                                                         cmd.dimensions[1]);
                return;
              } else {
                BRAHMA_DEBUG_TRACE("wrong dimensions : %d %d %d",
                                   cmd.dimensions.length(), cmd.dimensions[0],
                                   cmd.dimensions[1]);
                BRAHMA_DEBUG_TRACE("it should be : 2 %d nb_elem", ncmd);
                throw CORBA::BAD_PARAM();
              }
            } else {
              CORBA::ULong *data = (CORBA::ULong *)cmd.dimensions.get_buffer();

              cerr << "wrong dimensions: ";
              // Print in Normal order
              std::copy(data, data + data[0] + 1,
                        std::ostream_iterator<CORBA::ULong>(std::cout, ","));
              std::cout << "\n";
              throw CORBA::BAD_PARAM();
            }
          } else {
            std::ostringstream stm;
            stm << "Nothing to do with " << cmd_splited[0];
            throw stm.str().c_str();
          }
        } catch (char const *msg) {
          BRAHMA_DEBUG_TRACE("%s", msg);
        } catch (CORBA::BAD_PARAM &p) {
          std::ostringstream stm;
          stm << p;
          BRAHMA_DEBUG_TRACE("%s", stm.str().c_str());
        }
      } else {
        //         BRAHMA_DEBUG_TRACE("sutra_rtc_brahmaListenerImpl::on_data_available
        //         received a non-data sample. ");
      }
    }
  } catch (CORBA::Exception &e) {
    cerr << "Exception caught in read:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

// must also override:
template <typename T>
void sutra_rtc_brahmaListenerImpl<T>::on_requested_deadline_missed(
    DDS::DataReader_ptr reader,
    const DDS::RequestedDeadlineMissedStatus &status) noexcept(false) {
  //   BRAHMA_DEBUG_TRACE(
  //       "CommandDataReaderListenerImpl::on_requested_deadline_missed");
  //  cerr << "CommandDataReaderListenerImpl::on_requested_deadline_missed" <<
  //  endl;
}
template <typename T>
void sutra_rtc_brahmaListenerImpl<T>::on_requested_incompatible_qos(
    DDS::DataReader_ptr reader,
    const DDS::RequestedIncompatibleQosStatus &status) noexcept(false) {
  //   BRAHMA_DEBUG_TRACE(
  //       "CommandDataReaderListenerImpl::on_requested_incompatible_qos");
  //  cerr << "CommandDataReaderListenerImpl::on_requested_incompatible_qos" <<
  //  endl;
}
template <typename T>
void sutra_rtc_brahmaListenerImpl<T>::on_liveliness_changed(
    DDS::DataReader_ptr reader,
    const DDS::LivelinessChangedStatus &status) noexcept(false) {
  //   BRAHMA_DEBUG_TRACE(
  //       "CommandDataReaderListenerImpl::on_liveliness_changed");
  //  cerr << "CommandDataReaderListenerImpl::on_liveliness_changed" << endl;
}
template <typename T>
void sutra_rtc_brahmaListenerImpl<T>::on_subscription_matched(
    DDS::DataReader_ptr reader,
    const DDS::SubscriptionMatchedStatus &status) noexcept(false) {
  //   BRAHMA_DEBUG_TRACE(
  //       "CommandDataReaderListenerImpl::on_subscription_matched");
  //  cerr << "CommandDataReaderListenerImpl::on_subscription_matched" << endl;
}
template <typename T>
void sutra_rtc_brahmaListenerImpl<T>::on_sample_rejected(
    DDS::DataReader_ptr reader,
    const DDS::SampleRejectedStatus &status) noexcept(false) {
  //   BRAHMA_DEBUG_TRACE(
  //       "CommandDataReaderListenerImpl::on_sample_rejected");
  //  cerr << "CommandDataReaderListenerImpl::on_sample_rejected" << endl;
}
template <typename T>
void sutra_rtc_brahmaListenerImpl<T>::on_sample_lost(
    DDS::DataReader_ptr reader,
    const DDS::SampleLostStatus &status) noexcept(false) {
  //   BRAHMA_DEBUG_TRACE(
  //       "CommandDataReaderListenerImpl::on_sample_lost");
  //  cerr << "CommandDataReaderListenerImpl::on_sample_lost" << endl;
}

template class sutra_rtc_brahmaListenerImpl<float>;
#ifdef CAN_DO_HALF
template class sutra_rtc_brahmaListenerImpl<half>;
#endif

#endif /* USE_BRAHMA */
