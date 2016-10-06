#ifdef USE_BRAMA

#include "libBRAMATypeSupportImpl.h"
#include<sutra_rtc_bramaListenerImpl.h>
#include "dds/DCPS/Service_Participant.h"
#include "ace/streams.h"

#include "sutra_rtc_brama.h"

//Constructor
sutra_rtc_bramaListenerImpl::sutra_rtc_bramaListenerImpl(): rtc(0L) {

}

//Destructor
sutra_rtc_bramaListenerImpl::~sutra_rtc_bramaListenerImpl(void) {
}

// app-specific
void sutra_rtc_bramaListenerImpl::attach_rtc(sutra_rtc_brama *rtc_){
  rtc = rtc_;
}

void sutra_rtc_bramaListenerImpl::on_data_available(DDS::DataReader_ptr reader)
    throw (CORBA::SystemException) {
  try {
    BRAMA::CommandDataReader_var rtc_cmd_dr = BRAMA::CommandDataReader::_narrow(
        reader);

    if (CORBA::is_nil(rtc_cmd_dr.in())) {
      cerr << "CommandDataReaderListenerImpl:: " << "on_data_available:"
           << " _narrow failed." << endl;
      ACE_OS::exit(1);
    }
    int count = 0;
    while (true) {
      DDS::SampleInfo si;
      BRAMA::Command cmd;

//      ACE_Time_Value ace_wait(0, 250000);
//      while (!this->lastCmdApplied)
//        ACE_OS::sleep(ace_wait);
      {
        //ACE_Guard<ACE_Mutex> guard(this->lock_);
        DDS::ReturnCode_t status = rtc_cmd_dr->take_next_sample(cmd, si);
        if (status == DDS::RETCODE_OK) {
          ++count;
          {
            ACE_Guard<ACE_Mutex> guard(this->lock_);

            try {
              const string type(cmd.name.in());
              DEBUG_TRACE("get command type %s", type.c_str());

              std::vector<string> cmd_splited;
              carma_utils::split(cmd_splited, type, '_');

              if (cmd_splited[0] == "gain") {
                int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);
                DEBUG_TRACE("Updating mgain on controller %d", ncontrol);

                unsigned int ncmd = rtc->d_control[ncontrol]->nactu();
                CORBA::Float *data = (CORBA::Float*) cmd.data.get_buffer();

                if ((cmd.dimensions[0] != 1) && (cmd.dimensions[1] == ncmd)) {
                  BRAMA_DEBUG_TRACE("wrong dimensions : %d %d",
                                    cmd.dimensions[0], cmd.dimensions[1]);
                  BRAMA_DEBUG_TRACE("it should be : 1 %d", ncmd);
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
                      dynamic_cast<sutra_controller_generic *>(rtc->d_control[ncontrol]);
                  control->set_mgain(data);
                  return;
                } else {
                  BRAMA_DEBUG_TRACE(
                      "controller %d must be a ls or mv controller", ncontrol);
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
                  BRAMA_DEBUG_TRACE(
                      "controller %d must be a ls or mv controller", ncontrol);
                  throw CORBA::BAD_PARAM();
                }
              } else if (cmd_splited[0] == "openLoop") {
                int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);
                int openloop = carma_utils::from_string<int>(cmd_splited[2]);
                DEBUG_TRACE("Updating openloop num %d to %d", ncontrol,
                            openloop);
                rtc->d_control[ncontrol]->set_openloop(openloop);
              } else if (cmd_splited[0] == "CM") {
                int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);

                DEBUG_TRACE("Updating control matrix num %d", ncontrol);

                unsigned int nslope = rtc->d_control[ncontrol]->nslope();
                unsigned int ncmd = rtc->d_control[ncontrol]->nactu();

                if (cmd.dimensions[0] != 2 || cmd.dimensions[1] != ncmd
                    || cmd.dimensions[2] != nslope) {
                  BRAMA_DEBUG_TRACE("wrong dimensions : %d %d %d",
                                    cmd.dimensions[0], cmd.dimensions[1],
                                    cmd.dimensions[2]);
                  BRAMA_DEBUG_TRACE("it should be : 2 %d %d", ncmd, nslope);
                  throw CORBA::BAD_PARAM();
                }

                CORBA::Float *data = (CORBA::Float*) cmd.data.get_buffer();
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
                      dynamic_cast<sutra_controller_generic *>(rtc->d_control[ncontrol]);
                  control->set_cmat(data);
                  return;
                } else {
                  BRAMA_DEBUG_TRACE(
                      "controller %d must be a ls, mv or generic controller",
                      ncontrol);
                  throw CORBA::BAD_PARAM();
                }
              } else if (cmd_splited[0] == "E") {
                int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);

                DEBUG_TRACE("Updating E matrix num %d", ncontrol);

                unsigned int ncmd = rtc->d_control[ncontrol]->nactu();

                if (cmd.dimensions[0] != 2 || cmd.dimensions[1] != ncmd
                    || cmd.dimensions[2] != ncmd) {
                  BRAMA_DEBUG_TRACE("wrong dimensions : %d %d %d",
                                    cmd.dimensions[0], cmd.dimensions[1],
                                    cmd.dimensions[2]);
                  BRAMA_DEBUG_TRACE("it should be : 2 %d %d", ncmd, ncmd);
                  throw CORBA::BAD_PARAM();
                }

                CORBA::Float *data = (CORBA::Float*) cmd.data.get_buffer();
                if (rtc->d_control[ncontrol]->get_type() == "generic") {
                  sutra_controller_generic *control =
                      dynamic_cast<sutra_controller_generic *>(rtc->d_control[ncontrol]);
                  control->set_matE(data);
                  return;
                } else {
                  BRAMA_DEBUG_TRACE(
                      "controller %d must be a generic controller", ncontrol);
                  throw CORBA::BAD_PARAM();
                }
              } else if (cmd_splited[0] == "decayFactor") {
                int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);

                DEBUG_TRACE("Updating decayFactor num %d", ncontrol);

                unsigned int ncmd = rtc->d_control[ncontrol]->nactu();

                if (cmd.dimensions[0] != 1 || cmd.dimensions[1] != ncmd) {
                  BRAMA_DEBUG_TRACE("wrong dimensions : %d %d",
                                    cmd.dimensions[0], cmd.dimensions[1]);
                  BRAMA_DEBUG_TRACE("it should be : 1 %d", ncmd);
                  throw CORBA::BAD_PARAM();
                }

                CORBA::Float *data = (CORBA::Float*) cmd.data.get_buffer();
                if (rtc->d_control[ncontrol]->get_type() == "generic") {
                  sutra_controller_generic *control =
                      dynamic_cast<sutra_controller_generic *>(rtc->d_control[ncontrol]);
                  control->set_decayFactor(data);
                  return;
                } else {
                  BRAMA_DEBUG_TRACE(
                      "controller %d must be a generic controller", ncontrol);
                  throw CORBA::BAD_PARAM();
                }
              } else if (cmd_splited[0] == "commandLaw") {
                int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);
                string law = cmd_splited[2];

                DEBUG_TRACE("Updating commandLaw num %d to %s", ncontrol,
                            law.c_str());

                if (rtc->d_control[ncontrol]->get_type() == "generic") {
                  sutra_controller_generic *control =
                      dynamic_cast<sutra_controller_generic *>(rtc->d_control[ncontrol]);
                  control->set_commandlaw(law);
                  return;
                } else {
                  BRAMA_DEBUG_TRACE(
                      "controller %d must be a generic controller", ncontrol);
                  throw CORBA::BAD_PARAM();
                }
              } else if (cmd_splited[0] == "PertVolt") {
                int ncontrol = carma_utils::from_string<int>(cmd_splited[1]);
                DEBUG_TRACE("Updating perturbation voltages on controller %d",
                            ncontrol);

                unsigned int ncmd = rtc->d_control[ncontrol]->nactu();
                CORBA::Float *data = (CORBA::Float*) cmd.data.get_buffer();

                if (cmd.dimensions[0] == 1) {
                  if (cmd.dimensions[1] == ncmd) {
                	rtc->lockWFS();
                    rtc->d_control[ncontrol]->set_perturbcom(data, 1);
                	rtc->unlockWFS();
                    return;
                  } else {
                    BRAMA_DEBUG_TRACE("wrong dimensions : %d %d",
                                      cmd.dimensions[0], cmd.dimensions[1]);
                    BRAMA_DEBUG_TRACE("it should be : 1 %d", ncmd);
                    throw CORBA::BAD_PARAM();
                  }
                } else if (cmd.dimensions[0] == 2) {
                  if (cmd.dimensions[1] == ncmd) {
                  	rtc->lockWFS();
                    rtc->d_control[ncontrol]->set_perturbcom(data,
                                                             cmd.dimensions[2]);
                	rtc->unlockWFS();
                    return;
                  } else {
                    BRAMA_DEBUG_TRACE("wrong dimensions : %d %d %d",
                                      cmd.dimensions[0], cmd.dimensions[1],
                                      cmd.dimensions[2]);
                    BRAMA_DEBUG_TRACE("it should be : 2 %d nb_elem", ncmd);
                    throw CORBA::BAD_PARAM();
                  }
                } else {
                  CORBA::ULong *data =
                      (CORBA::ULong*) cmd.dimensions.get_buffer();

                  cerr << "wrong dimensions: ";
                  // Print in Normal order
                  std::copy(
                      data, data + data[0] + 1,
                      std::ostream_iterator<CORBA::ULong>(std::cout, ","));
                  std::cout << "\n";
                  throw CORBA::BAD_PARAM();
                }
              } else {
                throw "Unknown parameter";
              }
            } catch (char const*msg) {
              //BRAMA_DEBUG_TRACE("%s",msg);
            } catch (CORBA::BAD_PARAM &p) {
              std::ostringstream stm;
              stm << p;
              BRAMA_DEBUG_TRACE("%s", stm.str().c_str());
            }
          }
        } else if (status == DDS::RETCODE_NO_DATA) {
//          cerr << "INFO: Event reading complete after " << count << " samples." << endl;
//          cerr << "ERROR: reader received " << "DDS::RETCODE_NO_DATA!" << endl;
          break;
        } else {
          cerr << "ERROR: read BRAMACommand: Error: " << status << endl;
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
void sutra_rtc_bramaListenerImpl::on_requested_deadline_missed(
    DDS::DataReader_ptr reader,
    const DDS::RequestedDeadlineMissedStatus & status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_requested_deadline_missed" << endl;
}
void sutra_rtc_bramaListenerImpl::on_requested_incompatible_qos(
    DDS::DataReader_ptr reader,
    const DDS::RequestedIncompatibleQosStatus & status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_requested_incompatible_qos" << endl;
}
void sutra_rtc_bramaListenerImpl::on_liveliness_changed(
    DDS::DataReader_ptr reader, const DDS::LivelinessChangedStatus & status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_liveliness_changed" << endl;
}
void sutra_rtc_bramaListenerImpl::on_subscription_matched(
    DDS::DataReader_ptr reader, const DDS::SubscriptionMatchedStatus & status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_subscription_matched" << endl;
}
void sutra_rtc_bramaListenerImpl::on_sample_rejected(
    DDS::DataReader_ptr reader, const DDS::SampleRejectedStatus& status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_sample_rejected" << endl;
}
void sutra_rtc_bramaListenerImpl::on_sample_lost(
    DDS::DataReader_ptr reader, const DDS::SampleLostStatus& status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_sample_lost" << endl;
}

#endif /* USE_BRAMA */
