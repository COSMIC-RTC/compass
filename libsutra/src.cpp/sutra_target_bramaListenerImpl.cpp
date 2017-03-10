#ifdef USE_BRAMA

#include "libBRAMATypeSupportC.h"
#include "libBRAMATypeSupportImpl.h"
#include<sutra_target_bramaListenerImpl.h>
#include "dds/DCPS/Service_Participant.h"
#include "ace/streams.h"

#include "sutra_target_brama.h"

//Constructor
sutra_target_bramaListenerImpl::sutra_target_bramaListenerImpl(): target(0L) {

}

//Destructor
sutra_target_bramaListenerImpl::~sutra_target_bramaListenerImpl(void) {
}

// app-specific
void sutra_target_bramaListenerImpl::attach_target(sutra_target_brama *target_){
  target = target_;
}

void sutra_target_bramaListenerImpl::on_data_available(DDS::DataReader_ptr reader)
    throw (CORBA::SystemException) {
  try {
    BRAMA::CommandDataReader_var target_cmd_dr = BRAMA::CommandDataReader::_narrow(
        reader);

    if (CORBA::is_nil(target_cmd_dr.in())) {
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
                DEBUG_TRACE("Updating canapassExposure on target %d to %d", ntarget, subsample);

                target->set_subsample(ntarget, subsample);

//                if ((cmd.dimensions[0] != 1) && (cmd.dimensions[1] == ncmd)) {
//                  BRAMA_DEBUG_TRACE("wrong dimensions : %d %d",
//                                    cmd.dimensions[0], cmd.dimensions[1]);
//                  BRAMA_DEBUG_TRACE("it should be : 1 %d", ncmd);
//                  throw CORBA::BAD_PARAM();
//                }
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
void sutra_target_bramaListenerImpl::on_requested_deadline_missed(
    DDS::DataReader_ptr reader,
    const DDS::RequestedDeadlineMissedStatus & status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_requested_deadline_missed" << endl;
}
void sutra_target_bramaListenerImpl::on_requested_incompatible_qos(
    DDS::DataReader_ptr reader,
    const DDS::RequestedIncompatibleQosStatus & status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_requested_incompatible_qos" << endl;
}
void sutra_target_bramaListenerImpl::on_liveliness_changed(
    DDS::DataReader_ptr reader, const DDS::LivelinessChangedStatus & status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_liveliness_changed" << endl;
}
void sutra_target_bramaListenerImpl::on_subscription_matched(
    DDS::DataReader_ptr reader, const DDS::SubscriptionMatchedStatus & status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_subscription_matched" << endl;
}
void sutra_target_bramaListenerImpl::on_sample_rejected(
    DDS::DataReader_ptr reader, const DDS::SampleRejectedStatus& status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_sample_rejected" << endl;
}
void sutra_target_bramaListenerImpl::on_sample_lost(
    DDS::DataReader_ptr reader, const DDS::SampleLostStatus& status)
        throw (CORBA::SystemException) {
//  cerr << "CommandDataReaderListenerImpl::on_sample_lost" << endl;
}

#endif /* USE_BRAMA */
