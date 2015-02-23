/*
 * sutra_rtc_brama.h
 *
 *  Created on: Feb 17, 2015
 *      Author: sevin
 */


#ifndef SUTRA_RTC_BRAMA_H_
#define SUTRA_RTC_BRAMA_H_

#include<sutra_rtc.h>
#include<BRAMA_supervisor.h>

class sutra_rtc_brama: public sutra_rtc {
private:
  BRAMA_supervisor *brama;
  DDS::DataReaderListener_var cmd_listener;
  CommandDataReaderListenerImpl* cmd_listener_servant;
  DDS::DataReader_var cmd_dr;
  DDS::DataWriter_var superframe_base_dw;
  BRAMA::SuperFrameDataWriter_var superframe_dw;
  BRAMA::SuperFrame superframe;
  DDS::InstanceHandle_t superframe_handle;

  CORBA::Octet* buff_intensities;
  CORBA::Octet* buff_slopes;
  CORBA::Octet* buff_commands;

  CORBA::ULong* dims_intensities;
  CORBA::ULong* dims_slopes;
  CORBA::ULong* dims_commands;

  long framecounter;

public:
  sutra_rtc_brama(carma_context *context, ACE_TCHAR* name);
  ~sutra_rtc_brama();

  void initDDS();
  void publish();

};

#endif /* SUTRA_RTC_BRAMA_H_ */

