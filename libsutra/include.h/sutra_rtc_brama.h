/*
 * sutra_rtc_brama.h
 *
 *  Created on: Feb 17, 2015
 *      Author: sevin
 */


#ifndef SUTRA_RTC_BRAMA_H_
#define SUTRA_RTC_BRAMA_H_

#include<sutra_rtc.h>
#include<sutra_wfs.h>
#include<sutra_target.h>
#include<BRAMA_supervisor.h>
#include<sutra_rtc_bramaListenerImpl.h>


class sutra_rtc_brama: public sutra_rtc {
private:
  BRAMA_supervisor *brama;

  DDS::DataReaderListener_var cmd_listener;
  sutra_rtc_bramaListenerImpl* cmd_listener_servant;
  DDS::DataReader_var cmd_dr;

  DDS::DataWriter_var superframe_base_dw;
  BRAMA::SuperFrameDataWriter_var superframe_dw;
  DDS::InstanceHandle_t superframe_handle;

  DDS::DataWriter_var megaframe_base_dw;
  BRAMA::MegaFrameDataWriter_var megaframe_dw;
  DDS::InstanceHandle_t megaframe_handle;

  CORBA::Octet* buff_wfs;
  CORBA::Octet* buff_intensities;
  CORBA::Octet* buff_slopes;
  CORBA::Octet* buff_commands;
  CORBA::Octet* buff_target;

  CORBA::ULong* dims_wfs;
  CORBA::ULong* dims_intensities;
  CORBA::ULong* dims_slopes;
  CORBA::ULong* dims_commands;
  CORBA::ULong* dims_target;

  long framecounter;
  ACE_Mutex lock_;

  int wfs_size;
  sutra_sensors *wfs;
  int target_size;
  sutra_target *target;

  int nslp;
  int ncmd;
  int nvalid;

public:
  sutra_rtc_brama(carma_context *context, sutra_sensors *wfs, sutra_target *target, ACE_TCHAR* name);
  ~sutra_rtc_brama();

  void publish();
private:
  void allocateBuffers();

};

#endif /* SUTRA_RTC_BRAMA_H_ */
