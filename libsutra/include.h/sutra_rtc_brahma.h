/*
 * sutra_rtc_brahma.h
 *
 *  Created on: Feb 17, 2015
 *      Author: sevin
 */

#ifndef SUTRA_RTC_BRAHMA_H_
#define SUTRA_RTC_BRAHMA_H_

#include <BRAHMA_context.h>
#include <sutra_rtc.h>
#include <sutra_rtc_brahmaListenerImpl.h>
#include <sutra_target.h>
#include <sutra_wfs.h>

class sutra_rtc_brahma : public sutra_rtc {
 private:
  DDS::Subscriber_var sub;
  DDS::Publisher_var pub;

  DDS::DataReaderListener_var cmd_listener;
  sutra_rtc_brahmaListenerImpl *cmd_listener_servant;
  DDS::DataReader_var cmd_dr;

  DDS::DataWriter_var superframe_base_dw;
  BRAHMA::SuperFrameDataWriter_var superframe_dw;
  DDS::InstanceHandle_t superframe_handle;

  DDS::DataWriter_var megaframe_base_dw;
  BRAHMA::MegaFrameDataWriter_var megaframe_dw;
  DDS::InstanceHandle_t megaframe_handle;

  CORBA::Octet *buff_wfs;
  CORBA::Octet *buff_wfs_phase;
  CORBA::Octet *buff_intensities;
  CORBA::Octet *buff_slopes;
  CORBA::Octet *buff_commands;
  CORBA::Octet *buff_target;
  CORBA::Octet *buff_target_phase;

  CORBA::ULong *dims_wfs;
  CORBA::ULong *dims_wfs_phase;
  CORBA::ULong *dims_intensities;
  CORBA::ULong *dims_slopes;
  CORBA::ULong *dims_commands;
  CORBA::ULong *dims_target;
  CORBA::ULong *dims_target_phase;

  long framecounter;
  ACE_Mutex lock_;

  int wfs_size;
  int wfs_phase_size;
  sutra_sensors *wfs;
  int target_size;
  int target_phase_size;
  sutra_target *target;

  int nslp;
  int ncmd;
  int nvalid;

  int is_initialised;

 public:
  sutra_rtc_brahma(carma_context *context, sutra_sensors *wfs,
                   sutra_target *target, ACE_TCHAR *name);
  ~sutra_rtc_brahma();

  void publish();

 private:
  void allocateBuffers();
};

#endif /* SUTRA_RTC_BRAHMA_H_ */
