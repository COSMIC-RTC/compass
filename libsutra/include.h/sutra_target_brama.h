/*
 * sutra_target_brama.h
 *
 *  Created on: Feb 17, 2015
 *      Author: sevin
 */

#ifndef SUTRA_TARGET_BRAMA_H_
#define SUTRA_TARGET_BRAMA_H_

#include <BRAMA_context.h>
#include <sutra_target.h>
#include <sutra_target_bramaListenerImpl.h>

class sutra_target_brama : public sutra_target {
private:
  DDS::Subscriber_var sub;
  DDS::Publisher_var pub;
  DDS::DataReaderListener_var cmd_listener;
  sutra_target_bramaListenerImpl *cmd_listener_servant;
  DDS::DataReader_var cmd_dr;
  DDS::DataWriter_var frame_base_dw;
  BRAMA::FrameDataWriter_var frame_dw;
  DDS::InstanceHandle_t frame_handle;

  CORBA::Octet *buff_pixels;

  CORBA::ULong *dims_pixels;

  long framecounter;
  long samplecounter;
  int subsample;
  ACE_Mutex lock_;

  int is_initialised;

public:
  sutra_target_brama(carma_context *context, ACE_TCHAR *name,
                     sutra_telescope *d_tel, int subsample, int ntargets,
                     float *xpos, float *ypos, float *lambda, float *mag,
                     float zerop, long *sizes, int Npts, int device);
  ~sutra_target_brama();

  void set_subsample(int ntarget, int subsample);
  void publish();

private:
  void allocateBuffers();
};

#endif /* SUTRA_TARGET_BRAMA_H_ */
