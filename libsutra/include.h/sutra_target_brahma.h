/*
 * sutra_target_brahma.h
 *
 *  Created on: Feb 17, 2015
 *      Author: sevin
 */

#ifndef SUTRA_TARGET_BRAHMA_H_
#define SUTRA_TARGET_BRAHMA_H_

#include <BRAHMA_context.h>
#include <sutra_target.h>
#include <sutra_target_brahmaListenerImpl.h>

class sutra_target_brahma : public sutra_target {
 private:
  DDS::Subscriber_var sub;
  DDS::Publisher_var pub;
  DDS::DataReaderListener_var cmd_listener;
  sutra_target_brahmaListenerImpl *cmd_listener_servant;
  DDS::DataReader_var cmd_dr;
  DDS::DataWriter_var frame_base_dw;
  BRAHMA::FrameDataWriter_var frame_dw;
  DDS::InstanceHandle_t frame_handle;

  CORBA::Octet *buff_pixels;

  CORBA::ULong *dims_pixels;

  long framecounter;
  long samplecounter;
  int subsample;
  ACE_Mutex lock_;

  int is_initialised;

 public:
  sutra_target_brahma(carma_context *context, ACE_TCHAR *name,
                      sutra_telescope *d_tel, int subsample, int ntargets,
                      float *xpos, float *ypos, float *lambda, float *mag,
                      float zerop, long *sizes, int Npts, int device);
  ~sutra_target_brahma();

  void set_subsample(int ntarget, int subsample);
  void publish();

 private:
  void allocateBuffers();
};

#endif /* SUTRA_TARGET_BRAHMA_H_ */
