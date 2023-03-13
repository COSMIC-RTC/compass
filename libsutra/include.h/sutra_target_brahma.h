// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_target_brahma.h
//! \ingroup   libsutra
//! \class     SutraTargetBrahma
//! \brief     this class provides the target_brahma features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#ifndef SUTRA_TARGET_BRAHMA_H_
#define SUTRA_TARGET_BRAHMA_H_

#include <BRAHMA_context.h>
#include <sutra_target.h>
#include <sutra_target_brahmaListenerImpl.h>

class SutraTargetBrahma : public SutraTarget {
 private:
  DDS::Subscriber_var sub;
  DDS::Publisher_var pub;
  DDS::DataReaderListener_var cmd_listener;
  SutraTargetBrahmaListenerImpl *cmd_listener_servant;
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
  SutraTargetBrahma(CarmaContext *context, ACE_TCHAR *name,
                      SutraTelescope *d_tel, int subsample, int ntargets,
                      float *xpos, float *ypos, float *lambda, float *mag,
                      float zerop, long *sizes, int Npts, int device);
  ~SutraTargetBrahma();

  void set_subsample(int ntarget, int subsample);
  void publish();

 private:
  void allocate_buffers();
};

#endif /* SUTRA_TARGET_BRAHMA_H_ */
