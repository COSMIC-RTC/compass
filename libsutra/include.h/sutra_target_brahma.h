// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_target_brahma.h
//! \ingroup   libsutra
//! \class     SutraTargetBrahma
//! \brief     this class provides the target_brahma features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
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

  int64_t framecounter;
  int64_t samplecounter;
  int32_t subsample;
  ACE_Mutex lock_;

  int32_t is_initialised;

 public:
  SutraTargetBrahma(CarmaContext *context, ACE_TCHAR *name,
                      SutraTelescope *d_tel, int32_t subsample, int32_t ntargets,
                      float *xpos, float *ypos, float *lambda, float *mag,
                      float zerop, int64_t *sizes, int32_t Npts, int32_t device);
  ~SutraTargetBrahma();

  void set_subsample(int32_t ntarget, int32_t subsample);
  void publish();

 private:
  void allocate_buffers();
};

#endif /* SUTRA_TARGET_BRAHMA_H_ */
