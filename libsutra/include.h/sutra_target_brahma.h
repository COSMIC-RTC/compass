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

//! \file      sutra_target_brahma.h
//! \ingroup   libsutra
//! \class     SutraTargetBrahma
//! \brief     this class provides the target_brahma features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
