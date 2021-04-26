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

//! \file      SutraRtc_brahma.h
//! \ingroup   libsutra
//! \class     SutraRtcBrahma
//! \brief     this class provides the rtc_brahma features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef SUTRA_RTC_BRAHMA_H_
#define SUTRA_RTC_BRAHMA_H_

#include <BRAHMA_context.h>
#include <sutra_rtc.h>
#include <SutraRtcBrahmaListenerImpl.h>
#include <sutra_target.h>
#include <sutra_wfs.h>

template <typename T>
class SutraRtcBrahma : public SutraRtc<float, T, float> {
 private:
  DDS::Subscriber_var sub;
  DDS::Publisher_var pub;

  DDS::DataReaderListener_var cmd_listener;
  SutraRtcBrahmaListenerImpl<T> *cmd_listener_servant;
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
  SutraSensors *wfs;
  int target_size;
  int target_phase_size;
  SutraTarget *target;

  int nslp;
  int ncmd;
  int nvalid;

  int is_initialised;

 public:
  SutraRtcBrahma(CarmaContext *context, SutraSensors *wfs,
                   SutraTarget *target, ACE_TCHAR *name);
  ~SutraRtcBrahma();

  void publish();

 private:
  void allocate_buffers();
};

#endif /* SUTRA_RTC_BRAHMA_H_ */
