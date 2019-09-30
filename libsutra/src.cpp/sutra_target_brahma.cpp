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

//! \file      sutra_target_brahma.cpp
//! \ingroup   libsutra
//! \class     sutra_target_brahma
//! \brief     this class provides the target_brahma features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifdef USE_BRAHMA

#include <sutra_target_brahma.h>
#include <sutra_telescope.h>

sutra_target_brahma::sutra_target_brahma(carma_context *context,
                                         ACE_TCHAR *name,
                                         sutra_telescope *d_tel, int subsample_,
                                         int ntargets, float *xpos, float *ypos,
                                         float *lambda, float *mag, float zerop,
                                         long *sizes, int Npts, int device)
    : sutra_target(context, d_tel, ntargets, xpos, ypos, lambda, mag, zerop,
                   sizes, Npts, device) {
  DEBUG_TRACE("init %s", name);
  BRAHMA::BRAHMA_context brahma = BRAHMA::BRAHMA_context::get_instance(name);
  frame_handle = 0;
  framecounter = 0;
  samplecounter = 0;
  subsample = subsample_;
  is_initialised = 0;

  buff_pixels = NULL;
  dims_pixels = NULL;

  string topics[] = BRAHMA_TOPICS;

  if (!brahma.is_initialised()) {
    cerr << "brahma initialisation failed!" << endl;
    //    throw "brahma initialisation failed!";
    return;
  }

  try {
    // Create a subscriber for the command topic
    sub = brahma.create_subscriber();
    // Create a publisher for the frame topic
    pub = brahma.create_publisher();
    /*
        // Create an BRAHMA Command listener
        brahma.register_command_type(topics[BRAHMA::CommandType]);
        cmd_listener = (new sutra_target_brahmaListenerImpl);
        cmd_listener_servant =
            dynamic_cast<sutra_target_brahmaListenerImpl *>(cmd_listener.in());

        if (CORBA::is_nil(cmd_listener.in())) {
          throw "BRAHMA Command listener is nil.";
        }
        cmd_listener_servant->attach_target(this);

        cmd_dr = brahma.create_datareader(sub, topics[BRAHMA::CommandType],
                                          cmd_listener);
    */
    // Create an BRAHMA Frame writer
    brahma.register_frame_type(topics[BRAHMA::FrameType]);
    frame_base_dw = brahma.create_datawriter(pub, topics[BRAHMA::FrameType]);
    if (CORBA::is_nil(frame_base_dw.in())) {
      cerr << "create_datawriter for " << topics[BRAHMA::FrameType]
           << " failed." << endl;
      return;
    }
    frame_dw = BRAHMA::FrameDataWriter::_narrow(frame_base_dw.in());
    if (CORBA::is_nil(frame_dw.in())) {
      throw "FrameDataWriter could not be narrowed";
    }

    // BRAHMA::Frame zFrame;
    // frame_handle = frame_dw->register_instance(zFrame);

    is_initialised = 1;

  } catch (CORBA::Exception &e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

sutra_target_brahma::~sutra_target_brahma() {
  if (!is_initialised) {
    return;
  }

  if (buff_pixels) BRAHMA::Values::freebuf(buff_pixels);
  if (dims_pixels) BRAHMA::Dims::freebuf(dims_pixels);
}

void sutra_target_brahma::allocateBuffers() {
  if (!is_initialised) {
    return;
  }

  try {
    // TODO : handle targets with different supports...
    dims_pixels = BRAHMA::Dims::allocbuf(3);
    dims_pixels[0] = d_targets.size();
    dims_pixels[1] = d_targets[0]->d_image_le->getDims(1);
    dims_pixels[2] = d_targets[0]->d_image_le->getDims(2);

    buff_pixels = BRAHMA::Values::allocbuf(
        d_targets.size() * d_targets[0]->d_image_le->getNbElem() *
        sizeof(float));
  } catch (CORBA::Exception &e) {
    cerr << "Exception caught in main.cpp:" << endl << e << endl;
    ACE_OS::exit(1);
  }
}

void sutra_target_brahma::set_subsample(int ntarget, int subsample_) {
  if (!is_initialised) {
    return;
  }

  ACE_Guard<ACE_Mutex> guard(this->lock_);
  this->d_targets[ntarget]->reset_strehlmeter();
  this->samplecounter = 0;
  this->subsample = subsample_;
}

void sutra_target_brahma::publish() {
  if (!is_initialised) {
    cerr << "brahma not initialised!" << endl;
    return;
  }

  ACE_Guard<ACE_Mutex> guard(this->lock_);
  if (subsample <= 0 || samplecounter % subsample != 0 || samplecounter == 0) {
    samplecounter++;
    return;
  }

  if (buff_pixels == NULL) allocateBuffers();

  CORBA::Float *buff_pixels_servant = (CORBA::Float *)buff_pixels;

  long idx = 0;
  carma_obj<float> tmp_img(d_targets[0]->current_context,
                           d_targets[0]->d_image_le->getDims());
  for (size_t target = 0; target < d_targets.size(); target++) {
    d_targets[target]->comp_image(0, true);
    float flux = 1.0f;
    // d_targets[target]->zp * powf(10, -0.4 * d_targets[target]->mag);
    roll_mult<float>(tmp_img.getData(),
                     d_targets[target]->d_image_le->getData(),
                     d_targets[target]->d_image_le->getDims(1),
                     d_targets[target]->d_image_le->getDims(2), flux,
                     d_targets[target]->current_context->get_device(
                         d_targets[target]->device));
    tmp_img.device2host(buff_pixels_servant + idx);

    idx += d_targets[target]->d_image_le->getNbElem();
  }

  BRAHMA::Frame zFrame;
  zFrame.framecounter = framecounter;
  zFrame.timestamp = BRAHMA::get_timestamp();
  zFrame.source = CORBA::string_dup("BRAHMA TARGET");
  zFrame.dimensions = BRAHMA::Dims(3, 3, dims_pixels, 0);
  zFrame.data =
      BRAHMA::Values(idx * sizeof(float), idx * sizeof(float), buff_pixels, 0);
  zFrame.datatype = BRAHMA::BRAHMA_float32_t;
  zFrame.sizeofelements = sizeof(float);

  // cout << "Publishing zFrame: " << zFrame.framecounter << endl;
  DDS::ReturnCode_t ret = frame_dw->write(zFrame, frame_handle);

  if (ret != DDS::RETCODE_OK) {
    ACE_ERROR(
        (LM_ERROR, ACE_TEXT("(%P|%t)ERROR: frame write returned %d.\n"), ret));
    return;
  }

  for (size_t target = 0; target < d_targets.size(); target++) {
    d_targets[target]->reset_strehlmeter();
    samplecounter = 0;
  }

  framecounter++;
  // ACE_Time_Value ace_wait(0, 25);
  // ACE_OS::sleep(ace_wait);
}

#endif /* USE_BRAHMA */
