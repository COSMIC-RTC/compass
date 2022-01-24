// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      carma_streams.cpp
//! \ingroup   libcarma
//! \class     CarmaStreams
//! \brief     this class provides the stream features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.1
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <carma_streams.h>
#include <carma_utils.h>

CarmaStreams::CarmaStreams() {
  this->eventflags = 0;
  CarmaStreams(0);
}

CarmaStreams::CarmaStreams(unsigned int nb_streams) {
  // this->streams  = new vector<cudaStream_t>();
  for (unsigned int i = 0; i < nb_streams; i++) {
    add_stream();
  }

  // carma_safe_call(cudaEventCreateWithFlags(&(this->start_event),
  // cudaDeviceBlockingSync));
  // carma_safe_call(cudaEventCreateWithFlags(&(this->stop_event),
  // cudaDeviceBlockingSync)); cudaEventDefault
}

/*
 carma_stream::carma_stream(const carma_stream& src_carma_stream)
 {
 int nb_streams=src_carma_stream.get_nb_streams();
 carma_stream(nb_streams);
 }
 */

CarmaStreams::~CarmaStreams() {
  del_all_streams();
  this->streams.clear();

  // cudaEventDestroy(this->start_event);
  // cudaEventDestroy(this->stop_event);
}

int CarmaStreams::get_nb_streams() { return this->streams.size(); }

int CarmaStreams::add_stream() {
  cudaStream_t stream_tmp;

  carma_safe_call(cudaStreamCreate(&stream_tmp));
  this->streams.push_back(stream_tmp);

  cudaEvent_t event_tmp;
  carma_safe_call(cudaEventCreate(&event_tmp));
  this->events.push_back(event_tmp);

#if DEBUG
  printf("CARMA Stream created @ 0x%p\n", stream_tmp);
#endif

  return get_nb_streams();
}

int CarmaStreams::add_stream(int nb) {
  for (int stream = 0; stream < nb; stream++) add_stream();
  return get_nb_streams();
}

int CarmaStreams::del_stream() {
  if (streams.empty()) return 0;

#if DEBUG
  printf("CARMA Stream deleting @ 0x%p\n", this->streams.back());
#endif
  carma_safe_call(cudaStreamDestroy(this->streams.back()));
  this->streams.pop_back();

  carma_safe_call(cudaEventDestroy(this->events.back()));
  this->events.pop_back();

  return get_nb_streams();
}

int CarmaStreams::del_stream(int nb) {
  for (int stream = 0; stream < nb && !streams.empty(); stream++) del_stream();
  return get_nb_streams();
}

int CarmaStreams::del_all_streams() {
  while (!streams.empty()) del_stream();
  return get_nb_streams();
}

cudaStream_t CarmaStreams::get_stream(int stream) {
  return this->streams[stream];
}

cudaEvent_t CarmaStreams::get_event(int stream) {
  return this->events[stream];
}

int CarmaStreams::wait_event(int stream) {
  carma_safe_call(cudaEventSynchronize(this->events[stream]));
  return EXIT_SUCCESS;
}

int CarmaStreams::wait_stream(int stream) {
  carma_safe_call(cudaStreamSynchronize(this->streams[stream]));
  return EXIT_SUCCESS;
}

int CarmaStreams::wait_all_streams() {
  for (unsigned int stream = 0; stream < streams.size(); stream++)
    carma_safe_call(cudaStreamSynchronize(this->streams[stream]));
  return EXIT_SUCCESS;
}
