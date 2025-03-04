// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma_streams.cpp
//! \ingroup   libcarma
//! \class     CarmaStreams
//! \brief     this class provides the stream features to CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <carma_streams.hpp>
#include <carma_utils.hpp>

CarmaStreams::CarmaStreams() {
  this->eventflags = 0;
  CarmaStreams(0);
}

CarmaStreams::CarmaStreams(uint32_t nb_streams) {
  // this->streams  = new vector<cudaStream_t>();
  for (uint32_t i = 0; i < nb_streams; i++) {
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
 int32_t nb_streams=src_carma_stream.get_nb_streams();
 carma_stream(nb_streams);
 }
 */

CarmaStreams::~CarmaStreams() {
  del_all_streams();
  this->streams.clear();

  // cudaEventDestroy(this->start_event);
  // cudaEventDestroy(this->stop_event);
}

int32_t CarmaStreams::get_nb_streams() { return this->streams.size(); }

int32_t CarmaStreams::add_stream() {
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

int32_t CarmaStreams::add_stream(int32_t nb) {
  for (int32_t stream = 0; stream < nb; stream++) add_stream();
  return get_nb_streams();
}

int32_t CarmaStreams::del_stream() {
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

int32_t CarmaStreams::del_stream(int32_t nb) {
  for (int32_t stream = 0; stream < nb && !streams.empty(); stream++) del_stream();
  return get_nb_streams();
}

int32_t CarmaStreams::del_all_streams() {
  while (!streams.empty()) del_stream();
  return get_nb_streams();
}

cudaStream_t CarmaStreams::get_stream(int32_t stream) {
  return this->streams[stream];
}

cudaEvent_t CarmaStreams::get_event(int32_t stream) {
  return this->events[stream];
}

int32_t CarmaStreams::wait_event(int32_t stream) {
  carma_safe_call(cudaEventSynchronize(this->events[stream]));
  return EXIT_SUCCESS;
}

int32_t CarmaStreams::wait_stream(int32_t stream) {
  carma_safe_call(cudaStreamSynchronize(this->streams[stream]));
  return EXIT_SUCCESS;
}

int32_t CarmaStreams::wait_all_streams() {
  for (uint32_t stream = 0; stream < streams.size(); stream++)
    carma_safe_call(cudaStreamSynchronize(this->streams[stream]));
  return EXIT_SUCCESS;
}
