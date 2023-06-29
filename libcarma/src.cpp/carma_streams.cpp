// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_streams.cpp
//! \ingroup   libcarma
//! \class     CarmaStreams
//! \brief     this class provides the stream features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

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
