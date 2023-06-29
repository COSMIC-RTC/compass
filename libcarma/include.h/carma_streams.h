// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_streams.h
//! \ingroup   libcarma
//! \class     CarmaStreams
//! \brief     this class provides the stream features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

#ifndef _CARMA_STREAM_H_
#define _CARMA_STREAM_H_

#include <carma_utils.h>
#include <driver_types.h>
#include <vector>

class CarmaStreams {
 protected:
  std::vector<cudaStream_t> streams;
  std::vector<cudaEvent_t> events;
  int eventflags;

 public:
  CarmaStreams();
  CarmaStreams(unsigned int nb_streams);
  // carma_stream(const carma_stream& src_carma_stream);
  ~CarmaStreams();

  cudaStream_t get_stream(int stream);
  cudaEvent_t get_event(int stream);
  cudaStream_t operator[](int idx) { return get_stream(idx); }

  int get_nb_streams();
  int add_stream();
  int add_stream(int nb);
  int del_stream();
  int del_stream(int nb);
  int del_all_streams();
  int wait_event(int stream);
  int wait_stream(int stream);
  int wait_all_streams();
};

#endif  // _CARMA_STREAM_H_
