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
//! \version   5.5.0
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
  int32_t eventflags;

 public:
  CarmaStreams();
  CarmaStreams(uint32_t nb_streams);
  // carma_stream(const carma_stream& src_carma_stream);
  ~CarmaStreams();

  cudaStream_t get_stream(int32_t stream);
  cudaEvent_t get_event(int32_t stream);
  cudaStream_t operator[](int32_t idx) { return get_stream(idx); }

  int32_t get_nb_streams();
  int32_t add_stream();
  int32_t add_stream(int32_t nb);
  int32_t del_stream();
  int32_t del_stream(int32_t nb);
  int32_t del_all_streams();
  int32_t wait_event(int32_t stream);
  int32_t wait_stream(int32_t stream);
  int32_t wait_all_streams();
};

#endif  // _CARMA_STREAM_H_
