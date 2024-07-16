// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma_streams.hpp
//! \ingroup   libcarma
//! \class     CarmaStreams
//! \brief     this class provides the stream features to CarmaObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _CARMA_STREAM_H_
#define _CARMA_STREAM_H_

#include <carma_utils.hpp>
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
