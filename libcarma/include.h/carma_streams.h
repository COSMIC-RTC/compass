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

//! \file      carma_streams.h
//! \ingroup   libcarma
//! \class     CarmaStreams
//! \brief     this class provides the stream features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
