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

//! \file      sutra_telemetry.h
//! \ingroup   libsutra
//! \class     sutra_telemetry
//! \brief     this class provides the telemetry features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_TELEMETRY_H_
#define _SUTRA_TELEMETRY_H_

#include <carma.h>
#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_streams.h>
#include <map>

typedef std::pair<std::string, int> type_telemetry_pair;

class sutra_telemetry {
 protected:
  carma_streams *streams;
  std::map<type_telemetry_pair, carma_host_obj<float> *> objs;

 public:
  sutra_telemetry();
  sutra_telemetry(type_telemetry_pair obj, carma_host_obj<float> *host_obj,
                  unsigned int nbStreams);
  sutra_telemetry(std::string type_obj, int num_obj,
                  carma_host_obj<float> *host_obj, unsigned int nbStreams);
  // sutra_telemetry(const sutra_telemetry& src_sutra_telemetry);
  ~sutra_telemetry();

  // const carma_host_obj<float>& operator[](int idx) const {return
  // get_sutra_host_obj(idx);} ;

  int get_nbStreams();
  int add_stream();
  int add_stream(int nb);
  int del_stream();
  int del_stream(int nb);
  int del_all_streams();
  cudaStream_t get_cudaStream_t(int stream);

  int get_nbObjs();
  int add_obj(type_telemetry_pair obj, carma_host_obj<float> *host_obj);
  int add_obj(std::string type_obj, int num_obj,
              carma_host_obj<float> *host_obj);
  int del_obj(type_telemetry_pair obj);
  int del_obj(std::string type_obj, int num_obj);
  int del_all_objs();
  carma_host_obj<float> *get_carma_host_obj(type_telemetry_pair obj);
  carma_host_obj<float> *get_carma_host_obj(std::string type_obj, int num_obj);
  // carma_host_obj<float>* operator[](int idx) {return
  // get_carma_host_obj(idx);} ;
  int cpy_obj(std::string type_obj, int num_obj, carma_obj<float> *d_obj,
              cudaMemcpyKind flag);
  int wait_obj(std::string type_obj, int num_obj);

  int wait_stream(int stream);
  int wait_all_streams();

  /**< Memory transfer */
  int fill_from(std::string type_obj, int num_obj, float *data);
  int fill_into(std::string type_obj, int num_obj, float *data);
};

#endif  // _SUTRA_TELEMETRY_H_
