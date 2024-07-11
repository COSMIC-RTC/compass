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

//! \file      sutra_telemetry.cpp
//! \ingroup   libsutra
//! \class     SutraTelemetry
//! \brief     this class provides the telemetry features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_telemetry.hpp>

SutraTelemetry::SutraTelemetry() { streams = new CarmaStreams(); }

SutraTelemetry::SutraTelemetry(type_telemetry_pair obj,
                                 CarmaHostObj<float> *host_obj,
                                 uint32_t nb_streams) {
  add_obj(obj, host_obj);
  streams = new CarmaStreams(nb_streams);
  // carma_safe_call(cudaEventCreateWithFlags(&(this->start_event),
  // cudaDeviceBlockingSync)); carma_safe_call(
  // cudaEventCreateWithFlags(&(this->stop_event), cudaDeviceBlockingSync));
  // cudaEventDefault
}

SutraTelemetry::SutraTelemetry(std::string type_obj, int32_t num_obj,
                                 CarmaHostObj<float> *host_obj,
                                 uint32_t nb_streams) {
  add_obj(type_obj, num_obj, host_obj);
  streams = new CarmaStreams(nb_streams);
  // carma_safe_call(cudaEventCreateWithFlags(&(this->start_event),
  // cudaDeviceBlockingSync)); carma_safe_call(
  // cudaEventCreateWithFlags(&(this->stop_event), cudaDeviceBlockingSync));
  // cudaEventDefault
}

/*
 SutraTelemetry::SutraTelemetry(const SutraTelemetry& src_sutra_telemetry)
 {
 int32_t nb_streams=src_sutra_telemetry.get_nb_streams();
 SutraTelemetry(nb_streams);
 }
 */

SutraTelemetry::~SutraTelemetry() {
  for (std::map<type_telemetry_pair, CarmaHostObj<float> *>::iterator it =
           objs.begin();
       objs.end() != it; it++) {
    delete it->second;
  }
  this->objs.clear();
  delete this->streams;

  // cudaEventDestroy(this->start_event);
  // cudaEventDestroy(this->stop_event);
}

int32_t SutraTelemetry::get_nbObjs() { return this->objs.size(); }

int32_t SutraTelemetry::get_nb_streams() { return this->streams->get_nb_streams(); }

int32_t SutraTelemetry::add_stream() {
  this->streams->add_stream();
  return this->streams->get_nb_streams();
}

int32_t SutraTelemetry::add_stream(int32_t nb) {
  this->streams->add_stream(nb);
  return this->streams->get_nb_streams();
}

int32_t SutraTelemetry::del_stream() {
  this->streams->del_stream();
  return this->streams->get_nb_streams();
}

int32_t SutraTelemetry::del_stream(int32_t nb) {
  this->streams->del_stream(nb);
  return this->streams->get_nb_streams();
}

int32_t SutraTelemetry::del_all_streams() {
  this->streams->del_all_streams();
  return this->streams->get_nb_streams();
}

int32_t SutraTelemetry::add_obj(type_telemetry_pair obj,
                             CarmaHostObj<float> *host_obj) {
  this->objs[obj] = host_obj;
  return get_nbObjs();
}

int32_t SutraTelemetry::add_obj(std::string type_obj, int32_t num_obj,
                             CarmaHostObj<float> *host_obj) {
  this->objs[make_pair(type_obj, num_obj)] = host_obj;
  return get_nbObjs();
}

int32_t SutraTelemetry::del_obj(type_telemetry_pair obj) {
  delete this->objs[obj];
  this->objs.erase(obj);
  return get_nbObjs();
}

int32_t SutraTelemetry::del_obj(std::string type_obj, int32_t num_obj) {
  delete this->objs[make_pair(type_obj, num_obj)];
  this->objs.erase(make_pair(type_obj, num_obj));
  return get_nbObjs();
}

CarmaHostObj<float> *SutraTelemetry::get_CarmaHostObj(
    type_telemetry_pair obj) {
  return this->objs[obj];
}

CarmaHostObj<float> *SutraTelemetry::get_CarmaHostObj(std::string type_obj,
                                                           int32_t num_obj) {
  return this->objs[make_pair(type_obj, num_obj)];
}

int32_t SutraTelemetry::cpy_obj(std::string type_obj, int32_t num_obj,
                             CarmaObj<float> *d_obj, cudaMemcpyKind flag) {
  this->objs[make_pair(type_obj, num_obj)]->cpy_obj(d_obj, flag);
  return EXIT_SUCCESS;
}

/**< Memory transfer */
int32_t SutraTelemetry::fill_from(std::string type_obj, int32_t num_obj, float *data) {
  this->objs[make_pair(type_obj, num_obj)]->fill_from(data);
  return EXIT_SUCCESS;
}

int32_t SutraTelemetry::fill_into(std::string type_obj, int32_t num_obj, float *data) {
  this->objs[make_pair(type_obj, num_obj)]->fill_into(data);
  return EXIT_SUCCESS;
}

int32_t SutraTelemetry::wait_obj(std::string type_obj, int32_t num_obj) {
  this->objs[make_pair(type_obj, num_obj)]->wait_all_streams();
  return EXIT_SUCCESS;
}

cudaStream_t SutraTelemetry::get_cuda_stream(int32_t stream) {
  return this->streams->get_stream(stream);
}

int32_t SutraTelemetry::wait_stream(int32_t stream) {
  this->streams->wait_stream(stream);
  return EXIT_SUCCESS;
}

int32_t SutraTelemetry::wait_all_streams() {
  this->streams->wait_all_streams();
  return EXIT_SUCCESS;
}
