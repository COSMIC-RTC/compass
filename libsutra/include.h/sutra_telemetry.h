// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_telemetry.h
//! \ingroup   libsutra
//! \class     SutraTelemetry
//! \brief     this class provides the telemetry features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24

#ifndef _SUTRA_TELEMETRY_H_
#define _SUTRA_TELEMETRY_H_

#include <carma.h>
#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_streams.h>
#include <map>

typedef std::pair<std::string, int> type_telemetry_pair;

class SutraTelemetry {
 protected:
  CarmaStreams *streams;
  std::map<type_telemetry_pair, CarmaHostObj<float> *> objs;

 public:
  SutraTelemetry();
  SutraTelemetry(type_telemetry_pair obj, CarmaHostObj<float> *host_obj,
                  unsigned int nb_streams);
  SutraTelemetry(std::string type_obj, int num_obj,
                  CarmaHostObj<float> *host_obj, unsigned int nb_streams);
  // SutraTelemetry(const SutraTelemetry& src_sutra_telemetry);
  ~SutraTelemetry();

  // const CarmaHostObj<float>& operator[](int idx) const {return
  // get_sutra_host_obj(idx);} ;

  int get_nb_streams();
  int add_stream();
  int add_stream(int nb);
  int del_stream();
  int del_stream(int nb);
  int del_all_streams();
  cudaStream_t get_cuda_stream(int stream);

  int get_nbObjs();
  int add_obj(type_telemetry_pair obj, CarmaHostObj<float> *host_obj);
  int add_obj(std::string type_obj, int num_obj,
              CarmaHostObj<float> *host_obj);
  int del_obj(type_telemetry_pair obj);
  int del_obj(std::string type_obj, int num_obj);
  int del_all_objs();
  CarmaHostObj<float> *get_CarmaHostObj(type_telemetry_pair obj);
  CarmaHostObj<float> *get_CarmaHostObj(std::string type_obj, int num_obj);
  // CarmaHostObj<float>* operator[](int idx) {return
  // get_CarmaHostObj(idx);} ;
  int cpy_obj(std::string type_obj, int num_obj, CarmaObj<float> *d_obj,
              cudaMemcpyKind flag);
  int wait_obj(std::string type_obj, int num_obj);

  int wait_stream(int stream);
  int wait_all_streams();

  /**< Memory transfer */
  int fill_from(std::string type_obj, int num_obj, float *data);
  int fill_into(std::string type_obj, int num_obj, float *data);
};

#endif  // _SUTRA_TELEMETRY_H_
