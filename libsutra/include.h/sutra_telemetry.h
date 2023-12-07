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
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_TELEMETRY_H_
#define _SUTRA_TELEMETRY_H_

#include <carma.h>
#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_streams.h>
#include <map>

typedef std::pair<std::string, int32_t> type_telemetry_pair;

class SutraTelemetry {
 protected:
  CarmaStreams *streams;
  std::map<type_telemetry_pair, CarmaHostObj<float> *> objs;

 public:
  SutraTelemetry();
  SutraTelemetry(type_telemetry_pair obj, CarmaHostObj<float> *host_obj,
                  uint32_t nb_streams);
  SutraTelemetry(std::string type_obj, int32_t num_obj,
                  CarmaHostObj<float> *host_obj, uint32_t nb_streams);
  // SutraTelemetry(const SutraTelemetry& src_sutra_telemetry);
  ~SutraTelemetry();

  // const CarmaHostObj<float>& operator[](int32_t idx) const {return
  // get_sutra_host_obj(idx);} ;

  int32_t get_nb_streams();
  int32_t add_stream();
  int32_t add_stream(int32_t nb);
  int32_t del_stream();
  int32_t del_stream(int32_t nb);
  int32_t del_all_streams();
  cudaStream_t get_cuda_stream(int32_t stream);

  int32_t get_nbObjs();
  int32_t add_obj(type_telemetry_pair obj, CarmaHostObj<float> *host_obj);
  int32_t add_obj(std::string type_obj, int32_t num_obj,
              CarmaHostObj<float> *host_obj);
  int32_t del_obj(type_telemetry_pair obj);
  int32_t del_obj(std::string type_obj, int32_t num_obj);
  int32_t del_all_objs();
  CarmaHostObj<float> *get_CarmaHostObj(type_telemetry_pair obj);
  CarmaHostObj<float> *get_CarmaHostObj(std::string type_obj, int32_t num_obj);
  // CarmaHostObj<float>* operator[](int32_t idx) {return
  // get_CarmaHostObj(idx);} ;
  int32_t cpy_obj(std::string type_obj, int32_t num_obj, CarmaObj<float> *d_obj,
              cudaMemcpyKind flag);
  int32_t wait_obj(std::string type_obj, int32_t num_obj);

  int32_t wait_stream(int32_t stream);
  int32_t wait_all_streams();

  /**< Memory transfer */
  int32_t fill_from(std::string type_obj, int32_t num_obj, float *data);
  int32_t fill_into(std::string type_obj, int32_t num_obj, float *data);
};

#endif  // _SUTRA_TELEMETRY_H_
