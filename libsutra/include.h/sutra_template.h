// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_template.h
//! \ingroup   libsutra
//! \class     SutraTemplate
//! \brief     this class provides a class template to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_AOTEMPLATE_H_
#define _SUTRA_AOTEMPLATE_H_

#include <sutra_wfs.h>

class SutraTemplate {
 public:
  int32_t device;        // # device
  std::string type;  // a name for your data
  int64_t dim;          // # of elements

  CarmaObj<float> *d_data;  // the data
  CarmaObj<float> *d_res;   // the result

  CarmaContext *current_context;  // the context in which it has been created

 public:
  SutraTemplate(CarmaContext *context, const char *type, int64_t dim,
                   int32_t device);
  SutraTemplate(const SutraTemplate &aotemplate);
  ~SutraTemplate();

  int32_t fill_data(float *idata);
  int32_t fill_data();
  int32_t do_compute();
};
template <class T>
void comp_aotemplate(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t N);

#endif  // _SUTRA_AOTEMPLATE_H_
