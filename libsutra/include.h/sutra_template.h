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
//! \version   5.4.4
//! \date      2022/01/24

#ifndef _SUTRA_AOTEMPLATE_H_
#define _SUTRA_AOTEMPLATE_H_

#include <sutra_wfs.h>

class SutraTemplate {
 public:
  int device;        // # device
  std::string type;  // a name for your data
  long dim;          // # of elements

  CarmaObj<float> *d_data;  // the data
  CarmaObj<float> *d_res;   // the result

  CarmaContext *current_context;  // the context in which it has been created

 public:
  SutraTemplate(CarmaContext *context, const char *type, long dim,
                   int device);
  SutraTemplate(const SutraTemplate &aotemplate);
  ~SutraTemplate();

  int fill_data(float *idata);
  int fill_data();
  int do_compute();
};
template <class T>
void comp_aotemplate(int threads, int blocks, T *d_idata, T *d_odata, int N);

#endif  // _SUTRA_AOTEMPLATE_H_
