// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_template.hpp
//! \ingroup   libsutra
//! \class     SutraTemplate
//! \brief     this class provides a class template to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_AOTEMPLATE_H_
#define _SUTRA_AOTEMPLATE_H_

#include <sutra_wfs.hpp>

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
