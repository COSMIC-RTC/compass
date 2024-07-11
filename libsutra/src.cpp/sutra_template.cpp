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

//! \file      sutra_template.cpp
//! \ingroup   libsutra
//! \class     SutraTemplate
//! \brief     this class provides a class template to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_template.hpp>

SutraTemplate::SutraTemplate(CarmaContext *context, const char *type,
                                   int64_t dim, int32_t device) {
  // some inits
  this->current_context = context;
  this->dim = dim;
  this->type = type;
  this->device = device;

  // allocate data and result
  int64_t *dims_data1 = new int64_t[2];
  dims_data1[0] = 1;
  dims_data1[1] = dim;
  this->d_data = new CarmaObj<float>(context, dims_data1);
  this->d_res = new CarmaObj<float>(context, dims_data1);
  delete[] dims_data1;
}

SutraTemplate::~SutraTemplate() {
  // delete data and result
  delete this->d_data;
  delete this->d_res;
}

int32_t SutraTemplate::fill_data(float *idata) {
  // fill data with an external array
  this->d_data->host2device(idata);

  return EXIT_SUCCESS;
}

int32_t SutraTemplate::fill_data() {
  // fill data with random numbers
  this->d_data->init_prng();
  this->d_data->prng('N');

  return EXIT_SUCCESS;
}

int32_t SutraTemplate::do_compute() {
  // do computation on data and store in result
  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(current_context->get_device(device), this->dim,
                         nb_blocks, nb_threads);

  comp_aotemplate(nb_threads, nb_blocks, this->d_data->get_data(),
                  this->d_res->get_data(), this->d_data->get_nb_elements());

  return EXIT_SUCCESS;
}
