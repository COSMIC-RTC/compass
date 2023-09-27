// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_template.cpp
//! \ingroup   libsutra
//! \class     SutraTemplate
//! \brief     this class provides a class template to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_template.h>

SutraTemplate::SutraTemplate(CarmaContext *context, const char *type,
                                   long dim, int device) {
  // some inits
  this->current_context = context;
  this->dim = dim;
  this->type = type;
  this->device = device;

  // allocate data and result
  long *dims_data1 = new long[2];
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

int SutraTemplate::fill_data(float *idata) {
  // fill data with an external array
  this->d_data->host2device(idata);

  return EXIT_SUCCESS;
}

int SutraTemplate::fill_data() {
  // fill data with random numbers
  this->d_data->init_prng();
  this->d_data->prng('N');

  return EXIT_SUCCESS;
}

int SutraTemplate::do_compute() {
  // do computation on data and store in result
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(current_context->get_device(device), this->dim,
                         nb_blocks, nb_threads);

  comp_aotemplate(nb_threads, nb_blocks, this->d_data->get_data(),
                  this->d_res->get_data(), this->d_data->get_nb_elements());

  return EXIT_SUCCESS;
}
