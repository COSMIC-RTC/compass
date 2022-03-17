// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      sutra_template.cpp
//! \ingroup   libsutra
//! \class     SutraTemplate
//! \brief     this class provides a class template to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.1
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

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
