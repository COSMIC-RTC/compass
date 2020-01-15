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

//! \file      sutra_kl.cpp
//! \ingroup   libsutra
//! \class     sutra_kl
//! \brief     this class provides the kl features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_kl.h>

sutra_kl::sutra_kl(carma_context *context, long dim, long nr, long np, long nkl,
                   long nord, int device) {
  // some inits
  this->current_context = context;
  this->dim = dim;
  this->nr = nr;
  this->np = np;
  this->nkl = nkl;
  this->nord = nord;
  this->device = device;
  current_context->set_activeDevice(device, 1);

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = nkl;
  this->h_ord = new carma_host_obj<int>(dims_data1, MA_PAGELOCK);
  this->d_ord = new carma_obj<int>(context, dims_data1);

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = dim;
  dims_data2[2] = dim;
  this->d_cr = new carma_obj<float>(context, dims_data2);
  this->d_cp = new carma_obj<float>(context, dims_data2);

  dims_data2[1] = nr;
  dims_data2[2] = nkl;
  this->d_rabas = new carma_obj<float>(context, dims_data2);

  dims_data2[1] = np;
  dims_data2[2] = nord + 1;
  this->d_azbas = new carma_obj<float>(context, dims_data2);

  // delete[] dims_data1;
  // delete[] dims_data2;

  // Florian features
  // long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = nkl;
  this->d_evals = new carma_obj<float>(context, dims_data1);
  this->d_ord = new carma_obj<int>(context, dims_data1);

  // long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = dim;
  dims_data2[2] = dim;
  this->d_covmat = new carma_obj<float>(context, dims_data2);
  this->d_filter = new carma_obj<float>(context, dims_data2);
  dims_data2[1] = dim;
  dims_data2[2] = nkl;
  this->d_bas = new carma_obj<float>(context, dims_data2);

  delete[] dims_data1;
  delete[] dims_data2;

  if (nr == 0 && np == 0) {
    delete this->h_ord;
    delete this->d_ord;
    delete this->d_cp;
    delete this->d_cr;
    delete this->d_rabas;
    delete this->d_azbas;
  }
}

sutra_kl::~sutra_kl() {
  // delete current_context;
  current_context->set_activeDevice(device, 1);

  delete this->h_ord;
  delete this->d_ord;
  delete this->d_cp;
  delete this->d_cr;
  delete this->d_rabas;
  delete this->d_azbas;

  // Florian features
  delete this->d_covmat;
  delete this->d_filter;
  delete this->d_bas;
}

int sutra_kl::do_compute(float alpha, float ampli, float *odata, int nkl,
                         int size, int xoff, int yoff) {
  current_context->set_activeDevice(device, 1);
  // do computation on data and store in result
  int nord = this->h_ord->getData()[nkl] - 1;

  getkl(alpha, ampli, odata, &(this->d_rabas->getData()[(nkl) * this->nr]),
        &(this->d_azbas->getData()[nord * this->np]), this->d_cr->getData(),
        this->d_cp->getData(), this->nr, this->np, this->dim, size, xoff, yoff);

  return EXIT_SUCCESS;
}

int sutra_kl::do_compute(float ampli, float *odata, int nkl, int size, int xoff,
                         int yoff) {
  return do_compute(0.0f, ampli, odata, nkl, size, xoff, yoff);
}

int sutra_kl::do_compute(float *odata, int nkl, int size, int xoff, int yoff) {
  return do_compute(0.0f, 1.0f, odata, nkl, size, xoff, yoff);
}

int sutra_kl::do_combi(float *com, float *odata, int size, int xoff, int yoff) {
  current_context->set_activeDevice(device, 1);
  // do computation on data and store in result
  combikl(com, this->nkl, odata, this->d_rabas->getData(),
          this->d_ord->getData(), this->d_azbas->getData(),
          this->d_cr->getData(), this->d_cp->getData(), this->nr, this->np,
          this->dim, size, xoff, yoff);

  return EXIT_SUCCESS;
}

// Florian features
int sutra_kl::get_flokl() {
  current_context->set_activeDevice(device, 1);
  std::cout << "flag in function" << std::endl;
  cget_flokl(this->nkl, this->dim, this->d_covmat->getData(),
             this->d_filter->getData(), this->d_bas->getData());
  return EXIT_SUCCESS;
}
