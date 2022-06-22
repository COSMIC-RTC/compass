// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_kl.cpp
//! \ingroup   libsutra
//! \class     SutraKL
//! \brief     this class provides the kl features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <sutra_kl.h>

SutraKL::SutraKL(CarmaContext *context, long dim, long nr, long np, long nkl,
                   long nord, int device) {
  // some inits
  this->current_context = context;
  this->dim = dim;
  this->nr = nr;
  this->np = np;
  this->nkl = nkl;
  this->nord = nord;
  this->device = device;
  current_context->set_active_device(device, 1);

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = nkl;
  this->h_ord = new CarmaHostObj<int>(dims_data1, MA_PAGELOCK);
  this->d_ord = new CarmaObj<int>(context, dims_data1);

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = dim;
  dims_data2[2] = dim;
  this->d_cr = new CarmaObj<float>(context, dims_data2);
  this->d_cp = new CarmaObj<float>(context, dims_data2);

  dims_data2[1] = nr;
  dims_data2[2] = nkl;
  this->d_rabas = new CarmaObj<float>(context, dims_data2);

  dims_data2[1] = np;
  dims_data2[2] = nord + 1;
  this->d_azbas = new CarmaObj<float>(context, dims_data2);

  // delete[] dims_data1;
  // delete[] dims_data2;

  // Florian features
  // long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = nkl;
  this->d_evals = new CarmaObj<float>(context, dims_data1);
  this->d_ord = new CarmaObj<int>(context, dims_data1);

  // long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = dim;
  dims_data2[2] = dim;
  this->d_covmat = new CarmaObj<float>(context, dims_data2);
  this->d_filter = new CarmaObj<float>(context, dims_data2);
  dims_data2[1] = dim;
  dims_data2[2] = nkl;
  this->d_bas = new CarmaObj<float>(context, dims_data2);

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

SutraKL::~SutraKL() {
  // delete current_context;
  current_context->set_active_device(device, 1);

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

int SutraKL::do_compute(float alpha, float ampli, float *odata, int nkl,
                         int size, int xoff, int yoff) {
  current_context->set_active_device(device, 1);
  // do computation on data and store in result
  int nord = this->h_ord->get_data()[nkl] - 1;

  getkl(alpha, ampli, odata, &(this->d_rabas->get_data()[(nkl) * this->nr]),
        &(this->d_azbas->get_data()[nord * this->np]), this->d_cr->get_data(),
        this->d_cp->get_data(), this->nr, this->np, this->dim, size, xoff, yoff);

  return EXIT_SUCCESS;
}

int SutraKL::do_compute(float ampli, float *odata, int nkl, int size, int xoff,
                         int yoff) {
  return do_compute(0.0f, ampli, odata, nkl, size, xoff, yoff);
}

int SutraKL::do_compute(float *odata, int nkl, int size, int xoff, int yoff) {
  return do_compute(0.0f, 1.0f, odata, nkl, size, xoff, yoff);
}

int SutraKL::do_combi(float *com, float *odata, int size, int xoff, int yoff) {
  current_context->set_active_device(device, 1);
  // do computation on data and store in result
  combikl(com, this->nkl, odata, this->d_rabas->get_data(),
          this->d_ord->get_data(), this->d_azbas->get_data(),
          this->d_cr->get_data(), this->d_cp->get_data(), this->nr, this->np,
          this->dim, size, xoff, yoff);

  return EXIT_SUCCESS;
}

// Florian features
int SutraKL::get_flokl() {
  current_context->set_active_device(device, 1);
  std::cout << "flag in function" << std::endl;
  cget_flokl(this->nkl, this->dim, this->d_covmat->get_data(),
             this->d_filter->get_data(), this->d_bas->get_data());
  return EXIT_SUCCESS;
}
