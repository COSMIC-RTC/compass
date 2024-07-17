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

//! \file      sutra_kl.hpp
//! \ingroup   libsutra
//! \class     SutraKL
//! \brief     this class provides the kl features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_KL_H_
#define _SUTRA_KL_H_

#include <carma.hpp>
#include <carma_host_obj.hpp>
#include <carma_obj.hpp>

class SutraKL {
 public:
  int32_t device;  // # device
  int64_t dim;    // dim of final array
  int64_t nr;     // # radial points
  int64_t np;     // # of elements
  int64_t nkl;    // # of functions in the basis
  int64_t nord;   // # number of radial orders

  CarmaObj<float> *d_rabas;   // the radial array of the basis
  CarmaObj<float> *d_azbas;   // the azimuthal array of the basis
  CarmaHostObj<int32_t> *h_ord;  // the radial orders of the basis
  CarmaObj<int32_t> *d_ord;       // the radial orders of the basis
  CarmaObj<float> *d_cr;      //
  CarmaObj<float> *d_cp;      //

  // Florian features
  CarmaObj<float> *d_covmat;
  CarmaObj<float> *d_filter;
  CarmaObj<float> *d_bas;
  CarmaObj<float> *d_evals;

  CarmaContext *current_context;  // the context in which it has been created
 public:
  SutraKL(CarmaContext *context, int64_t dim, int64_t nr, int64_t np, int64_t nkl,
           int64_t nord, int32_t device);
  ~SutraKL();

  int32_t do_compute(float alpha, float ampli, float *odata, int32_t nkl, int32_t size,
                 int32_t xoff, int32_t yoff);
  int32_t do_compute(float ampli, float *odata, int32_t nkl, int32_t size, int32_t xoff,
                 int32_t yoff);
  int32_t do_compute(float *odata, int32_t nkl, int32_t size, int32_t xoff, int32_t yoff);
  int32_t do_combi(float *com, float *odata, int32_t size, int32_t xoff, int32_t yoff);

  // Florian features
  int32_t get_flokl();
};
int32_t getkl(float alpha, float ampli, float *d_odata, float *rabas, float *azbas,
          float *cr, float *cp, int32_t nr, int32_t np, int32_t nx, int32_t Nx, int32_t xoff,
          int32_t yoff);
int32_t getkl(float ampli, float *d_odata, float *rabas, float *azbas, float *cr,
          float *cp, int32_t nr, int32_t np, int32_t nx, int32_t Nx, int32_t xoff, int32_t yoff);
int32_t getkl(float *d_odata, float *rabas, float *azbas, float *cr, float *cp,
          int32_t nr, int32_t np, int32_t nx, int32_t Nx, int32_t xoff, int32_t yoff);
int32_t combikl(float *com, int32_t nkl, float *d_odata, float *rabas, int32_t *h_ord,
            float *azbas, float *cr, float *cp, int32_t nr, int32_t np, int32_t nx, int32_t Nx,
            int32_t xoff, int32_t yoff);
int32_t cget_flokl(int64_t nkl, int64_t dim, float *covmat, float *filter, float *bas);
// template <class T> void comp_kl(int32_t threads, int32_t blocks, T *d_idata, T
// *d_odata, int32_t N);

#endif  // _SUTRA_KL_H_
