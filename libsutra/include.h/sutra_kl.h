// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_kl.h
//! \ingroup   libsutra
//! \class     SutraKL
//! \brief     this class provides the kl features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24

#ifndef _SUTRA_KL_H_
#define _SUTRA_KL_H_

#include <carma.h>
#include <carma_host_obj.h>
#include <carma_obj.h>

class SutraKL {
 public:
  int device;  // # device
  long dim;    // dim of final array
  long nr;     // # radial points
  long np;     // # of elements
  long nkl;    // # of functions in the basis
  long nord;   // # number of radial orders

  CarmaObj<float> *d_rabas;   // the radial array of the basis
  CarmaObj<float> *d_azbas;   // the azimuthal array of the basis
  CarmaHostObj<int> *h_ord;  // the radial orders of the basis
  CarmaObj<int> *d_ord;       // the radial orders of the basis
  CarmaObj<float> *d_cr;      //
  CarmaObj<float> *d_cp;      //

  // Florian features
  CarmaObj<float> *d_covmat;
  CarmaObj<float> *d_filter;
  CarmaObj<float> *d_bas;
  CarmaObj<float> *d_evals;

  CarmaContext *current_context;  // the context in which it has been created
 public:
  SutraKL(CarmaContext *context, long dim, long nr, long np, long nkl,
           long nord, int device);
  ~SutraKL();

  int do_compute(float alpha, float ampli, float *odata, int nkl, int size,
                 int xoff, int yoff);
  int do_compute(float ampli, float *odata, int nkl, int size, int xoff,
                 int yoff);
  int do_compute(float *odata, int nkl, int size, int xoff, int yoff);
  int do_combi(float *com, float *odata, int size, int xoff, int yoff);

  // Florian features
  int get_flokl();
};
int getkl(float alpha, float ampli, float *d_odata, float *rabas, float *azbas,
          float *cr, float *cp, int nr, int np, int nx, int Nx, int xoff,
          int yoff);
int getkl(float ampli, float *d_odata, float *rabas, float *azbas, float *cr,
          float *cp, int nr, int np, int nx, int Nx, int xoff, int yoff);
int getkl(float *d_odata, float *rabas, float *azbas, float *cr, float *cp,
          int nr, int np, int nx, int Nx, int xoff, int yoff);
int combikl(float *com, int nkl, float *d_odata, float *rabas, int *h_ord,
            float *azbas, float *cr, float *cp, int nr, int np, int nx, int Nx,
            int xoff, int yoff);
int cget_flokl(long nkl, long dim, float *covmat, float *filter, float *bas);
// template <class T> void comp_kl(int threads, int blocks, T *d_idata, T
// *d_odata, int N);

#endif  // _SUTRA_KL_H_
