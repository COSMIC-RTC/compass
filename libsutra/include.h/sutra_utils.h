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

//! \file      sutra_utils.h
//! \ingroup   libsutra
//! \class     sutra_utils
//! \brief     this file provides utilities to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_UTILS_H_
#define _SUTRA_UTILS_H_
#include <carma.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>

int compute_nmaxhr(long nvalid);
int cfillrealp(cuFloatComplex *d_odata, float *d_idata, int N,
               CarmaDevice *device);
int cgetrealp(float *d_odata, cuFloatComplex *d_idata, int N,
              CarmaDevice *device);
int abs2(float *d_odata, cuFloatComplex *d_idata, int N, CarmaDevice *device);
int abs2(float *d_odata, cuFloatComplex *d_idata, int N, float fact,
         CarmaDevice *device);
int abs2c(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
          CarmaDevice *device);
int convolve(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
             CarmaDevice *device);
int convolve_modulate(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int mod,
                      int N, CarmaDevice *device);
int subap_norm(float *d_odata, float *d_idata, float *fact, float *norm,
               float nphot, int n, int N, CarmaDevice *device);
int fillindx(float *d_odata, float *d_idata, int *indx, float alpha, float beta,
             int N, CarmaDevice *device);
int fillindx(float *d_odata, float *d_idata, int *indx, float alpha, int N,
             CarmaDevice *device);
int fillindx(float *d_odata, float *d_idata, int *indx, int N,
             CarmaDevice *device);
int fillarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
              CarmaDevice *device);
int fillarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
              int dir, CarmaDevice *device);
int getarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
             CarmaDevice *device);
template <class T>
int addai(T *d_odata, T *i_data, int i, int sgn, int N, CarmaDevice *device);
int subap_norm_async(float *d_odata, float *d_idata, float *fact, float *norm,
                     float nphot, int n, int N, CarmaStreams *streams,
                     CarmaDevice *device);
// templates
template <class T>
int roll(T *idata, int N, int M, int nim, CarmaDevice *device);
template <class T>
int roll(T *idata, int N, int M, CarmaDevice *device);
template <class T>
int roll_mult(T *odata, T *idata, int N, int M, T alpha, CarmaDevice *device);
template <class T>
int sutra_invgene(CarmaObj<T> *imat, CarmaObj<T> *cmat,
                  CarmaObj<T> *eigenvals, CarmaObj<T> *mod2act,
                  CarmaObj<T> *mes2mod, int nfilt);
template <class T>
int remove_avg(T *data, int N, CarmaDevice *device);
template <class T>
int mult_vect(T *d_data, T *scale, int N, CarmaDevice *device);
template <class T>
int mult_vect(T *d_data, T *scale, T gain, int N, CarmaDevice *device);
template <class T>
int mult_vect(T *d_data, T gain, int N, CarmaDevice *device);

int mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
             CarmaDevice *device);
int mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
             CarmaDevice *device, CarmaStreams *streams);
int mult_int(float *o_data, float *i_data, float gain, int N,
             CarmaDevice *device);
int add_md(float *o_matrix, float *i_matrix, float *i_vector, int N,
           CarmaDevice *device);

#endif  // _SUTRA_UTILS_H_
