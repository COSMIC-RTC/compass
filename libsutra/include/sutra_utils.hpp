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

//! \file      sutra_utils.hpp
//! \ingroup   libsutra
//! \class     sutra_utils
//! \brief     this file provides utilities to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_UTILS_H_
#define _SUTRA_UTILS_H_
#include <carma.hpp>
#include <carma_obj.hpp>
#include <carma_sparse_obj.hpp>

int32_t compute_nmaxhr(int64_t nvalid);
int32_t cfillrealp(cuFloatComplex *d_odata, float *d_idata, int32_t N,
               CarmaDevice *device);
int32_t cgetrealp(float *d_odata, cuFloatComplex *d_idata, int32_t N,
              CarmaDevice *device);
int32_t abs2(float *d_odata, cuFloatComplex *d_idata, int32_t N, CarmaDevice *device);
int32_t abs2(float *d_odata, cuFloatComplex *d_idata, int32_t N, float fact,
         CarmaDevice *device);
int32_t abs2c(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t N,
          CarmaDevice *device);
int32_t convolve(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t N,
             CarmaDevice *device);
int32_t convolve_modulate(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t mod,
                      int32_t N, CarmaDevice *device);
int32_t subap_norm(float *d_odata, float *d_idata, float *fact, float *norm,
               float nphot, int32_t n, int32_t N, CarmaDevice *device);
int32_t fillindx(float *d_odata, float *d_idata, int32_t *indx, float alpha, float beta,
             int32_t N, CarmaDevice *device);
int32_t fillindx(float *d_odata, float *d_idata, int32_t *indx, float alpha, int32_t N,
             CarmaDevice *device);
int32_t fillindx(float *d_odata, float *d_idata, int32_t *indx, int32_t N,
             CarmaDevice *device);
int32_t fillarr2d(float *d_odata, float *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
              CarmaDevice *device);
int32_t fillarr2d(float *d_odata, float *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
              int32_t dir, CarmaDevice *device);
int32_t getarr2d(float *d_odata, float *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
             CarmaDevice *device);
template <class T>
int32_t addai(T *d_odata, T *i_data, int32_t i, int32_t sgn, int32_t N, CarmaDevice *device);
int32_t subap_norm_async(float *d_odata, float *d_idata, float *fact, float *norm,
                     float nphot, int32_t n, int32_t N, CarmaStreams *streams,
                     CarmaDevice *device);
// templates
template <class T>
int32_t roll(T *idata, int32_t N, int32_t M, int32_t nim, CarmaDevice *device);
template <class T>
int32_t roll(T *idata, int32_t N, int32_t M, CarmaDevice *device);
template <class T>
int32_t roll_mult(T *odata, T *idata, int32_t N, int32_t M, T alpha, CarmaDevice *device);
template <class T>
int32_t sutra_invgene(CarmaObj<T> *imat, CarmaObj<T> *cmat,
                  CarmaObj<T> *eigenvals, CarmaObj<T> *mod2act,
                  CarmaObj<T> *mes2mod, int32_t nfilt);
template <class T>
int32_t remove_avg(T *data, int32_t N, CarmaDevice *device);
template <class T>
int32_t mult_vect(T *d_data, T *scale, int32_t N, CarmaDevice *device);
template <class T>
int32_t mult_vect(T *d_data, T *scale, T gain, int32_t N, CarmaDevice *device);
template <class T>
int32_t mult_vect(T *d_data, T gain, int32_t N, CarmaDevice *device);

int32_t mult_int(float *o_data, float *i_data, float *scale, float gain, int32_t N,
             CarmaDevice *device);
int32_t mult_int(float *o_data, float *i_data, float *scale, float gain, int32_t N,
             CarmaDevice *device, CarmaStreams *streams);
int32_t mult_int(float *o_data, float *i_data, float gain, int32_t N,
             CarmaDevice *device);
int32_t add_md(float *o_matrix, float *i_matrix, float *i_vector, int32_t N,
           CarmaDevice *device);

#endif  // _SUTRA_UTILS_H_
