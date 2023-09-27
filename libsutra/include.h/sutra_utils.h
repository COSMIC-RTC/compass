// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_utils.h
//! \ingroup   libsutra
//! \class     sutra_utils
//! \brief     this file provides utilities to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

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
