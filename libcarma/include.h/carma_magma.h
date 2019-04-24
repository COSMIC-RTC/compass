/**
 * \file carma_magma.h
 *
 * \class carma_magma
 *
 * \ingroup libcarma
 *
 * \brief this class provides wrappers to the magma functions
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 1.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _CARMA_MAGMA_H_
#define _CARMA_MAGMA_H_

#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>

// MAGMA functions
int magma_disabled();
// template <class T>
// int carma_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
//               carma_obj<T> *mod2act, carma_obj<T> *mes2mod);
template <class T>
int carma_syevd(char jobz, carma_obj<T> *mat, carma_host_obj<T> *eigenvals);
// template <class T, int method>
// int carma_syevd(char jobz, carma_obj<T> *mat, carma_host_obj<T> *eigenvals);
template <class T>
int carma_syevd_m(long ngpu, char jobz, long N, T *mat, T *eigenvals);
template <class T>
int carma_syevd_m(long ngpu, char jobz, carma_host_obj<T> *mat,
                  carma_host_obj<T> *eigenvals);
template <class T>
int carma_syevd_m(long ngpu, char jobz, carma_host_obj<T> *mat,
                  carma_host_obj<T> *eigenvals, carma_host_obj<T> *U);
template <class T>
int carma_getri(carma_obj<T> *d_iA);
template <class T>
int carma_potri(carma_obj<T> *d_iA);
template <class T>
int carma_potri_m(long num_gpus, carma_host_obj<T> *h_A, carma_obj<T> *d_iA);

// MAGMA functions (direct access)
template <class T>
int carma_syevd(char jobz, long N, T *mat, T *eigenvals);
// template <class T, int method>
// int carma_syevd(char jobz, long N, T *mat, T *eigenvals);
template <class T>
int carma_syevd_m(long ngpu, char jobz, long N, T *mat, T *eigenvals);
// template <class T>
// int carma_potri_m(long num_gpus, long N, T *h_A, T *d_iA);

// CULA functions
template <class T>
int carma_cula_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
                   carma_obj<T> *mod2act, carma_obj<T> *mes2mod);

#ifdef CAN_DO_HALF
int custom_half_axpy(half alpha, half *source, int incx, int incy, int N,
                     half *dest, carma_device *device);
#endif

#endif  // _CARMA_MAGMA_H_
