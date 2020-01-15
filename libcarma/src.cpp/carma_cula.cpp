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

//! \file      carma_cula.cpp
//! \ingroup   libcarma
//! \brief     this file provides wrappers to the cula functions
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <carma_cula.h>

#ifdef DEBUG
#define CULA_TRACE(fmt, args...) fprintf(stderr, fmt, ##args)
#else
#define CULA_TRACE(fmt, args...) /* */
#endif

#ifndef max
#define max(a, b) (((a) < (b)) ? (b) : (a))
#endif
#ifndef min
#define min(a, b) (((a) > (b)) ? (b) : (a))
#endif

#ifdef USE_CULA

#include <cula.hpp>

/** These templates are used to select the proper Iamax executable from T_data*/
template <class T>
culaStatus carma_culaDevice_sgesvd(int m, int n, T *mat, T *eigenvals, T *U,
                                   T *VT);
/**< Generic template for Iamax executable selection */
template <>
culaStatus carma_culaDevice_sgesvd<float>(int m, int n, float *mat,
                                          float *eigenvals, float *U,
                                          float *VT) {
  return culaDeviceSgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
  // magma_sgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork,
  // &info); lapackf77_sgesvd("A", "A", &m, &n, mat, &m, eigenvals, U, &m, VT,
  // &n, h_work, &lwork, &info);
}
template <>
culaStatus carma_culaDevice_sgesvd<double>(int m, int n, double *mat,
                                           double *eigenvals, double *U,
                                           double *VT) {
#ifdef _FULL_CULA
  return culaDeviceDgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
#else
  printf(
      "you have to buy CULA for using culaDeviceDgesvd (or compile with "
      "-DFULL_CULA)\n");
  return culaNoError;
#endif
  // magma_dgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork,
  // &info);
}

/** These templates are used to select the proper Iamax executable from T_data*/
template <class T>
culaStatus carma_cula_sgesvd(int m, int n, T *mat, T *eigenvals, T *U, T *VT);
/**< Generic template for Iamax executable selection */
template <>
culaStatus carma_cula_sgesvd<float>(int m, int n, float *mat, float *eigenvals,
                                    float *U, float *VT) {
  return culaSgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
  // magma_sgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork,
  // &info); lapackf77_sgesvd("A", "A", &m, &n, mat, &m, eigenvals, U, &m, VT,
  // &n, h_work, &lwork, &info);
}
template <>
culaStatus carma_cula_sgesvd<double>(int m, int n, double *mat,
                                     double *eigenvals, double *U, double *VT) {
#ifdef _FULL_CULA
  return culaDgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
#else
  printf(
      "you have to buy CULA for using culaDgesvd (or compile with "
      "-DFULL_CULA)\n");
  return culaNoError;
#endif
  // magma_dgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork,
  // &info);
}

template <class T>
int carma_cula_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
                   carma_obj<T> *mod2act, carma_obj<T> *mes2mod) {
  int n = imat->getDims(2);  // number of rows
  int m = imat->getDims(1);  // number of cols

  carma_obj<T> tmp(imat);
  carma_culaDevice_sgesvd<T>(m, n, tmp, *eigenvals, *mod2act, *mes2mod);
  return EXIT_SUCCESS;
}

template int carma_cula_svd<float>(caObjS *imat, caObjS *eigenvals,
                                   caObjS *mod2act, caObjS *mes2mod);
template int carma_cula_svd<double>(caObjD *imat, caObjD *eigenvals,
                                    caObjD *mod2act, caObjD *mes2mod);

template <class T>
int carma_cula_svd(carma_host_obj<T> *imat, carma_host_obj<T> *eigenvals,
                   carma_host_obj<T> *mod2act, carma_host_obj<T> *mes2mod) {
  int n = imat->getDims(2);  // number of rows
  int m = imat->getDims(1);  // number of cols

  carma_host_obj<T> tmp(imat, MA_PAGELOCK);
  carma_cula_sgesvd<T>(m, n, tmp, *eigenvals, *mes2mod, *mod2act);

  return EXIT_SUCCESS;
}

template int carma_cula_svd<float>(carma_host_obj<float> *imat,
                                   carma_host_obj<float> *eigenvals,
                                   carma_host_obj<float> *mod2act,
                                   carma_host_obj<float> *mes2mod);
template int carma_cula_svd<double>(carma_host_obj<double> *imat,
                                    carma_host_obj<double> *eigenvals,
                                    carma_host_obj<double> *mod2act,
                                    carma_host_obj<double> *mes2mod);
#else
// #warning "CULA will not be used"
template <class T>
int carma_cula_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
                   carma_obj<T> *mod2act, carma_obj<T> *mes2mod) {
  CULA_TRACE("!!!!!! CULA not used !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_cula_svd<float>(carma_obj<float> *imat,
                                   carma_obj<float> *eigenvals,
                                   carma_obj<float> *mod2act,
                                   carma_obj<float> *mes2mod);
template int carma_cula_svd<double>(carma_obj<double> *imat,
                                    carma_obj<double> *eigenvals,
                                    carma_obj<double> *mod2act,
                                    carma_obj<double> *mes2mod);

template <class T_data>
int carma_cula_svd(carma_host_obj<T_data> *imat,
                   carma_host_obj<T_data> *eigenvals,
                   carma_host_obj<T_data> *mod2act,
                   carma_host_obj<T_data> *mes2mod) {
  CULA_TRACE("!!!!!! CULA not used !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_cula_svd<float>(carma_host_obj<float> *imat,
                                   carma_host_obj<float> *eigenvals,
                                   carma_host_obj<float> *mod2act,
                                   carma_host_obj<float> *mes2mod);
template int carma_cula_svd<double>(carma_host_obj<double> *imat,
                                    carma_host_obj<double> *eigenvals,
                                    carma_host_obj<double> *mod2act,
                                    carma_host_obj<double> *mes2mod);

#endif
