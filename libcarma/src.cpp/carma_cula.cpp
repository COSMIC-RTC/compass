/*
 * carma_cula.cpp
 *
 *  Created on: Apr 12, 2012
 *      Author: sevin
 */

#include<carma_obj.h>
#include<carma_host_obj.h>

#ifdef DEBUG
#define CULA_TRACE(fmt, args...) fprintf(stderr, fmt, ## args)
#else
#define CULA_TRACE(fmt, args...) /* */
#endif

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif
#ifndef min
#define min(a,b)  (((a)>(b))?(b):(a))
#endif

#ifdef USE_CULA

#include<cula.hpp>

/** These templates are used to select the proper Iamax executable from T_data*/
template<class T> culaStatus carma_culaDevice_sgesvd(int m, int n, T *mat, T *eigenvals, T *U, T*VT);
/**< Generic template for Iamax executable selection */
template<> culaStatus carma_culaDevice_sgesvd<float>(int m, int n, float *mat, float *eigenvals, float *U, float *VT)
{
  culaStatus info=culaNoError;
  info = culaDeviceSgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
  //magma_sgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  //lapackf77_sgesvd("A", "A", &m, &n, mat, &m, eigenvals, U, &m, VT, &n, h_work, &lwork, &info);

  return info;
}
template<> culaStatus carma_culaDevice_sgesvd<double>(int m, int n, double *mat, double *eigenvals, double *U, double *VT)
{
  culaStatus info=culaNoError;
#ifdef _FULL_CULA
  info = culaDeviceDgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
#else
  printf("you have to buy CULA for using culaDeviceDgesvd (or compile with -DFULL_CULA)\n");
#endif
  //magma_dgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  return info;
}

/** These templates are used to select the proper Iamax executable from T_data*/
template<class T> culaStatus carma_cula_sgesvd(int m, int n, T *mat, T *eigenvals, T *U, T*VT);
/**< Generic template for Iamax executable selection */
template<> culaStatus carma_cula_sgesvd<float>(int m, int n, float *mat, float *eigenvals, float *U, float *VT)
{
  culaStatus info=culaNoError;
  info = culaSgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
  //magma_sgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  //lapackf77_sgesvd("A", "A", &m, &n, mat, &m, eigenvals, U, &m, VT, &n, h_work, &lwork, &info);

  return info;
}
template<> culaStatus carma_cula_sgesvd<double>(int m, int n, double *mat, double *eigenvals, double *U, double *VT)
{
  culaStatus info=culaNoError;
#ifdef _FULL_CULA
  info = culaDgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
#else
  printf("you have to buy CULA for using culaDgesvd (or compile with -DFULL_CULA)\n");
#endif
  //magma_dgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  return info;
}

template <class T> int carma_cula_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals, carma_obj<T> *mod2act, carma_obj<T> *mes2mod)
{

  int n = imat->getDims(2); // number of rows
  int m = imat->getDims(1);// number of cols

    carma_obj<T> tmp(imat);
  carma_culaDevice_sgesvd<T>(m, n, tmp, *eigenvals, *mod2act, *mes2mod);
  return EXIT_SUCCESS;
}

template int carma_cula_svd<float>(caObjS *imat, caObjS *eigenvals, caObjS *mod2act, caObjS *mes2mod);
template int carma_cula_svd<double>(caObjD *imat, caObjD *eigenvals, caObjD *mod2act, caObjD *mes2mod);

template <class T> int carma_cula_svd(carma_host_obj<T> *imat, carma_host_obj<T> *eigenvals, carma_host_obj<T> *mod2act, carma_host_obj<T> *mes2mod)
{
  int n = imat->getDims(2); // number of rows
  int m = imat->getDims(1);// number of cols

  carma_host_obj<T> tmp(imat, MA_PAGELOCK);
  carma_cula_sgesvd<T>(m, n, tmp, *eigenvals, *mes2mod, *mod2act);

  return EXIT_SUCCESS;
}

template int carma_cula_svd<float>(carma_host_obj<float> *imat, carma_host_obj<float> *eigenvals, carma_host_obj<float> *mod2act, carma_host_obj<float> *mes2mod);
template int carma_cula_svd<double>(carma_host_obj<double> *imat, carma_host_obj<double> *eigenvals, carma_host_obj<double> *mod2act, carma_host_obj<double> *mes2mod);
#else
#warning "CULA will not be used"
template<class T>
int carma_cula_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
    carma_obj<T> *mod2act, carma_obj<T> *mes2mod) {
  CULA_TRACE("!!!!!! CULA not used !!!!!!\n");
  return EXIT_FAILURE;
}
template int
carma_cula_svd<float>(carma_obj<float> *imat, carma_obj<float> *eigenvals,
    carma_obj<float> *mod2act, carma_obj<float> *mes2mod);
template int
carma_cula_svd<double>(carma_obj<double> *imat, carma_obj<double> *eigenvals,
    carma_obj<double> *mod2act, carma_obj<double> *mes2mod);

template<class T_data>
int carma_cula_svd(carma_host_obj<T_data> *imat,
    carma_host_obj<T_data> *eigenvals, carma_host_obj<T_data> *mod2act,
    carma_host_obj<T_data> *mes2mod) {
  CULA_TRACE("!!!!!! CULA not used !!!!!!\n");
  return EXIT_FAILURE;
}
template int
carma_cula_svd<float>(carma_host_obj<float> *imat,
    carma_host_obj<float> *eigenvals, carma_host_obj<float> *mod2act,
    carma_host_obj<float> *mes2mod);
template int
carma_cula_svd<double>(carma_host_obj<double> *imat,
    carma_host_obj<double> *eigenvals, carma_host_obj<double> *mod2act,
    carma_host_obj<double> *mes2mod);

#endif
