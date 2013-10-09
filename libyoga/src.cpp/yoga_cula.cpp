/*
 * yoga_cula.cpp
 *
 *  Created on: Apr 12, 2012
 *      Author: sevin
 */

#include<cula.hpp>
#include<yoga_host_obj.h>

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif
#ifndef min
#define min(a,b)  (((a)>(b))?(b):(a))
#endif


/** These templates are used to select the proper Iamax executable from T_data*/
template<class T> culaStatus yoga_culaDevice_sgesvd(int m, int n, T *mat, T *eigenvals, T *U, T*VT);
/**< Generic template for Iamax executable selection */
template<> culaStatus yoga_culaDevice_sgesvd<float>(int m, int n, float *mat, float *eigenvals, float *U, float *VT)
{
	culaStatus  info=culaNoError;
	info = culaDeviceSgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
  //magma_sgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  //lapackf77_sgesvd("A", "A", &m, &n, mat, &m, eigenvals, U, &m, VT, &n, h_work, &lwork, &info);

  return info;
}
template<> culaStatus yoga_culaDevice_sgesvd<double>(int m, int n, double *mat, double *eigenvals, double *U, double *VT)
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
template<class T> culaStatus yoga_cula_sgesvd(int m, int n, T *mat, T *eigenvals, T *U, T*VT);
/**< Generic template for Iamax executable selection */
template<> culaStatus yoga_cula_sgesvd<float>(int m, int n, float *mat, float *eigenvals, float *U, float *VT)
{
	culaStatus  info=culaNoError;
    info = culaSgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n);
  //magma_sgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  //lapackf77_sgesvd("A", "A", &m, &n, mat, &m, eigenvals, U, &m, VT, &n, h_work, &lwork, &info);

  return info;
}
template<> culaStatus yoga_cula_sgesvd<double>(int m, int n, double *mat, double *eigenvals, double *U, double *VT)
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

template <class T> int yoga_cula_svd(yoga_obj<T> *imat, yoga_obj<T> *eigenvals, yoga_obj<T> *mod2act, yoga_obj<T> *mes2mod)
{

  int m = imat->getDims(1);
  int n = imat->getDims(2);

  yoga_obj<T> *tmp = new yoga_obj<T>(imat);
  yoga_culaDevice_sgesvd(m, n, tmp->getData(), eigenvals->getData(), mod2act->getData(), mes2mod->getData());
  delete tmp;

  /* Fonctionne pas...
  yoga_obj<T> tmp(imat);
  yoga_culaDevice_sgesvd(m, n, tmp.getData(), eigenvals->getData(), mod2act->getData(), mes2mod->getData());
  */
  return EXIT_SUCCESS;
}

template int yoga_cula_svd<float>(yObjS *imat, yObjS *eigenvals, yObjS *mod2act, yObjS *mes2mod);
template int yoga_cula_svd<double>(yObjD *imat, yObjD *eigenvals, yObjD *mod2act, yObjD *mes2mod);


template <class T> int yoga_cula_svd(yoga_host_obj<T> *imat, yoga_host_obj<T> *eigenvals, yoga_host_obj<T> *mod2act, yoga_host_obj<T> *mes2mod)
{
  int m = imat->getDims(1);
  int n = imat->getDims(2);

  yoga_host_obj<T> *tmp = new yoga_host_obj<T>(imat, MA_PAGELOCK);
  yoga_cula_sgesvd(m, n, tmp->getData(), eigenvals->getData(), mes2mod->getData(), mod2act->getData());
  delete tmp;

  /* Fonctionne pas...
  yoga_host_obj<T> tmp(imat, MA_PAGELOCK);
  yoga_cula_sgesvd(m, n, tmp.getData(), eigenvals->getData(), mes2mod->getData(), mod2act->getData());
  */

  return EXIT_SUCCESS;
}

template int yoga_cula_svd<float>(yoga_host_obj<float> *imat, yoga_host_obj<float> *eigenvals, yoga_host_obj<float> *mod2act, yoga_host_obj<float> *mes2mod);
template int yoga_cula_svd<double>(yoga_host_obj<double> *imat, yoga_host_obj<double> *eigenvals, yoga_host_obj<double> *mod2act, yoga_host_obj<double> *mes2mod);
