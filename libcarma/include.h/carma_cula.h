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
 * \version 4.2.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _CARMA_CULA_H_
#define _CARMA_CULA_H_

#include <carma_host_obj.h>
#include <carma_obj.h>

template <class T>
int carma_cula_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals,
                   carma_obj<T> *mod2act, carma_obj<T> *mes2mod);

template <class T_data>
int carma_cula_svd(carma_host_obj<T_data> *imat,
                   carma_host_obj<T_data> *eigenvals,
                   carma_host_obj<T_data> *mod2act,
                   carma_host_obj<T_data> *mes2mod);

#endif  // _CARMA_CULA_H_
