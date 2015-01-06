//kp_carma_tools.h
//some tools to convert carma_obj <-> kp_matrix, kp_vector


#ifndef __SEGER__KP_CARMA_TOOLS_H__
#define __SEGER__KP_CARMA_TOOLS_H__

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include "kp_KFPP.h"
#include <carma_obj.h>

#include "kp_matrix.h"
#include "kp_vector.h"
#include "kp_cu_vector.h"
#include <typeinfo>

template<typename T, typename U> void kp_carma_host_obj_to_kp_matrix(carma_host_obj<T>& ch, kp_matrix<U>& k);
template<typename T, typename U> void kp_carma_obj_to_kp_matrix(carma_obj<T>& c, kp_matrix<U>& k);
//void kp_carma_obj_to_kp_matrix(carma_obj<int>& c, kp_matrix& k);
template<typename T, typename U> void kp_carma_host_obj_to_kp_vector(carma_host_obj<T>& ch, kp_vector<U>& k);
template<typename T, typename U> void kp_carma_obj_to_kp_vector(carma_obj<T>& c, kp_vector<U>& k);
//void kp_carma_obj_to_kp_vector(carma_obj<int>& c, kp_vector& k);
template<typename T> size_t kp_carma_obj_get_dim1(carma_obj<float>& M);
template<typename T> size_t kp_carma_obj_get_dim2(carma_obj<float>& M);
template<typename T, typename U> void kp_kp_vector_to_carma_obj(const kp_vector<T>& k, carma_obj<U>& c);
template<typename T, typename U> void kp_kp_cu_vector_to_carma_obj(const kp_cu_vector<T>& cu_k, carma_obj<U>& c);
template<typename T, typename U> void kp_carma_obj_to_kp_cu_vector(const carma_obj<T>& c, kp_cu_vector<U>& cu_k);

#include "kp_carma_tools.hpp"

#endif
