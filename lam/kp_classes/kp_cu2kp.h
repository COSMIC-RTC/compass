
#ifndef __KP_CU2KP_H__
#define __KP_CU2KP_H__

#include <carma_obj.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <typeinfo>
#include <vector>
#include "cublas_v2.h"
#include "kernels.h"
#include "kp_KFPP.h"
#include "kp_cu_cublas.h"
#include "kp_cu_cuda.h"

using namespace std;

#include "kp_cu_smatrix.h"
#include "kp_smatrix.h"
#include "kp_vector.h"

template <class real>
class kp_vector;
template <class real>
class kp_matrix;
template <class real>
class kp_smatrix;
template <class real>
class kp_cu_vector;
template <class real>
class kp_cu_matrix;
template <class real>
class kp_cu_smatrix;
// class kp_vector_double;
// class kp_matrix_double;
// class kp_smatrix_double;
// class kp_cu_vector_double;
// class kp_cu_matrix_double;
// class kp_cu_smatrix_double;

template <typename T, typename U>
void kp_cu2kp_smatrix(kp_smatrix<T>& A, const kp_cu_smatrix<U>& M);

template <typename T, typename U>
void kp_cu2kp_matrix(kp_matrix<T>& A, const kp_cu_matrix<U>& M);

template <typename T, typename U>
void kp_cu2kp_vector(kp_vector<T>& u, const kp_cu_vector<U>& v);

template <typename real>
void init_from_kp_cu_vector2kp_vector(kp_vector<real>& u,
                                      const kp_cu_vector<real>& v,
                                      int ind_begin, int ind_end);

#include "kp_cu2kp.hpp"

#endif
