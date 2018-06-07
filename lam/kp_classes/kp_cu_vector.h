// kp_cu_vector.h
// Kalman project vector

#ifndef __SEGER__KP_CU_VECTOR_H__
#define __SEGER__KP_CU_VECTOR_H__

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

#include "kp_cu_matrix.h"
#include "kp_matrix.h"
#include "kp_vector.h"

template <typename real>
class kp_cu_matrix;

using namespace std;

template <typename real>
class kp_cu_vector {
  template <typename T>
  friend class kp_cu_vector;

#ifndef KP_WITH_CARMA
 public:
  kp_cu_vector(int s_);

  template <typename T>
  kp_cu_vector(const kp_vector<T>&);
  kp_cu_vector(const int vec_size, real val);

  real* getData() const { return d_cu; }
  real* getData(int idx) const { return &(d_cu[idx]); }
  int size() const { return s; };

  void resize(int s_);

 private:
  void _create(int s_);

#else
 public:
  kp_cu_vector(carma_context* context_);
  kp_cu_vector(carma_context* context_, const kp_vector<real>& v);
  kp_cu_vector(carma_context* context_, const int vec_size, real val);

  real* getData(){return carma_object.getData()};
  real* getData(int idx){return carma_object.getData(idx)};

  int size() const { return carma_object.getDims(1); };

 private:
  carma_obj<real>* carma_object;
  void _create(carma_context* context_, int s_);

#endif

 private:
  real* d_cu;
  int s;
  void _clear();

 public:
  // each method exists twice : one if KP_WITH_CARMA is not defined ; another
  // one if KP_WITH_CARMA is defined
  kp_cu_vector();
  template <typename T>
  kp_cu_vector(const kp_cu_vector<T>&);
  kp_cu_vector(const kp_cu_vector<real>&);
  template <typename T>
  void operator=(const kp_cu_vector<T>& v);
  void operator=(const kp_cu_vector<real>& v);
  template <typename T>
  void operator=(const kp_vector<T>& v);
  void init_from_matrix_column(const kp_cu_matrix<real>& M, int col);
  void init_from_matrix_column(const kp_matrix<real>& M, int col);
  // MATLAB code this(ind_begin + 1, int_end)
  void init_from_vector(const kp_cu_vector<real>&, int ind_begin, int int_end);
  void init_from_vector(const kp_vector<real>&, int ind_begin, int int_end);
  // if KP_WITH_CARMA is defined, handle (1st argument) is not used. This is the
  // handle from carma_object which is used
  void gemv(cublasHandle_t handle, char op_A, real alpha,
            const kp_cu_matrix<real>& A, const kp_cu_vector<real>& x,
            real beta);

  // same methods (do not depend on whether KP_WITH_CARMA is defined or not)
  ~kp_cu_vector();
  void zeros();
  void operator+=(const kp_cu_vector<real>& v);
  void operator-=(const kp_cu_vector<real>& v);
  void operator*=(const kp_cu_vector<real>& v);
  void operator/=(const kp_cu_vector<real>& v);
  void operator+=(real val);
  void operator-=(real val);
  void operator*=(real val);
  void operator/=(real val);
  // in idx we have indexes of v which we should use
  // MATLAB code: this = v(idx)
  void init_from_idx(const kp_vector<real>& v, const vector<int>& idx);
  void init_from_idx(const kp_cu_vector<real>& v, const vector<int>& idx);
  // subv is vector which was created by init_from_rowidx
  // MATLAB code: this(idx) = subv
  void set_from_subvector(const kp_vector<real>& subv, const vector<int>& idx);
  void set_from_subvector(const kp_cu_vector<real>& subv,
                          const vector<int>& idx);
};

template <typename real>
void kp_cu_inverse(kp_cu_vector<real>& v);

template <typename real>
void sqrt(kp_cu_vector<real>& v);

#include "kp_cu_vector.hpp"

#endif
