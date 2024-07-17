// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma_sparse_host_obj.hpp
//! \ingroup   libcarma
//! \class     CarmaSparseHostObj
//! \brief     this class provides wrappers to the generic carma sparse host object
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24


#ifndef CARMASPARSEHOSTOBJ_H_
#define CARMASPARSEHOSTOBJ_H_
#include "carma_host_obj.hpp"

template <class T_data>
class CarmaSparseObj;

template <class T_data>
class CarmaSparseHostObj {
 public:
  CarmaSparseHostObj();
  CarmaSparseHostObj(CarmaSparseObj<T_data>& sm);
  CarmaSparseHostObj(CarmaSparseHostObj<T_data>& sm);
  CarmaSparseHostObj(const int64_t* dims, T_data* M, char order);
  virtual ~CarmaSparseHostObj();

  // delete all arrays and create new for nnz=new_nnz
  void resize(int32_t new_nnz, int32_t dim1_, int32_t dim2_);

  void operator=(CarmaSparseObj<T_data>& M);
  void operator=(CarmaSparseHostObj<T_data>& M);

  void init_from_matrix(const int64_t* dims, T_data* M, char major_dim);
  void copy_into_matrix(T_data* M, char major_dim);
  void resize2row_major();
  void resize2col_major();
  char get_major_dim() const { return major_dim; }

  /**< General Utilities */
  operator T_data*() { return h_data; }
  T_data* operator[](int32_t index) { return &h_data[index]; }
  T_data* get_data() { return h_data; }
  T_data* get_data(int32_t index) { return &h_data[index]; }
  const int64_t* get_dims() { return dims_data; }
  int64_t get_dims(int32_t i) { return dims_data[i]; }
  int32_t get_nonzero_elem() { return nz_elem; }

 private:
  void _create(int32_t nnz_, int32_t dim1_, int32_t dim2_);  // create new arrays
  void _clear();                                 // clear arrays

  char major_dim;  // U - undefined
  // R - row major
  // C - col major

  MemAlloc malloc_type;  ///< type of host alloc

 public:
  int64_t dims_data[3];  ///< dimensions of the array
  int32_t nz_elem;        ///< number of elements in the array

  // see Sparse Matrix Storage Formats (coordinate format)
  // ONE-BASED!!!!
  T_data* h_data;
  int32_t* rowind;
  int32_t* colind;
};

// Multiply sparce matrix by vector
// y := alpha*A*x + betta * y
template <class T_data>
void carma_gemv(T_data alpha, CarmaSparseHostObj<T_data>* A,
                CarmaHostObj<T_data>* x, T_data betta,
                CarmaHostObj<T_data>* y,
                void (*ptr_coomv)(char* transa, int64_t* m, int64_t* k, T_data* alpha,
                                  char* matdescra, T_data* val, int32_t* rowind,
                                  int32_t* colind, int32_t* nnz, T_data* x,
                                  T_data* beta, T_data* y));

// Multiply sparce matrix by dense matrix
// C := alpha*op(A)*B + betta * y
// op_A could be N (do nothing) or T (transpose)
template <class T_data>
void carma_gemm(char op_A, T_data alpha, CarmaSparseHostObj<T_data>* A,
                CarmaHostObj<T_data>* B, T_data betta,
                CarmaHostObj<T_data>* C);

// Multiply sparce matrix by dense matrix
// C := alpha*A*B + betta * y
template <class T_data>
inline void carma_gemm(T_data alpha, CarmaSparseHostObj<T_data>* A,
                       CarmaHostObj<T_data>* B, T_data betta,
                       CarmaHostObj<T_data>* C) {
  kp_gemm('N', alpha, A, B, betta, C);
}

#endif /* CARMASPARSEHOSTOBJ_H_ */
