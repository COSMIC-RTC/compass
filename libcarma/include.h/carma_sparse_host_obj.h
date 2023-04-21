// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_sparse_host_obj.h
//! \ingroup   libcarma
//! \class     CarmaSparseHostObj
//! \brief     this class provides wrappers to the generic carma sparse host object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24


#ifndef CARMASPARSEHOSTOBJ_H_
#define CARMASPARSEHOSTOBJ_H_
#include "carma_host_obj.h"

template <class T_data>
class CarmaSparseObj;

template <class T_data>
class CarmaSparseHostObj {
 public:
  CarmaSparseHostObj();
  CarmaSparseHostObj(CarmaSparseObj<T_data>& sm);
  CarmaSparseHostObj(CarmaSparseHostObj<T_data>& sm);
  CarmaSparseHostObj(const long* dims, T_data* M, char order);
  virtual ~CarmaSparseHostObj();

  // delete all arrays and create new for nnz=new_nnz
  void resize(int new_nnz, int dim1_, int dim2_);

  void operator=(CarmaSparseObj<T_data>& M);
  void operator=(CarmaSparseHostObj<T_data>& M);

  void init_from_matrix(const long* dims, T_data* M, char major_dim);
  void copy_into_matrix(T_data* M, char major_dim);
  void resize2row_major();
  void resize2col_major();
  char get_major_dim() const { return major_dim; }

  /**< General Utilities */
  operator T_data*() { return h_data; }
  T_data* operator[](int index) { return &h_data[index]; }
  T_data* get_data() { return h_data; }
  T_data* get_data(int index) { return &h_data[index]; }
  const long* get_dims() { return dims_data; }
  long get_dims(int i) { return dims_data[i]; }
  int get_nonzero_elem() { return nz_elem; }

 private:
  void _create(int nnz_, int dim1_, int dim2_);  // create new arrays
  void _clear();                                 // clear arrays

  char major_dim;  // U - undefined
  // R - row major
  // C - col major

  MemAlloc malloc_type;  ///< type of host alloc

 public:
  long dims_data[3];  ///< dimensions of the array
  int nz_elem;        ///< number of elements in the array

  // see Sparse Matrix Storage Formats (coordinate format)
  // ONE-BASED!!!!
  T_data* h_data;
  int* rowind;
  int* colind;
};

// Multiply sparce matrix by vector
// y := alpha*A*x + betta * y
template <class T_data>
void carma_gemv(T_data alpha, CarmaSparseHostObj<T_data>* A,
                CarmaHostObj<T_data>* x, T_data betta,
                CarmaHostObj<T_data>* y,
                void (*ptr_coomv)(char* transa, long* m, long* k, T_data* alpha,
                                  char* matdescra, T_data* val, int* rowind,
                                  int* colind, int* nnz, T_data* x,
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
