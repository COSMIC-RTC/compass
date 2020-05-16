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

//! \file      carma_sparse_obj.h
//! \ingroup   libcarma
//! \class     CarmaSparseObj
//! \brief     this class provides wrappers to the generic carma sparse object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#ifndef CARMA_SPARSE_OBJ_H_
#define CARMA_SPARSE_OBJ_H_
#include <cusparse_v2.h>
#include "carma_obj.h"

template <class T_data>
class CarmaSparseHostObj;

#ifndef USE_MAGMA_SPARSE
typedef void *magma_d_sparse_matrix;
typedef void *magma_s_sparse_matrix;
#else
#include "magmasparse.h"
#endif

template <class T_data>
class CarmaSparseObj {
 public:
  long dims_data[3];  ///< dimensions of the array
  int nz_elem;        ///< number of elements in the array
  int device;         ///< device where the CarmaObj is allocate
  CarmaContext *current_context;

  // ZERO-BASED INDEXING CSR-FORMAT
  T_data *d_data;  // nz_elem elements
  int *d_rowind;   // dim1+1  elements
  int *d_colind;   // nz_elem elements
  cusparseMatDescr_t descr;

  char major_dim;
  std::string format;
  int block_dim;  // blockdim for BSR format

  // Magma stuff
  union {
    magma_d_sparse_matrix d_sparse_mat;
    magma_s_sparse_matrix s_sparse_mat;
  };

 private:
 public:
  CarmaSparseObj(CarmaContext *current_context);
  CarmaSparseObj(CarmaObj<T_data> *M);
  CarmaSparseObj(CarmaSparseObj<T_data> *M);
  CarmaSparseObj(CarmaContext *current_context,
                   CarmaSparseHostObj<T_data> *M);
  CarmaSparseObj(CarmaContext *current_context, const long *dims, T_data *M,
                   bool load_from_host);
  CarmaSparseObj(CarmaContext *current_context, const long *dims,
                   T_data *values, int *colind, int *rowind, int nz,
                   bool load_from_host);
  virtual ~CarmaSparseObj();

  void operator=(CarmaSparseObj<T_data> &M);
  void operator=(CarmaSparseHostObj<T_data> &M);

  void resize(int nnz_, int dim1_, int dim2_);
  void init_from_transpose(CarmaSparseObj<T_data> *M);
  bool is_column_major();
  char get_major_dim() const { return major_dim; }
  void set_majorDim(char c) { major_dim = c; }

  /**< General Utilities */
  operator T_data *() { return d_data; }
  T_data *operator[](int index) { return &d_data[index]; }
  T_data *get_data() { return d_data; }
  T_data *get_data(int index) { return &d_data[index]; }
  const long *get_dims() { return dims_data; }
  long get_dims(int i) { return dims_data[i]; }
  int get_nonzero_elem() { return nz_elem; }
  CarmaContext *get_context() { return current_context; }

  int get_device() { return device; }

  void sparse_to_host(int *h_rowInd, int *h_colInd, T_data *h_data);

 private:
  void _create(int nnz_, int dim1_, int dim2_);
  void _clear();

  template <cusparseStatus_t CUSPARSEAPI (*ptr_nnz)(
                cusparseHandle_t handle, cusparseDirection_t dirA, int m, int n,
                const cusparseMatDescr_t descrA, const T_data *A, int lda,
                int *nnzPerRowCol, int *nnzTotalDevHostPtr),
            cusparseStatus_t CUSPARSEAPI (*ptr_dense2csr)(
                cusparseHandle_t handle, int m, int n,
                const cusparseMatDescr_t descrA, const T_data *A, int lda,
                const int *nnzPerRow, T_data *csrValA, int *csrRowPtrA,
                int *csrColIndA)>
  void init_carma_sparse_obj(CarmaContext *current_context, const long *dims,
                             T_data *M, bool load_from_host);
};

template <class T_data>
cusparseStatus_t carma_gemv(cusparseHandle_t handle, char op_A, T_data alpha,
                            CarmaSparseObj<T_data> *A, T_data *x, T_data beta,
                            T_data *y);

template <class T_data>
cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, T_data alpha,
                            CarmaSparseObj<T_data> *A, CarmaObj<T_data> *B,
                            T_data beta, CarmaObj<T_data> *C);

template <class T_data>
cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, char op_B,
                            CarmaSparseObj<T_data> *A,
                            CarmaSparseObj<T_data> *B,
                            CarmaSparseObj<T_data> *C);

template <class T_data>
cusparseStatus_t carma_csr2dense(CarmaSparseObj<T_data> *src, T_data *dest);

template <class T_data>
cusparseStatus_t carma_csr2bsr(CarmaSparseObj<T_data> *src, int block_dim,
                               CarmaSparseObj<T_data> *dest);

template <class T_data>
cusparseStatus_t carma_bsr2csr(CarmaSparseObj<T_data> *src,
                               CarmaSparseObj<T_data> *dest);

template <class T_data>
int carma_kgemv(CarmaSparseObj<T_data> *A, T_data alpha,
                const T_data *__restrict x, T_data beta, T_data *y);
#endif /* CARMA_SPARSE_OBJ_H_ */
