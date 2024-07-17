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

//! \file      carma_sparse_obj.hpp
//! \ingroup   libcarma
//! \class     CarmaSparseObj
//! \brief     this class provides wrappers to the generic carma sparse object
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24


#ifndef CARMA_SPARSE_OBJ_H_
#define CARMA_SPARSE_OBJ_H_
#include <cusparse_v2.h>
#include "carma_obj.hpp"

template <class T_data>
class CarmaSparseHostObj;

template <class T_data>
class CarmaSparseObj {
 public:
  int64_t dims_data[3];  ///< dimensions of the array
  int32_t nz_elem;        ///< number of elements in the array
  int32_t device;         ///< device where the CarmaObj is allocate
  CarmaContext *current_context;

  // ZERO-BASED INDEXING CSR-FORMAT
  T_data *d_data;  // nz_elem elements
  int32_t *d_rowind;   // dim1+1  elements
  int32_t *d_colind;   // nz_elem elements
  cusparseMatDescr_t descr;
#if CUDA_VERSION >= 11000
  cusparseSpMatDescr_t sp_descr; // New cuSPARSE API from CUDA 11
  cusparseDnMatDescr_t dn_descr; // New cuSPARSE API from CUDA 11
#endif

  char major_dim;
  std::string format;
  int32_t block_dim;  // blockdim for BSR format

 public:
  CarmaSparseObj(CarmaContext *current_context);
  CarmaSparseObj(CarmaObj<T_data> *M);
  CarmaSparseObj(CarmaSparseObj<T_data> *M);
  CarmaSparseObj(CarmaContext *current_context,
                   CarmaSparseHostObj<T_data> *M);
  CarmaSparseObj(CarmaContext *current_context, const int64_t *dims, T_data *M,
                   bool load_from_host);
  CarmaSparseObj(CarmaContext *current_context, const int64_t *dims,
                   T_data *values, int32_t *colind, int32_t *rowind, int32_t nz,
                   bool load_from_host);
  virtual ~CarmaSparseObj();

  void operator=(CarmaSparseObj<T_data> &M);
  void operator=(CarmaSparseHostObj<T_data> &M);

  void resize(int32_t nnz_, int32_t dim1_, int32_t dim2_);
  bool is_column_major();
  char get_major_dim() const { return major_dim; }
  void set_majorDim(char c) { major_dim = c; }

  /**< General Utilities */
  operator T_data *() { return d_data; }
  T_data *operator[](int32_t index) { return &d_data[index]; }
  T_data *get_data() { return d_data; }
  T_data *get_data(int32_t index) { return &d_data[index]; }
  const int64_t *get_dims() { return dims_data; }
  int64_t get_dims(int32_t i) { return dims_data[i]; }
  int32_t get_nonzero_elem() { return nz_elem; }
  CarmaContext *get_context() { return current_context; }

  int32_t get_device() { return device; }

  void sparse_to_host(int32_t *h_rowInd, int32_t *h_colInd, T_data *h_data);

#if CUDA_VERSION >= 11000
  void transpose();
  cudaDataType_t get_data_type() {
    cudaDataType_t data_type;
    if (std::is_same<T_data, float>::value)
      data_type = CUDA_R_32F;
    else if (std::is_same<T_data, double>::value)
      data_type = CUDA_R_64F;
    else
      std::cerr << "Unsupported data type" << std::endl;
    return data_type;
  }
#endif

  void allocate(int32_t nnz, int32_t dim1, int32_t dim2);

 private:
  void _create(int32_t nnz_, int32_t dim1_, int32_t dim2_);
  void _clear();

#if CUDA_VERSION < 12000
  template <cusparseStatus_t CUSPARSEAPI (*ptr_nnz)(
                cusparseHandle_t handle, cusparseDirection_t dirA, int32_t m, int32_t n,
                const cusparseMatDescr_t descrA, const T_data *A, int32_t lda,
                int32_t *nnzPerRowCol, int32_t *nnzTotalDevHostPtr),
            cusparseStatus_t CUSPARSEAPI (*ptr_dense2csr)(
                cusparseHandle_t handle, int32_t m, int32_t n,
                const cusparseMatDescr_t descrA, const T_data *A, int32_t lda,
                const int32_t *nnzPerRow, T_data *csrValA, int32_t *csrRowPtrA,
                int32_t *csrColIndA)>
  void init_carma_sparse_obj(CarmaContext *current_context, const int64_t *dims,
                             T_data *M, bool load_from_host);
# else
  void init_carma_sparse_obj(CarmaContext *current_context, const int64_t *dims,
                             T_data *M, bool load_from_host);

#endif
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
cusparseStatus_t carma_csr2bsr(CarmaSparseObj<T_data> *src, int32_t block_dim,
                               CarmaSparseObj<T_data> *dest);

template <class T_data>
cusparseStatus_t carma_bsr2csr(CarmaSparseObj<T_data> *src,
                               CarmaSparseObj<T_data> *dest);

template <class T_data>
int32_t carma_kgemv(CarmaSparseObj<T_data> *A, T_data alpha,
                const T_data *__restrict x, T_data beta, T_data *y);
#endif /* CARMA_SPARSE_OBJ_H_ */
