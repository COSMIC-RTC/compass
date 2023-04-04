// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_sparse_obj.h
//! \ingroup   libcarma
//! \class     CarmaSparseObj
//! \brief     this class provides wrappers to the generic carma sparse object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.2
//! \date      2022/01/24


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
#if CUDA_VERSION >= 11000
  cusparseSpMatDescr_t sp_descr; // New cuSPARSE API from CUDA 11
  cusparseDnMatDescr_t dn_descr; // New cuSPARSE API from CUDA 11
#endif

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

  void allocate(int nnz, int dim1, int dim2);

 private:
  void _create(int nnz_, int dim1_, int dim2_);
  void _clear();

#if CUDA_VERSION < 12000
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
# else
  void init_carma_sparse_obj(CarmaContext *current_context, const long *dims,
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
cusparseStatus_t carma_csr2bsr(CarmaSparseObj<T_data> *src, int block_dim,
                               CarmaSparseObj<T_data> *dest);

template <class T_data>
cusparseStatus_t carma_bsr2csr(CarmaSparseObj<T_data> *src,
                               CarmaSparseObj<T_data> *dest);

template <class T_data>
int carma_kgemv(CarmaSparseObj<T_data> *A, T_data alpha,
                const T_data *__restrict x, T_data beta, T_data *y);
#endif /* CARMA_SPARSE_OBJ_H_ */
