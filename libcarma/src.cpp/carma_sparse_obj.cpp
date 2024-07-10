// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_sparse_obj.cpp
//! \ingroup   libcarma
//! \class     CarmaSparseObj
//! \brief     this class provides wrappers to the generic carma sparse object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24


#include "carma_sparse_obj.hpp"
#include "carma_sparse_host_obj.hpp"
#include "carma_timer.hpp"

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaContext *current_context) {
  this->current_context = current_context;
  _create(0, 0, 0);
}

#if CUDA_VERSION < 12000
template <class T_data>
template <cusparseStatus_t CUSPARSEAPI (*ptr_nnz)(
              cusparseHandle_t handle, cusparseDirection_t dirA, int32_t m, int32_t n,
              const cusparseMatDescr_t descrA, const T_data *A, int32_t lda,
              int32_t *nnzPerRowCol, int32_t *nnzTotalDevHostPtr),
          cusparseStatus_t CUSPARSEAPI (*ptr_dense2csr)(
              cusparseHandle_t handle, int32_t m, int32_t n,
              const cusparseMatDescr_t descrA, const T_data *A, int32_t lda,
              const int32_t *nnzPerRow, T_data *csrValA, int32_t *csrRowPtrA,
              int32_t *csrColIndA)>
void CarmaSparseObj<T_data>::init_carma_sparse_obj(
    CarmaContext *current_context, const int64_t *dims, T_data *M,
    bool load_from_host) {
  _create(0, 0, 0);
  this->current_context = current_context;
  device = current_context->get_active_device();
  cusparseHandle_t handle = current_context->get_cusparse_handle();
  T_data *d_M;
  if (load_from_host) {
    cudaMalloc((void **)&d_M, dims[1] * dims[2] * sizeof(T_data));
    cudaMemcpy(d_M, M, dims[1] * dims[2] * sizeof(T_data),
               cudaMemcpyHostToDevice);
  } else {
    d_M = M;
  }
  int32_t *nnzPerRow = NULL, nnzTotalDevHostPtr = 0;
  cudaMalloc((void **)&nnzPerRow, dims[1] * sizeof(int32_t));
  // nz_elem=0;
  if (dims[1] == 1) {
    int32_t *tmp_colind = NULL;
    cudaMalloc((void **)&tmp_colind, dims[2] * sizeof(int32_t));
    find_nnz(d_M, tmp_colind, dims[2], nnzPerRow, nnzTotalDevHostPtr,
             current_context->get_device(device));
    if (!nnzTotalDevHostPtr) {
      DEBUG_TRACE("Warning : empty CarmaObj cannot be sparsed");
      return;
    }
    resize(nnzTotalDevHostPtr, dims[1], dims[2]);
    fill_sparse_vect(d_M, tmp_colind, this->d_data, this->d_colind,
                     this->d_rowind, this->nz_elem,
                     current_context->get_device(device));
    cudaFree(tmp_colind);
  } else {
    ptr_nnz(handle, CUSPARSE_DIRECTION_ROW, dims[1], dims[2], descr, d_M,
            dims[1], nnzPerRow, &nnzTotalDevHostPtr);
    resize(nnzTotalDevHostPtr, dims[1], dims[2]);
    ptr_dense2csr(handle, dims[1], dims[2], descr, d_M, dims[1], nnzPerRow,
                  this->d_data, this->d_rowind, this->d_colind);
  }
  format = "CSR";

#if CUDA_VERSION >= 11000
  carma_check_cusparse_status(cusparseCreateCsr(&sp_descr, dims[1], dims[2], nnzTotalDevHostPtr,
                                      this->d_rowind, this->d_colind, this->d_data,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, this->get_data_type()));
#endif

  cudaFree(nnzPerRow);
  if (load_from_host) {
    cudaFree(d_M);
  }
}

#else
template <class T_data>
void CarmaSparseObj<T_data>::init_carma_sparse_obj(
    CarmaContext *current_context, const int64_t *dims, T_data *M,
    bool load_from_host) {
  _create(0, 0, 0);
  this->current_context = current_context;
  device = current_context->get_active_device();
  cusparseHandle_t handle = current_context->get_cusparse_handle();
  T_data *d_M;
  if (load_from_host) {
    cudaMalloc((void **)&d_M, dims[1] * dims[2] * sizeof(T_data));
    cudaMemcpy(d_M, M, dims[1] * dims[2] * sizeof(T_data),
               cudaMemcpyHostToDevice);
  } else {
    d_M = M;
  }

  cudaMalloc((void**)&(this->d_rowind), (dims[1] + 1)*sizeof(T_data));

  carma_check_cusparse_status(cusparseCreateDnMat(&(this->dn_descr), dims[1], dims[2], dims[1], d_M,
                                        this->get_data_type(), CUSPARSE_ORDER_COL));
  carma_check_cusparse_status(cusparseCreateCsr(&(this->sp_descr), dims[1], dims[2], 0,
                                      this->d_rowind, NULL, NULL,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, this->get_data_type()));
  void* d_buffer = NULL;
  size_t bufferSize = 0;
  carma_check_cusparse_status(cusparseDenseToSparse_bufferSize(handle, this->dn_descr, this->sp_descr,
                                    CUSPARSE_DENSETOSPARSE_ALG_DEFAULT,
                                    &bufferSize));
  cudaMalloc(&d_buffer, bufferSize);
  carma_check_cusparse_status(cusparseDenseToSparse_analysis(handle, this->dn_descr, this->sp_descr,
                                  CUSPARSE_DENSETOSPARSE_ALG_DEFAULT,
                                  d_buffer));
  int64_t num_rows_tmp, num_cols_tmp, nnz;
  carma_check_cusparse_status(cusparseSpMatGetSize(this->sp_descr, &num_rows_tmp, &num_cols_tmp,
                                         &nnz));
  this->nz_elem = nnz;
  cudaMalloc((void**) &d_colind, nnz * sizeof(int32_t));
  cudaMalloc((void**) &d_data, nnz * sizeof(T_data));
  carma_check_cusparse_status(cusparseCsrSetPointers(this->sp_descr, this->d_rowind, this->d_colind,
                                           this->d_data));
  carma_check_cusparse_status(cusparseDenseToSparse_convert(handle, this->dn_descr, this->sp_descr,
                                        CUSPARSE_DENSETOSPARSE_ALG_DEFAULT,
                                        d_buffer));

  cusparseDestroyDnMat(this->dn_descr);
  cudaFree(d_buffer);
}
#endif
template <class T_data>
void CarmaSparseObj<T_data>::sparse_to_host(int32_t *h_rowInd, int32_t *h_colInd,
                                              T_data *h_data) {
  cudaMemcpy(h_rowInd, this->d_rowind, (this->dims_data[1] + 1) * sizeof(int32_t),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_colInd, this->d_colind, (this->nz_elem) * sizeof(int32_t),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_data, this->d_data, (this->nz_elem) * sizeof(T_data),
             cudaMemcpyDeviceToHost);
}

template <>
void CarmaSparseObj<double>::sparse_to_host(int32_t *h_rowInd, int32_t *h_colInd,
                                              double *h_data) {
  cudaMemcpy(h_rowInd, this->d_rowind, (this->dims_data[1] + 1) * sizeof(int32_t),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_colInd, this->d_colind, (this->nz_elem) * sizeof(int32_t),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_data, this->d_data, (this->nz_elem) * sizeof(double),
             cudaMemcpyDeviceToHost);
}

template <>
void CarmaSparseObj<float>::sparse_to_host(int32_t *h_rowInd, int32_t *h_colInd,
                                             float *h_data) {
  cudaMemcpy(h_rowInd, this->d_rowind, (this->dims_data[1] + 1) * sizeof(int32_t),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_colInd, this->d_colind, (this->nz_elem) * sizeof(int32_t),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_data, this->d_data, (this->nz_elem) * sizeof(float),
             cudaMemcpyDeviceToHost);
}

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaContext *current_context,
                                           const int64_t *dims, T_data *M,
                                           bool load_from_host) {
  _create(0, 0, 0);
}

#if CUDA_VERSION < 12000
template <>
CarmaSparseObj<float>::CarmaSparseObj(CarmaContext *current_context,
                                          const int64_t *dims, float *M,
                                          bool load_from_host) {
  init_carma_sparse_obj<cusparseSnnz, cusparseSdense2csr>(current_context, dims,
                                                          M, load_from_host);
}
template <>
CarmaSparseObj<double>::CarmaSparseObj(CarmaContext *current_context,
                                           const int64_t *dims, double *M,
                                           bool load_from_host) {
  init_carma_sparse_obj<cusparseDnnz, cusparseDdense2csr>(current_context, dims,
                                                          M, load_from_host);
}
#else
template <>
CarmaSparseObj<float>::CarmaSparseObj(CarmaContext *current_context,
                                          const int64_t *dims, float *M,
                                          bool load_from_host) {
  init_carma_sparse_obj(current_context, dims, M, load_from_host);
}
template <>
CarmaSparseObj<double>::CarmaSparseObj(CarmaContext *current_context,
                                           const int64_t *dims, double *M,
                                           bool load_from_host) {
  init_carma_sparse_obj(current_context, dims, M, load_from_host);
}

#endif
template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaContext *current_context,
                                           const int64_t *dims, T_data *values,
                                           int32_t *colind, int32_t *rowind, int32_t nz,
                                           bool load_from_host) {
  _create(nz, dims[1], dims[2]);
  this->current_context = current_context;
  device = current_context->get_active_device();
  this->format = "CSR";

  if (load_from_host) {
    cudaMemcpy(this->d_data, values, nz * sizeof(T_data),
               cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_colind, colind, nz * sizeof(int32_t),
               cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_rowind, rowind, (dims[1] + 1) * sizeof(int32_t),
               cudaMemcpyHostToDevice);
  } else {
    cudaMemcpy(this->d_data, values, nz * sizeof(T_data),
               cudaMemcpyDeviceToDevice);
    cudaMemcpy(this->d_colind, colind, nz * sizeof(int32_t),
               cudaMemcpyDeviceToDevice);
    cudaMemcpy(this->d_rowind, rowind, (dims[1] + 1) * sizeof(int32_t),
               cudaMemcpyDeviceToDevice);
  }
#if CUDA_VERSION >= 11000
  carma_check_cusparse_status(cusparseCreateCsr(&sp_descr, dims[1], dims[2], nz,
                                      this->d_rowind, this->d_colind, this->d_data,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, this->get_data_type()));
#endif
}
// template CarmaSparseObj<float>::CarmaSparseObj(
//     CarmaContext *current_context, const int64_t *dims, float *values,
//     int32_t *colind, int32_t *rowind, int32_t nz, bool load_from_host);
// template CarmaSparseObj<double>::CarmaSparseObj(
//     CarmaContext *current_context, const int64_t *dims, double *values,
//     int32_t *colind, int32_t *rowind, int32_t nz, bool load_from_host);

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaObj<T_data> *M) {
  _create(0, 0, 0);
}
#if CUDA_VERSION < 12000
template <>
CarmaSparseObj<float>::CarmaSparseObj(CarmaObj<float> *M) {
  init_carma_sparse_obj<cusparseSnnz, cusparseSdense2csr>(
      M->get_context(), M->get_dims(), M->get_data(), false);
}
template <>
CarmaSparseObj<double>::CarmaSparseObj(CarmaObj<double> *M) {
  init_carma_sparse_obj<cusparseDnnz, cusparseDdense2csr>(
      M->get_context(), M->get_dims(), M->get_data(), false);
}
#else
template <>
CarmaSparseObj<float>::CarmaSparseObj(CarmaObj<float> *M) {
  init_carma_sparse_obj(
      M->get_context(), M->get_dims(), M->get_data(), false);
}
template <>
CarmaSparseObj<double>::CarmaSparseObj(CarmaObj<double> *M) {
  init_carma_sparse_obj(
      M->get_context(), M->get_dims(), M->get_data(), false);
}
#endif
template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaSparseObj<T_data> *M) {
  this->current_context = M->current_context;

  _create(M->nz_elem, M->get_dims(1), M->get_dims(2));

  cudaMemcpy(d_data, M->d_data, nz_elem * sizeof(T_data),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_rowind, M->d_rowind, (M->get_dims(1) + 1) * sizeof(int32_t),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_colind, M->d_colind, nz_elem * sizeof(int32_t),
             cudaMemcpyDeviceToDevice);

  major_dim = M->major_dim;
  format = M->format;

  cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M->descr));
  cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M->descr));
  cusparseSetMatIndexBase(descr, cusparseGetMatIndexBase(M->descr));
  cusparseSetMatType(descr, cusparseGetMatType(M->descr));
#if CUDA_VERSION >= 11000
  carma_check_cusparse_status(cusparseCreateCsr(&sp_descr, M->get_dims(1), M->get_dims(2), M->nz_elem,
                                      this->d_rowind, this->d_colind, this->d_data,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, M->get_data_type()));
#endif
}

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaContext *current_context,
                                           CarmaSparseHostObj<T_data> *M) {
  this->device = current_context->get_active_device();
  this->current_context = current_context;

  _create(M->nz_elem, M->get_dims(1), M->get_dims(2));

  cudaMemcpy(d_data, M->h_data, nz_elem * sizeof(T_data),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_rowind, M->rowind, (M->get_dims(1) + 1) * sizeof(int32_t),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_colind, M->colind, nz_elem * sizeof(int32_t),
             cudaMemcpyHostToDevice);

  major_dim = M->get_major_dim();
  format = "CSR";
#if CUDA_VERSION >= 11000
  carma_check_cusparse_status(cusparseCreateCsr(&sp_descr, M->get_dims(1), M->get_dims(2), M->nz_elem,
                                      this->d_rowind, this->d_colind, this->d_data,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, this->get_data_type()));
#endif

}

template <class T_data>
void CarmaSparseObj<T_data>::operator=(CarmaSparseObj<T_data> &M) {
  resize(M.nz_elem, M.get_dims(1), M.get_dims(2));
  cudaMemcpy(d_data, M.d_data, nz_elem * sizeof(T_data),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_rowind, M.d_rowind, (M.get_dims(1) + 1) * sizeof(int32_t),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_colind, M.d_colind, nz_elem * sizeof(int32_t),
             cudaMemcpyDeviceToDevice);

  major_dim = M.major_dim;
  this->format = M.format;

  cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M.descr));
  cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M.descr));
  cusparseSetMatIndexBase(descr, cusparseGetMatIndexBase(M.descr));
  cusparseSetMatType(descr, cusparseGetMatType(M.descr));
#if CUDA_VERSION >= 11000
  carma_check_cusparse_status(cusparseCreateCsr(&sp_descr, M.get_dims(1), M.get_dims(2), M.nz_elem,
                                      this->d_rowind, this->d_colind, this->d_data,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, this->get_data_type()));
#endif

}

template <class T_data>
void CarmaSparseObj<T_data>::operator=(CarmaSparseHostObj<T_data> &M) {
  resize(M.nz_elem, M.get_dims(1), M.get_dims(2));
  cudaMemcpy(d_data, M.h_data, nz_elem * sizeof(T_data),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_rowind, M.rowind, (M.get_dims(1) + 1) * sizeof(int32_t),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_colind, M.colind, nz_elem * sizeof(int32_t), cudaMemcpyHostToDevice);

  major_dim = M.get_major_dim();
  this->format = "CSR";
#if CUDA_VERSION >= 11000
  carma_check_cusparse_status(cusparseCreateCsr(&sp_descr, M.get_dims(1), M.get_dims(2), M.nz_elem,
                                      this->d_rowind, this->d_colind, this->d_data,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, this->get_data_type()));
#endif
}

template <class T_data>
void CarmaSparseObj<T_data>::resize(int32_t nz_elem_, int32_t dim1_, int32_t dim2_) {
    _clear();
    _create(nz_elem_, dim1_, dim2_);
}

template <class T_data>
void CarmaSparseObj<T_data>::allocate(int32_t nz_elem, int32_t dim1, int32_t dim2) {
  if (d_data != NULL) {
    carma_check_msg(cudaFree(d_data));
    carma_check_msg(cudaFree(d_rowind));
    carma_check_msg(cudaFree(d_colind));
    carma_check_cusparse_status(cusparseDestroyMatDescr(descr));
  }
  _create(nz_elem, dim1, dim2);
#if CUDA_VERSION >= 11000
  if(sp_descr != 0)
      carma_check_cusparse_status(cusparseDestroySpMat(sp_descr));

  carma_check_cusparse_status(cusparseCreateCsr(&sp_descr, dim1, dim2, nz_elem,
                                      this->d_rowind, this->d_colind, this->d_data,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, this->get_data_type()));
#endif
}

template <class T_data>
void CarmaSparseObj<T_data>::_create(int32_t nz_elem_, int32_t dim1_, int32_t dim2_) {
  cusparseStatus_t status;
  nz_elem = nz_elem_;
  dims_data[0] = 2;
  dims_data[1] = dim1_;
  dims_data[2] = dim2_;

  if (nz_elem > 0) {
    cudaMalloc((void **)&d_data, nz_elem * sizeof(T_data));
    cudaMalloc((void **)&d_rowind, (dim1_ + 1) * sizeof(int32_t));
    cudaMalloc((void **)&d_colind, nz_elem * sizeof(int32_t));
  } else {
    d_data = NULL;
    d_rowind = d_colind = NULL;
  }

  major_dim = 'U';
  format = "CSR";

  carma_check_cusparse_status(cusparseCreateMatDescr(&descr));
#if CUDA_VERSION >= 11000
  sp_descr = 0;
#endif
  cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  // cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_NON_UNIT);
  // cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_LOWER);
  cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
}

template <class T_data>
void CarmaSparseObj<T_data>::_clear() {
  cusparseStatus_t status;
  // DEBUG_TRACE("clear %p : d_data %p d_rowind %p d_colind %p", this, d_data,
  // d_rowind, d_colind);
  if (nz_elem > 0) {
    if (d_data == NULL || d_rowind == NULL || d_colind == NULL) {
      DEBUG_TRACE("Error | CarmaSparseObj<T_data>::_clear | double clear");
      throw "Error | CarmaSparseObj<T_data>::_clear | double clear";
    }
    cudaFree(d_data);
    cudaFree(d_rowind);
    cudaFree(d_colind);
  }
  d_data = NULL;
  d_rowind = NULL;
  d_colind = NULL;
  nz_elem = 0;
  dims_data[0] = 2;
  dims_data[1] = 0;
  dims_data[2] = 0;
  major_dim = 'U';

  //  DEBUG_TRACE("clear %p : d_data %p d_rowind %p d_colind %p", this, d_data,
  //  d_rowind, d_colind);
  carma_check_cusparse_status(cusparseDestroyMatDescr(descr));
#if CUDA_VERSION >= 11000
  if(sp_descr != 0)
    carma_check_cusparse_status(cusparseDestroySpMat(sp_descr));
#endif
  //  DEBUG_TRACE("clear %p : d_data %p d_rowind %p d_colind %p", this, d_data,
  //  d_rowind, d_colind);

  descr = 0;

}
#if CUDA_VERSION >= 11000
template <class T_data>
void CarmaSparseObj<T_data>::transpose() {

  T_data *csc_val;
  int32_t *csc_rows;
  int32_t *csc_cols;

  cudaMalloc((void**)&csc_val, nz_elem * sizeof(T_data));
  cudaMalloc((void**)&csc_rows, nz_elem * sizeof(int32_t));
  cudaMalloc((void**)&csc_cols, (this->get_dims(2) + 1) * sizeof(int32_t));

  size_t buffer_size = 0;
  void *d_buffer;

  cusparseCsr2cscEx2_bufferSize(current_context->get_cusparse_handle(),
                                this->get_dims(1), this->get_dims(2), nz_elem,
                                d_data, d_rowind, d_colind,
                                csc_val, csc_cols, csc_rows,
                                this->get_data_type(), CUSPARSE_ACTION_NUMERIC,
                                CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1,
                                &buffer_size);

  cudaMalloc((void**) &d_buffer, buffer_size);

  cusparseCsr2cscEx2(current_context->get_cusparse_handle(),
                      this->get_dims(1), this->get_dims(2), nz_elem,
                                d_data, d_rowind, d_colind,
                                csc_val, csc_cols, csc_rows,
                                this->get_data_type(), CUSPARSE_ACTION_NUMERIC,
                                CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1,
                                d_buffer);

  resize(this->nz_elem, this->get_dims(2), this->get_dims(1));

  cudaMemcpy(d_data, csc_val, nz_elem * sizeof(T_data),
              cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_colind, csc_rows, nz_elem * sizeof(int32_t),
              cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_rowind, csc_cols, (this->get_dims(1) + 1) * sizeof(int32_t),
              cudaMemcpyDeviceToDevice);

  cudaFree(csc_val);
  cudaFree(csc_rows);
  cudaFree(csc_cols);
  cudaFree(d_buffer);

  if (this->major_dim == 'C')
    major_dim = 'R';
  else if (this->major_dim == 'R')
    major_dim = 'C';
  else
    major_dim = 'U';
#if CUDA_VERSION >= 11000
  if(sp_descr != 0) {
    carma_check_cusparse_status(cusparseDestroySpMat(sp_descr));
  }
  sp_descr = 0;
  carma_check_cusparse_status(cusparseCreateCsr(&sp_descr, this->get_dims(1), this->get_dims(2), this->nz_elem,
                                      this->d_rowind, this->d_colind, this->d_data,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, this->get_data_type()));
#endif

}
#endif

template <class T_data>
bool CarmaSparseObj<T_data>::is_column_major() {
  bool colMajor = true;
  /* TODO: demerdé ça
    int32_t colm1 = 0;
    int32_t rowm1 = 0;

    CarmaSparseHostObj<T_data> A_tmp;
    kp_cu2kp_smatrix(A_tmp, *this);

    for (int32_t i = 0; i < nz_elem; i++) {

      if (A_tmp.d_colind[i] == colm1) {
        colMajor = colMajor && (A_tmp.d_rowind[i] > rowm1);
        rowm1 = A_tmp.d_rowind[i];
      } else {
        rowm1 = A_tmp.d_rowind[i];
        colMajor = colMajor && (A_tmp.d_colind[i] > colm1);
        colm1 = A_tmp.d_colind[i];
      }
    }
  */
  return colMajor;
}

template <class T_data>
CarmaSparseObj<T_data>::~CarmaSparseObj<T_data>() {

  _clear();
}

template class CarmaSparseObj<float>;
template class CarmaSparseObj<double>;

// #define EXPLICITE_TEMPLATE(T_data)                                         \
//   template CarmaSparseObj<T_data>::CarmaSparseObj(                     \
//       CarmaContext *current_context);                                     \
//   template CarmaSparseObj<T_data>::CarmaSparseObj(                     \
//       CarmaSparseObj<T_data> *M);                                        \
//   template CarmaSparseObj<T_data>::CarmaSparseObj(                     \
//       CarmaContext *current_context, CarmaSparseHostObj<T_data> *M);   \
//   template CarmaSparseObj<T_data>::~CarmaSparseObj<T_data>();          \
//   template void CarmaSparseObj<T_data>::operator=(                       \
//       CarmaSparseObj<T_data> &M);                                        \
//   template void CarmaSparseObj<T_data>::operator=(                       \
//       CarmaSparseHostObj<T_data> &M);                                   \
//   template void CarmaSparseObj<T_data>::resize(int32_t nz_elem_, int32_t dim1_,  \
//                                                  int32_t dim2_);               \
//   template void CarmaSparseObj<T_data>::_create(int32_t nz_elem_, int32_t dim1_, \
//                                                   int32_t dim2_);              \
//   template void CarmaSparseObj<T_data>::_clear();                        \
//   template void CarmaSparseObj<T_data>::transpose(             \
//       CarmaSparseObj<T_data> *M);                                        \
//   template bool CarmaSparseObj<T_data>::is_column_major();

// EXPLICITE_TEMPLATE(double);
// EXPLICITE_TEMPLATE(float);

// #undef EXPLICITE_TEMPLATE
