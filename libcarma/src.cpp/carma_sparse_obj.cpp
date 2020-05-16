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

//! \file      carma_sparse_obj.cpp
//! \ingroup   libcarma
//! \class     CarmaSparseObj
//! \brief     this class provides wrappers to the generic carma sparse object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#include "carma_sparse_obj.h"
#include "carma_sparse_host_obj.h"
#include "carma_timer.h"

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaContext *current_context) {
  this->current_context = current_context;
  _create(0, 0, 0);
}

template <class T_data>
template <cusparseStatus_t CUSPARSEAPI (*ptr_nnz)(
              cusparseHandle_t handle, cusparseDirection_t dirA, int m, int n,
              const cusparseMatDescr_t descrA, const T_data *A, int lda,
              int *nnzPerRowCol, int *nnzTotalDevHostPtr),
          cusparseStatus_t CUSPARSEAPI (*ptr_dense2csr)(
              cusparseHandle_t handle, int m, int n,
              const cusparseMatDescr_t descrA, const T_data *A, int lda,
              const int *nnzPerRow, T_data *csrValA, int *csrRowPtrA,
              int *csrColIndA)>
void CarmaSparseObj<T_data>::init_carma_sparse_obj(
    CarmaContext *current_context, const long *dims, T_data *M,
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
  int *nnzPerRow = NULL, nnzTotalDevHostPtr = 0;
  cudaMalloc((void **)&nnzPerRow, dims[1] * sizeof(int));
  // nz_elem=0;
  if (dims[1] == 1) {
    int *tmp_colind = NULL;
    cudaMalloc((void **)&tmp_colind, dims[2] * sizeof(int));
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
#ifndef USE_MAGMA_SPARSE
  d_sparse_mat = s_sparse_mat = 0L;
#endif

  cudaFree(nnzPerRow);
  if (load_from_host) {
    cudaFree(d_M);
  }
}

template <class T_data>
void CarmaSparseObj<T_data>::sparse_to_host(int *h_rowInd, int *h_colInd,
                                              T_data *h_data) {
  cudaMemcpy(h_rowInd, this->d_rowind, (this->dims_data[1] + 1) * sizeof(int),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_colInd, this->d_colind, (this->nz_elem) * sizeof(int),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_data, this->d_data, (this->nz_elem) * sizeof(T_data),
             cudaMemcpyDeviceToHost);
}

template <>
void CarmaSparseObj<double>::sparse_to_host(int *h_rowInd, int *h_colInd,
                                              double *h_data) {
  cudaMemcpy(h_rowInd, this->d_rowind, (this->dims_data[1] + 1) * sizeof(int),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_colInd, this->d_colind, (this->nz_elem) * sizeof(int),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_data, this->d_data, (this->nz_elem) * sizeof(double),
             cudaMemcpyDeviceToHost);
}

template <>
void CarmaSparseObj<float>::sparse_to_host(int *h_rowInd, int *h_colInd,
                                             float *h_data) {
  cudaMemcpy(h_rowInd, this->d_rowind, (this->dims_data[1] + 1) * sizeof(int),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_colInd, this->d_colind, (this->nz_elem) * sizeof(int),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(h_data, this->d_data, (this->nz_elem) * sizeof(float),
             cudaMemcpyDeviceToHost);
}

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaContext *current_context,
                                           const long *dims, T_data *M,
                                           bool load_from_host) {
  _create(0, 0, 0);
}
template <>
CarmaSparseObj<float>::CarmaSparseObj(CarmaContext *current_context,
                                          const long *dims, float *M,
                                          bool load_from_host) {
  init_carma_sparse_obj<cusparseSnnz, cusparseSdense2csr>(current_context, dims,
                                                          M, load_from_host);
}
template <>
CarmaSparseObj<double>::CarmaSparseObj(CarmaContext *current_context,
                                           const long *dims, double *M,
                                           bool load_from_host) {
  init_carma_sparse_obj<cusparseDnnz, cusparseDdense2csr>(current_context, dims,
                                                          M, load_from_host);
}

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaContext *current_context,
                                           const long *dims, T_data *values,
                                           int *colind, int *rowind, int nz,
                                           bool load_from_host) {
  _create(nz, dims[1], dims[2]);
  this->current_context = current_context;
  device = current_context->get_active_device();
  this->format = "CSR";

  if (load_from_host) {
    cudaMemcpy(this->d_data, values, nz * sizeof(T_data),
               cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_colind, colind, nz * sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_rowind, rowind, (dims[1] + 1) * sizeof(int),
               cudaMemcpyHostToDevice);
  } else {
    cudaMemcpy(this->d_data, values, nz * sizeof(T_data),
               cudaMemcpyDeviceToDevice);
    cudaMemcpy(this->d_colind, colind, nz * sizeof(int),
               cudaMemcpyDeviceToDevice);
    cudaMemcpy(this->d_rowind, rowind, (dims[1] + 1) * sizeof(int),
               cudaMemcpyDeviceToDevice);
  }
}
template CarmaSparseObj<float>::CarmaSparseObj(
    CarmaContext *current_context, const long *dims, float *values,
    int *colind, int *rowind, int nz, bool load_from_host);
template CarmaSparseObj<double>::CarmaSparseObj(
    CarmaContext *current_context, const long *dims, double *values,
    int *colind, int *rowind, int nz, bool load_from_host);

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaObj<T_data> *M) {
  _create(0, 0, 0);
}
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

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaSparseObj<T_data> *M) {
  this->current_context = M->current_context;

  _create(M->nz_elem, M->get_dims(1), M->get_dims(2));

  cudaMemcpy(d_data, M->d_data, nz_elem * sizeof(T_data),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_rowind, M->d_rowind, (M->get_dims(1) + 1) * sizeof(int),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_colind, M->d_colind, nz_elem * sizeof(int),
             cudaMemcpyDeviceToDevice);

  major_dim = M->major_dim;
  format = M->format;

  d_sparse_mat = M->d_sparse_mat;
  s_sparse_mat = M->s_sparse_mat;

  cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M->descr));
  cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M->descr));
  cusparseSetMatIndexBase(descr, cusparseGetMatIndexBase(M->descr));
  cusparseSetMatType(descr, cusparseGetMatType(M->descr));
}

template <class T_data>
CarmaSparseObj<T_data>::CarmaSparseObj(CarmaContext *current_context,
                                           CarmaSparseHostObj<T_data> *M) {
  this->device = current_context->get_active_device();
  this->current_context = current_context;

  _create(M->nz_elem, M->get_dims(1), M->get_dims(2));

  cudaMemcpy(d_data, M->h_data, nz_elem * sizeof(T_data),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_rowind, M->rowind, (M->get_dims(1) + 1) * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_colind, M->colind, nz_elem * sizeof(int),
             cudaMemcpyHostToDevice);

  major_dim = M->get_major_dim();
  format = "CSR";
}

template <class T_data>
void CarmaSparseObj<T_data>::operator=(CarmaSparseObj<T_data> &M) {
  resize(M.nz_elem, M.get_dims(1), M.get_dims(2));
  cudaMemcpy(d_data, M.d_data, nz_elem * sizeof(T_data),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_rowind, M.d_rowind, (M.get_dims(1) + 1) * sizeof(int),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_colind, M.d_colind, nz_elem * sizeof(int),
             cudaMemcpyDeviceToDevice);

  major_dim = M.major_dim;
  this->format = M.format;

  d_sparse_mat = M.d_sparse_mat;
  s_sparse_mat = M.s_sparse_mat;

  cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M.descr));
  cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M.descr));
  cusparseSetMatIndexBase(descr, cusparseGetMatIndexBase(M.descr));
  cusparseSetMatType(descr, cusparseGetMatType(M.descr));
}

template <class T_data>
void CarmaSparseObj<T_data>::operator=(CarmaSparseHostObj<T_data> &M) {
  resize(M.nz_elem, M.get_dims(1), M.get_dims(2));
  cudaMemcpy(d_data, M.h_data, nz_elem * sizeof(T_data),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_rowind, M.rowind, (M.get_dims(1) + 1) * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_colind, M.colind, nz_elem * sizeof(int), cudaMemcpyHostToDevice);

  major_dim = M.get_major_dim();
  this->format = "CSR";
}

template <class T_data>
void CarmaSparseObj<T_data>::resize(int nz_elem_, int dim1_, int dim2_) {
  if (nz_elem != nz_elem_) {
    _clear();
    _create(nz_elem_, dim1_, dim2_);
  } else {
    dims_data[0] = 2;
    dims_data[1] = dim1_;
    dims_data[2] = dim2_;
  }
}

template <class T_data>
void CarmaSparseObj<T_data>::_create(int nz_elem_, int dim1_, int dim2_) {
  cusparseStatus_t status;
  nz_elem = nz_elem_;
  dims_data[0] = 2;
  dims_data[1] = dim1_;
  dims_data[2] = dim2_;

  if (nz_elem > 0) {
    cudaMalloc((void **)&d_data, nz_elem * sizeof(T_data));
    cudaMalloc((void **)&d_rowind, (dim1_ + 1) * sizeof(int));
    cudaMalloc((void **)&d_colind, nz_elem * sizeof(int));
  } else {
    d_data = NULL;
    d_rowind = d_colind = NULL;
  }

  major_dim = 'U';
  format = "CSR";

  status = cusparseCreateMatDescr(&descr);
  if (status != CUSPARSE_STATUS_SUCCESS) {
    DEBUG_TRACE(
        "Error | CarmaSparseObj<T_data>::_create | Matrix descriptor "
        "initialization failed");
    throw "Error | CarmaSparseObj<T_data>::_create | Matrix descriptor initialization failed";
    // exit(EXIT_FAILURE);
  }

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
  status = cusparseDestroyMatDescr(descr);
  //  DEBUG_TRACE("clear %p : d_data %p d_rowind %p d_colind %p", this, d_data,
  //  d_rowind, d_colind);

  descr = 0;

  if (status != CUSPARSE_STATUS_SUCCESS) {
    DEBUG_TRACE(
        "Error | CarmaSparseObj<T_data>::_clear | Matrix descriptor "
        "destruction failed");
    throw "Error | CarmaSparseObj<T_data>::_clear | Matrix descriptor destruction failed";
    // exit(EXIT_FAILURE);
  }
}

template <class T_data>
void CarmaSparseObj<T_data>::init_from_transpose(
    CarmaSparseObj<T_data> *M) {
  resize(M->nz_elem, M->get_dims(1), M->get_dims(2));

  cudaMemcpy(d_data, M->d_data, nz_elem * sizeof(T_data),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_rowind, M->d_rowind, (M->get_dims(1) + 1) * sizeof(int),
             cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_colind, M->d_colind, nz_elem * sizeof(int),
             cudaMemcpyDeviceToDevice);

  d_sparse_mat = M->d_sparse_mat;
  s_sparse_mat = M->s_sparse_mat;

  if (M->major_dim == 'C')
    major_dim = 'R';
  else if (M->major_dim == 'R')
    major_dim = 'C';
  else
    major_dim = 'U';
}

template <class T_data>
bool CarmaSparseObj<T_data>::is_column_major() {
  bool colMajor = true;
  /* TODO: demerdé ça
    int colm1 = 0;
    int rowm1 = 0;

    CarmaSparseHostObj<T_data> A_tmp;
    kp_cu2kp_smatrix(A_tmp, *this);

    for (int i = 0; i < nz_elem; i++) {

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
#ifdef USE_MAGMA_SPARSE
  carma_magma_sparse_free<T_data>(this);
#endif

  _clear();
}

#define EXPLICITE_TEMPLATE(T_data)                                         \
  template CarmaSparseObj<T_data>::CarmaSparseObj(                     \
      CarmaContext *current_context);                                     \
  template CarmaSparseObj<T_data>::CarmaSparseObj(                     \
      CarmaSparseObj<T_data> *M);                                        \
  template CarmaSparseObj<T_data>::CarmaSparseObj(                     \
      CarmaContext *current_context, CarmaSparseHostObj<T_data> *M);   \
  template CarmaSparseObj<T_data>::~CarmaSparseObj<T_data>();          \
  template void CarmaSparseObj<T_data>::operator=(                       \
      CarmaSparseObj<T_data> &M);                                        \
  template void CarmaSparseObj<T_data>::operator=(                       \
      CarmaSparseHostObj<T_data> &M);                                   \
  template void CarmaSparseObj<T_data>::resize(int nz_elem_, int dim1_,  \
                                                 int dim2_);               \
  template void CarmaSparseObj<T_data>::_create(int nz_elem_, int dim1_, \
                                                  int dim2_);              \
  template void CarmaSparseObj<T_data>::_clear();                        \
  template void CarmaSparseObj<T_data>::init_from_transpose(             \
      CarmaSparseObj<T_data> *M);                                        \
  template bool CarmaSparseObj<T_data>::is_column_major();

EXPLICITE_TEMPLATE(double);
EXPLICITE_TEMPLATE(float);

#undef EXPLICITE_TEMPLATE
