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

//! \file      carma_cusparse.cpp
//! \ingroup   libcarma
//! \brief     this file provides the cusparse features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <carma_cusparse.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>
#include <string>

#define carma_check_cusparse_status(status) \
  carma_checkCusparseStatus_v2(status, __LINE__, __FILE__)

cusparseStatus_t carma_checkCusparseStatus_v2(cusparseStatus_t status, int line,
                                              std::string file) {
  /**< Generic CUSPARSE check status routine */
  switch (status) {
    case CUSPARSE_STATUS_SUCCESS:
      return status;
    case CUSPARSE_STATUS_NOT_INITIALIZED:
      std::cerr << "Cusparse error : The CUSPARSE library was not initialized."
                << std::endl;
      break;
    case CUSPARSE_STATUS_ALLOC_FAILED:
      std::cerr << "Cusparse error : Resource allocation failed inside the "
                   "CUSPARSE library. !!!!!"
                << std::endl;
      break;
    case CUSPARSE_STATUS_INVALID_VALUE:
      std::cerr << "Cusparse error : An unsupported value or parameter was "
                   "passed to the function."
                << std::endl;
      break;
    case CUSPARSE_STATUS_ARCH_MISMATCH:
      std::cerr << "Cusparse error : The function requires a feature absent "
                   "from the device architecture."
                << std::endl;
      break;
    case CUSPARSE_STATUS_MAPPING_ERROR:
      std::cerr << "Cusparse error : An access to GPU memory space failed, "
                   "which is usually caused by a failure to bind a texture."
                << std::endl;
      break;
    case CUSPARSE_STATUS_EXECUTION_FAILED:
      std::cerr << "Cusparse error : The function requires a feature absent "
                   "from the device architecture."
                << std::endl;
      break;
    case CUSPARSE_STATUS_INTERNAL_ERROR:
      std::cerr << "Cusparse error : An internal CUSPARSE operation failed."
                << std::endl;
      break;
    case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
      std::cerr << "Cusparse error : The matrix type is not supported by this "
                   "function."
                << std::endl;
      break;
// Define not supported status for pre-6.0 compatibility.
#if CUDA_VERSION >= 6000
    case CUSPARSE_STATUS_ZERO_PIVOT:
      std::cerr << "Cusparse error : An entry of the matrix is either "
                   "structural zero or numerical zero (singular block)"
                << std::endl;
      break;
#endif
  }
  std::cerr << "Cusparse error in " << file << "@" << line << std::endl;
  return status;
}

cusparseStatus_t carma_init_cusparse(cusparseHandle_t *cusparse_handle) {
  /**< Generic CUSPARSE init routine */
  return carma_check_cusparse_status(cusparseCreate(cusparse_handle));
}

cusparseStatus_t carma_shutdown_cusparse(cusparseHandle_t cusparse_handle) {
  /**< Generic CUSPARSE shutdown routine */
  return carma_check_cusparse_status(cusparseDestroy(cusparse_handle));
}

cusparseOperation_t carma_char2cusparse_operation(char operation) {
  switch (operation) {
    case 't':
    case 'T':
      return CUSPARSE_OPERATION_TRANSPOSE;
    case 'c':
    case 'C':
      return CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    default:
      return CUSPARSE_OPERATION_NON_TRANSPOSE;
  }
}

/*
 * _____ _____ __  __ ____  _        _  _____ _____ ____
 *|_   _| ____|  \/  |  _ \| |      / \|_   _| ____/ ___|
 *  | | |  _| | |\/| | |_) | |     / _ \ | | |  _| \___ \
 *  | | | |___| |  | |  __/| |___ / ___ \| | | |___ ___) |
 *  |_| |_____|_|  |_|_|   |_____/_/   \_\_| |_____|____/
 *
 */

template <class T_data,
          cusparseStatus_t CUSPARSEAPI (*csrmv)(
              cusparseHandle_t handle, cusparseOperation_t transA, int m, int n,
              int nnz, const T_data *alpha, const cusparseMatDescr_t descrA,
              const T_data *csrValA, const int *csrRowPtrA,
              const int *csrColIndA, const T_data *x, const T_data *beta,
              T_data *y)>
cusparseStatus_t carma_gemv(cusparseHandle_t handle, char op_A, T_data alpha,
                            CarmaSparseObj<T_data> *A, T_data *x, T_data beta,
                            T_data *y) {
  cusparseOperation_t trans = carma_char2cusparse_operation(op_A);

  return carma_check_cusparse_status(
      csrmv(handle, trans, A->get_dims(1), A->get_dims(2), A->nz_elem, &alpha,
            A->descr, A->d_data, A->d_rowind, A->d_colind, x, &beta, y));
}

template <>
cusparseStatus_t carma_gemv<double>(cusparseHandle_t handle, char op_A,
                                    double alpha, CarmaSparseObj<double> *A,
                                    double *x, double beta, double *y) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_gemv needs a CSR matrix as input");
  }
  return carma_gemv<double, cusparseDcsrmv>(handle, op_A, alpha, A, x, beta, y);
}
template <>
cusparseStatus_t carma_gemv<float>(cusparseHandle_t handle, char op_A,
                                   float alpha, CarmaSparseObj<float> *A,
                                   float *x, float beta, float *y) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_gemv needs a CSR matrix as input");
  }
  return carma_gemv<float, cusparseScsrmv>(handle, op_A, alpha, A, x, beta, y);
}

template <class T_data,
          cusparseStatus_t CUSPARSEAPI (*csrmm)(
              cusparseHandle_t handle, cusparseOperation_t transA, int m, int n,
              int k, int nnz, const T_data *alpha,
              const cusparseMatDescr_t descrA, const T_data *csrValA,
              const int *csrRowPtrA, const int *csrColIndA, const T_data *B,
              int ldb, const T_data *beta, T_data *C, int ldc)>
cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, T_data alpha,
                            CarmaSparseObj<T_data> *A, CarmaObj<T_data> *B,
                            T_data beta, CarmaObj<T_data> *C) {
  cusparseOperation_t transa = carma_char2cusparse_operation(op_A);
  cusparseStatus_t status;

  status = carma_check_cusparse_status(
      csrmm(handle, transa, A->get_dims(1), B->get_dims(2), A->get_dims(2),
            A->nz_elem, &alpha, A->descr, A->d_data, A->d_rowind, A->d_colind,
            *B, B->get_dims(1), &beta, *C, C->get_dims(1)));

  if (status != CUSPARSE_STATUS_SUCCESS) {
    std::cerr
        << "Error | carma_gemm (sparse) | Matrix-matrix multiplication failed"
        << std::endl;
    throw "Error | carma_gemm (sparse) | Matrix-matrix multiplication failed";
    // exit(EXIT_FAILURE);
  }
  return status;
}

template <>
cusparseStatus_t carma_gemm<float>(cusparseHandle_t handle, char op_A,
                                   float alpha, CarmaSparseObj<float> *A,
                                   CarmaObj<float> *B, float beta,
                                   CarmaObj<float> *C) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_gemm needs a CSR matrix as input");
  }
  return carma_gemm<float, cusparseScsrmm>(handle, op_A, alpha, A, B, beta, C);
}

template <>
cusparseStatus_t carma_gemm<double>(cusparseHandle_t handle, char op_A,
                                    double alpha, CarmaSparseObj<double> *A,
                                    CarmaObj<double> *B, double beta,
                                    CarmaObj<double> *C) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_gemm needs a CSR matrix as input");
  }
  return carma_gemm<double, cusparseDcsrmm>(handle, op_A, alpha, A, B, beta, C);
}

template <class T_data,
          cusparseStatus_t csrgemm(
              cusparseHandle_t handle, cusparseOperation_t transA,
              cusparseOperation_t transB, int m, int n, int k,
              const cusparseMatDescr_t descrA, const int nnzA,
              const T_data *csrValA, const int *csrRowPtrA,
              const int *csrColIndA, const cusparseMatDescr_t descrB,
              const int nnzB, const T_data *csrValB, const int *csrRowPtrB,
              const int *csrColIndB, const cusparseMatDescr_t descrC,
              T_data *csrValC, const int *csrRowPtrC, int *csrColIndC)>
cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, char op_B,
                            CarmaSparseObj<T_data> *A,
                            CarmaSparseObj<T_data> *B,
                            CarmaSparseObj<T_data> *C) {
  cusparseOperation_t transA = carma_char2cusparse_operation(op_A);
  cusparseOperation_t transB = carma_char2cusparse_operation(op_B);
  cusparseStatus_t status;

  const int m = (op_A == 't' ? A->get_dims(2) : A->get_dims(1));
  const int n = (op_B == 't' ? B->get_dims(1) : B->get_dims(2));
  const int k = (op_A == 't' ? A->get_dims(1) : A->get_dims(2));

  int nnzC = 0;
  // nnzTotalDevHostPtr points to host memory
  int *nnzTotalDevHostPtr = &nnzC;
  cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
  int *csrRowPtrC;
  cudaMalloc(reinterpret_cast<void **>(&csrRowPtrC), sizeof(int) * (m + 1));
  status = carma_check_cusparse_status(cusparseXcsrgemmNnz(
      handle, transA, transB, m, n, k, A->descr, A->nz_elem, A->d_rowind,
      A->d_colind, B->descr, B->nz_elem, B->d_rowind, B->d_colind, C->descr,
      csrRowPtrC, nnzTotalDevHostPtr));
  // if (status != CUSPARSE_STATUS_SUCCESS) {
  //    cerr << "Error | carma_gemm (sparse) | Matrix-matrix multiplication
  //    failed"
  //        << endl;
  //    throw "Error | carma_gemm (sparse) | Matrix-matrix multiplication
  //    failed"; exit(EXIT_FAILURE);
  //  }
  if (NULL != nnzTotalDevHostPtr) {
    nnzC = *nnzTotalDevHostPtr;
  } else {
    int baseC = 0;
    cudaMemcpy(&nnzC, csrRowPtrC + m, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&baseC, csrRowPtrC, sizeof(int), cudaMemcpyDeviceToHost);
    nnzC -= baseC;
  }
  if (nnzC > 0) {
    C->resize(nnzC, m, n);
    cudaMemcpy(C->d_rowind, csrRowPtrC, sizeof(int) * (m + 1),
               cudaMemcpyDeviceToDevice);
    status = carma_check_cusparse_status(csrgemm(
        handle, transA, transB, m, n, k, A->descr, A->nz_elem, A->d_data,
        A->d_rowind, A->d_colind, B->descr, B->nz_elem, B->d_data, B->d_rowind,
        B->d_colind, C->descr, C->d_data, C->d_rowind, C->d_colind));
    if (status != CUSPARSE_STATUS_SUCCESS) {
      std::cerr
          << "Error | carma_gemm (sparse) | Matrix-matrix multiplication failed"
          << std::endl;
      throw "Error | carma_gemm (sparse) | Matrix-matrix multiplication failed";
      // exit(EXIT_FAILURE);
    }
  }
  cudaFree(csrRowPtrC);
  return status;
}

template <>
cusparseStatus_t carma_gemm<float>(cusparseHandle_t handle, char op_A,
                                   char op_B, CarmaSparseObj<float> *A,
                                   CarmaSparseObj<float> *B,
                                   CarmaSparseObj<float> *C) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_gemm needs a CSR matrix as input");
  }
  return carma_gemm<float, cusparseScsrgemm>(handle, op_A, op_B, A, B, C);
}

template <>
cusparseStatus_t carma_gemm<double>(cusparseHandle_t handle, char op_A,
                                    char op_B, CarmaSparseObj<double> *A,
                                    CarmaSparseObj<double> *B,
                                    CarmaSparseObj<double> *C) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_gemm needs a CSR matrix as input");
  }
  return carma_gemm<double, cusparseDcsrgemm>(handle, op_A, op_B, A, B, C);
}

template <class T_data,
          cusparseStatus_t csr2dense(
              cusparseHandle_t handle, int m, int n,
              const cusparseMatDescr_t descrA, const T_data *csrValA,
              const int *csrRowPtrA, const int *csrColIndA, T_data *B, int ldb)>
cusparseStatus_t carma_csr2dense_gen(CarmaSparseObj<T_data> *A, T_data *B) {
  cusparseStatus_t status;
  cusparseHandle_t handle = A->current_context->get_cusparse_handle();

  status = csr2dense(handle, A->dims_data[1], A->dims_data[2], A->descr,
                     A->d_data, A->d_rowind, A->d_colind, B, A->dims_data[1]);
  if (status != CUSPARSE_STATUS_SUCCESS) {
    std::cerr << "Error | carma_csr2dense (sparse) | csr2dense failed"
              << std::endl;
    throw "Error | carma_csr2dense (sparse) | csr2dense failed";
  }
  return status;
}

template <>
cusparseStatus_t carma_csr2dense(CarmaSparseObj<float> *A, float *B) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_csr2dense needs a CSR matrix as input");
  }
  return carma_csr2dense_gen<float, cusparseScsr2dense>(A, B);
}

template <>
cusparseStatus_t carma_csr2dense(CarmaSparseObj<double> *A, double *B) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_csr2dense needs a CSR matrix as input");
  }
  return carma_csr2dense_gen<double, cusparseDcsr2dense>(A, B);
}

template <class T_data,
          cusparseStatus_t (*csr2bsr)(
              cusparseHandle_t handle, cusparseDirection_t dir, int m, int n,
              const cusparseMatDescr_t descrA, const T_data *csrValA,
              const int *csrRowPtrA, const int *csrColIndA, int block_dim,
              const cusparseMatDescr_t descrC, T_data *bsrValC, int *bsrRowPtrC,
              int *bsrColIndC)>
cusparseStatus_t carma_csr2bsr_gen(CarmaSparseObj<T_data> *A, int block_dim,
                                   CarmaSparseObj<T_data> *B) {
  // Given CSR format (csrRowPtrA, csrcolIndA, csrValA) and
  // blocks of BSR format are stored in column-major order.
  if (B->nz_elem > 0) {
    cudaFree(B->d_rowind);
    cudaFree(B->d_colind);
    cudaFree(B->d_data);
  }

  cusparseDirection_t dir = CUSPARSE_DIRECTION_COLUMN;
  cusparseHandle_t handle = A->current_context->get_cusparse_handle();
  int m = A->dims_data[1];
  int n = A->dims_data[2];
  int mb = (m + block_dim - 1) / block_dim;
  cudaMalloc(reinterpret_cast<void **>(&B->d_rowind), sizeof(int) * (mb + 1));

  int nnzb;
  // nnzTotalDevHostPtr points to host memory
  cusparseXcsr2bsrNnz(handle, dir, m, n, A->descr, A->d_rowind, A->d_colind,
                      block_dim, B->descr, B->d_rowind, &nnzb);
  B->dims_data[0] = 2;
  B->dims_data[1] = m;
  B->dims_data[2] = n;
  B->nz_elem = nnzb;
  B->format = "BSR";
  B->block_dim = block_dim;
  cudaMalloc(reinterpret_cast<void **>(&B->d_colind), sizeof(int) * B->nz_elem);
  cudaMalloc(reinterpret_cast<void **>(&B->d_data),
             sizeof(T_data) * (B->block_dim * B->block_dim) * B->nz_elem);
  return csr2bsr(handle, dir, m, n, A->descr, A->d_data, A->d_rowind,
                 A->d_colind, B->block_dim, B->descr, B->d_data, B->d_rowind,
                 B->d_colind);
}

template <>
cusparseStatus_t carma_csr2bsr(CarmaSparseObj<float> *src, int block_dim,
                               CarmaSparseObj<float> *dest) {
  if (src->format != "CSR") {
    DEBUG_TRACE("carma_csr2bsr needs a CSR matrix as input");
  }
  return carma_csr2bsr_gen<float, cusparseScsr2bsr>(src, block_dim, dest);
}

template <>
cusparseStatus_t carma_csr2bsr(CarmaSparseObj<double> *src, int block_dim,
                               CarmaSparseObj<double> *dest) {
  if (src->format != "CSR") {
    DEBUG_TRACE("carma_csr2bsr needs a CSR matrix as input");
  }
  return carma_csr2bsr_gen<double, cusparseDcsr2bsr>(src, block_dim, dest);
}

template <class T_data,
          cusparseStatus_t (*bsr2csr)(
              cusparseHandle_t handle, cusparseDirection_t dir, int mb, int nb,
              const cusparseMatDescr_t descrA, const T_data *bsrValA,
              const int *bsrRowPtrA, const int *bsrColIndA, int block_dim,
              const cusparseMatDescr_t descrC, T_data *csrValC, int *csrRowPtrC,
              int *csrColIndC)>
cusparseStatus_t carma_bsr2csr_gen(CarmaSparseObj<T_data> *A,
                                   CarmaSparseObj<T_data> *B) {
  // Given BSR format (bsrRowPtrA, bsrcolIndA, bsrValA) and
  // blocks of BSR format are stored in column-major order.
  if (B->nz_elem > 0) {
    cudaFree(B->d_rowind);
    cudaFree(B->d_colind);
    cudaFree(B->d_data);
  }

  cusparseDirection_t dir = CUSPARSE_DIRECTION_COLUMN;
  cusparseHandle_t handle = A->current_context->get_cusparse_handle();
  int m = A->dims_data[1];
  int n = A->dims_data[2];
  int mb = (m + A->block_dim - 1) / A->block_dim;
  int nb = (n + A->block_dim - 1) / A->block_dim;
  int nnzb = A->nz_elem;                       // number of blocks
  int nnz = nnzb * A->block_dim * A->block_dim;  // number of elements

  B->dims_data[0] = 2;
  B->dims_data[1] = m;
  B->dims_data[2] = n;
  B->nz_elem = nnz;
  B->format = "CSR";
  cudaMalloc(reinterpret_cast<void **>(&B->d_rowind), sizeof(int) * (m + 1));
  cudaMalloc(reinterpret_cast<void **>(&B->d_colind), sizeof(int) * nnz);
  cudaMalloc(reinterpret_cast<void **>(&B->d_data), sizeof(T_data) * nnz);
  return bsr2csr(handle, dir, mb, nb, A->descr, A->d_data, A->d_rowind,
                 A->d_colind, A->block_dim, B->descr, B->d_data, B->d_rowind,
                 B->d_colind);
}

template <>
cusparseStatus_t carma_bsr2csr(CarmaSparseObj<float> *src,
                               CarmaSparseObj<float> *dest) {
  if (src->format != "BSR") {
    DEBUG_TRACE("carma_bsr2csr needs a BSR matrix as input");
  }
  return carma_bsr2csr_gen<float, cusparseSbsr2csr>(src, dest);
}

template <>
cusparseStatus_t carma_bsr2csr(CarmaSparseObj<double> *src,
                               CarmaSparseObj<double> *dest) {
  if (src->format != "BSR") {
    DEBUG_TRACE("carma_bsr2csr needs a BSR matrix as input");
  }
  return carma_bsr2csr_gen<double, cusparseDbsr2csr>(src, dest);
}
