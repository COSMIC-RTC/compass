// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_cusparse.cpp
//! \ingroup   libcarma
//! \brief     this file provides the cusparse features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <carma_cusparse.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>
#include <string>

cusparseStatus_t carma_check_cusparse_status_v2(cusparseStatus_t status, int line,
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
    case CUSPARSE_STATUS_NOT_SUPPORTED:
      std::cerr << "The operation or data type combination is currently not supported by the function" << std::endl;
// Define not supported status for pre-6.0 compatibility.
#if CUDA_VERSION >= 6000
    case CUSPARSE_STATUS_ZERO_PIVOT:
      std::cerr << "Cusparse error : An entry of the matrix is either "
                   "structural zero or numerical zero (singular block)"
                << std::endl;
      break;
#endif
#if CUDA_VERSION >= 11000
    case CUSPARSE_STATUS_INSUFFICIENT_RESOURCES:
      std::cerr << "The resources for the computation, such as GPU global or shared memory, are not sufficient to complete the operation" << std::endl;
#endif
  }
  std::cerr << "Cusparse error in " << file << "@" << line << " : " << cusparseGetErrorString(status) << std::endl;
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

namespace detail {
  template<typename T> struct Sparse;

  // Specialization of Sparse with float routine.
  template<> struct Sparse<float>
  {
#if CUDA_VERSION < 11000
    constexpr static auto mv()    {  return cusparseScsrmv; }
    constexpr static auto mm()    { return cusparseScsrmm; }
    constexpr static auto gemm()  { return cusparseScsrgemm; }
    constexpr static auto dense() { return cusparseScsr2dense; }
#endif
    constexpr static auto bsr()   { return cusparseScsr2bsr; }
    constexpr static auto csr()   { return cusparseSbsr2csr; }
  };

  // Specialization of Sparse with double routine.
  template<> struct Sparse<double>
  {
#if CUDA_VERSION < 11000
    constexpr static auto mv()    { return cusparseDcsrmv; }
    constexpr static auto mm()    { return cusparseDcsrmm; }
    constexpr static auto gemm()  { return cusparseDcsrgemm; }
    constexpr static auto dense() { return cusparseDcsr2dense; }
#endif
    constexpr static auto bsr()   { return cusparseDcsr2bsr; }
    constexpr static auto csr()   { return cusparseDbsr2csr; }
  };
}

// Type dependant function binded to corresponding routine.
#if CUDA_VERSION < 11000
template<typename T> constexpr auto sparse_mv =    detail::Sparse<T>::mv();
template<typename T> constexpr auto sparse_mm =    detail::Sparse<T>::mm();
template<typename T> constexpr auto sparse_gemm =  detail::Sparse<T>::gemm();
template<typename T> constexpr auto sparse_dense = detail::Sparse<T>::dense();
#endif
template<typename T> constexpr auto sparse_bsr =   detail::Sparse<T>::bsr();
template<typename T> constexpr auto sparse_csr =   detail::Sparse<T>::csr();

#if CUDA_VERSION >= 11000
  template <class T_data>
  cusparseStatus_t carma_gemv(cusparseHandle_t handle, char op_A, T_data alpha,
                            CarmaSparseObj<T_data> *A, T_data *x, T_data beta,
                            T_data *y) {
    cusparseStatus_t status;
    cusparseDnVecDescr_t vec_x, vec_y;
    cusparseSpMatDescr_t matA;
    CarmaSparseObj<T_data> *At = nullptr;
    int64_t A_num_rows, A_num_cols, A_nnz;
    void*  d_buffer    = NULL;
    size_t buffer_size = 0;
    cusparseOperation_t trans_A = CUSPARSE_OPERATION_NON_TRANSPOSE;
    if(op_A == 't') { // Have to transpose manually
      At = new CarmaSparseObj<T_data>(A);
      At->transpose();
      matA = At->sp_descr;
    }
    else {
      matA = A->sp_descr;
    }
    carma_check_cusparse_status(cusparseSpMatGetSize(matA, &A_num_rows, &A_num_cols, &A_nnz));
    carma_check_cusparse_status(cusparseCreateDnVec(&vec_x, A_num_cols, x, A->get_data_type()));
    carma_check_cusparse_status(cusparseCreateDnVec(&vec_y, A_num_rows, y, A->get_data_type()));
    carma_check_cusparse_status(cusparseSpMV_bufferSize(
                                 handle, trans_A,
                                 &alpha, matA, vec_x, &beta, vec_y, A->get_data_type(),
                                 CUSPARSE_SPMV_ALG_DEFAULT, &buffer_size));
    cudaMalloc(&d_buffer, buffer_size);
    status = carma_check_cusparse_status(cusparseSpMV(handle, trans_A,
                                 &alpha, matA, vec_x, &beta, vec_y, A->get_data_type(),
                                 CUSPARSE_SPMV_ALG_DEFAULT, d_buffer));
    carma_check_cusparse_status(cusparseDestroyDnVec(vec_x));
    carma_check_cusparse_status(cusparseDestroyDnVec(vec_y));
    carma_check_msg(cudaFree(d_buffer));
    if(At != nullptr)
      delete At;

    return status;
  }

  template <class T_data>
  cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, T_data alpha,
                              CarmaSparseObj<T_data> *A, CarmaObj<T_data> *B,
                              T_data beta, CarmaObj<T_data> *C) {
    cusparseStatus_t status;
    cusparseDnMatDescr_t mat_b, mat_c;
    cusparseOperation_t trans_A = carma_char2cusparse_operation(op_A);

    void*  d_buffer    = NULL;
    size_t buffer_size = 0;
    carma_check_cusparse_status(cusparseCreateDnMat(&mat_b, B->get_dims(1), B->get_dims(2), B->get_dims(1),
                                                    B->get_data(), A->get_data_type(), CUSPARSE_ORDER_COL));
    carma_check_cusparse_status(cusparseCreateDnMat(&mat_c, C->get_dims(1), C->get_dims(2), C->get_dims(1),
                                                    C->get_data(), A->get_data_type(), CUSPARSE_ORDER_COL));
    carma_check_cusparse_status(cusparseSpMM_bufferSize(
                                 handle, trans_A, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, A->sp_descr, mat_b, &beta, mat_c, A->get_data_type(),
                                 CUSPARSE_SPMM_ALG_DEFAULT, &buffer_size));
    carma_check_msg(cudaMalloc(&d_buffer, buffersize));
    status = carma_check_cusparse_status(cusparseSpMM(handle, trans_A, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, A->sp_descr, mat_b, &beta, mat_c, A->get_data_type(),
                                 CUSPARSE_SPMM_ALG_DEFAULT, d_buffer));
    carma_check_cusparse_status(cusparseDestroyDnMat(mat_b));
    carma_check_cusparse_status(cusparseDestroyDnMat(mat_c));
    carma_check_msg(cudaFree(d_buffer));

    return status;
  }
  template <class T_data>
  cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, char op_B,
                            CarmaSparseObj<T_data> *A,
                            CarmaSparseObj<T_data> *B,
                            CarmaSparseObj<T_data> *C) {
    cusparseStatus_t status;
    cusparseSpMatDescr_t matA, matB, matC;

    int64_t A_num_rows, A_num_cols, A_nnz;
    int64_t B_num_rows, B_num_cols, B_nnz;
    cusparseOperation_t trans_A = CUSPARSE_OPERATION_NON_TRANSPOSE; // Only one supported
    cusparseOperation_t trans_B = CUSPARSE_OPERATION_NON_TRANSPOSE; // Only one supported
    CarmaSparseObj<T_data> *At = nullptr;
    CarmaSparseObj<T_data> *Bt = nullptr;
    if(op_A == 't') { // Have to transpose manually
      At = new CarmaSparseObj<T_data>(A);
      At->transpose();
      matA = At->sp_descr;
    }
    else {
      matA = A->sp_descr;
    }
    if(op_B == 't') { // Have to transpose manually
      Bt = new CarmaSparseObj<T_data>(B);
      Bt->transpose();
      matB = Bt->sp_descr;

    }
    else {
      matB = B->sp_descr;
    }
    cusparseSpMatGetSize(matA, &A_num_rows, &A_num_cols, &A_nnz);
    cusparseSpMatGetSize(matB, &B_num_rows, &B_num_cols, &B_nnz);

    void* d_buffer1 = NULL, *d_buffer2 = NULL;
    size_t buffer_size1 = 0, buffer_size2 = 0;
    T_data alpha = T_data(1);
    T_data beta = T_data(0);
    // Init C descriptor --> C must have been created with CarmaSparseObj(context)
    carma_check_cusparse_status(cusparseCreateCsr(&(C->sp_descr), A_num_rows, B_num_cols, 0,
                                      NULL, NULL, NULL,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, A->get_data_type()));
    // SpGEMM Computation
    cusparseSpGEMMDescr_t spgemm_desc;
    carma_check_cusparse_status(cusparseSpGEMM_createDescr(&spgemm_desc));
    // ask bufferSize1 bytes for external memory
    carma_check_cusparse_status(cusparseSpGEMM_workEstimation(handle, trans_A, trans_B,
                                  &alpha, matA, matB, &beta, C->sp_descr,
                                  A->get_data_type(), CUSPARSE_SPGEMM_DEFAULT,
                                  spgemm_desc, &buffer_size1, NULL));
    cudaMalloc((void**) &d_buffer1, buffer_size1);
    // inspect the matrices A and B to understand the memory requiremnent for
    // the next step
    carma_check_cusparse_status(cusparseSpGEMM_workEstimation(handle, trans_A, trans_B,
                                  &alpha, matA, matB, &beta, C->sp_descr,
                                  A->get_data_type(), CUSPARSE_SPGEMM_DEFAULT,
                                  spgemm_desc, &buffer_size1, d_buffer1));
    // ask bufferSize2 bytes for external memory
    carma_check_cusparse_status(cusparseSpGEMM_compute(handle, trans_A, trans_B,
                                  &alpha, matA, matB, &beta, C->sp_descr,
                                  A->get_data_type(), CUSPARSE_SPGEMM_DEFAULT,
                                  spgemm_desc, &buffer_size2, NULL));

    cudaMalloc((void**) &d_buffer2, buffer_size2);

    // compute the intermediate product of A * B
    carma_check_cusparse_status(cusparseSpGEMM_compute(handle, trans_A, trans_B,
                            &alpha, matA, matB, &beta, C->sp_descr,
                            A->get_data_type(), CUSPARSE_SPGEMM_DEFAULT,
                            spgemm_desc, &buffer_size2, d_buffer2));
    // get matrix C non-zero entries C_num_nnz1
    int64_t C_num_rows1, C_num_cols1, C_num_nnz1;
    cusparseSpMatGetSize(C->sp_descr, &C_num_rows1, &C_num_cols1, &C_num_nnz1);
    // allocate matrix C and update descriptor
    C->allocate(C_num_nnz1, A_num_rows, B_num_cols);

    // copy the final products to the matrix C
    cusparseSpGEMM_copy(handle, trans_A, trans_B,
                        &alpha, matA, matB, &beta, C->sp_descr,
                        A->get_data_type(), CUSPARSE_SPGEMM_DEFAULT, spgemm_desc);
    // destroy matrix/vector descriptors
    status = carma_check_cusparse_status(cusparseSpGEMM_destroyDescr(spgemm_desc));
    cudaFree(d_buffer1);
    cudaFree(d_buffer2);
    if(At != nullptr)
      delete At;
    if(Bt != nullptr)
      delete Bt;

    return status;
  }

#else
  template <class T_data>
  cusparseStatus_t carma_gemv(cusparseHandle_t handle, char op_A, T_data alpha,
                              CarmaSparseObj<T_data> *A, T_data *x, T_data beta,
                              T_data *y) {
    if (A->format != "CSR") {
      DEBUG_TRACE("carma_gemv needs a CSR matrix as input");
    }

    cusparseOperation_t trans = carma_char2cusparse_operation(op_A);

    return carma_check_cusparse_status(
        sparse_mv<T_data>(handle, trans, A->get_dims(1), A->get_dims(2), A->nz_elem, &alpha,
              A->descr, A->d_data, A->d_rowind, A->d_colind, x, &beta, y));
  }

  template <class T_data>
  cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, T_data alpha,
                              CarmaSparseObj<T_data> *A, CarmaObj<T_data> *B,
                              T_data beta, CarmaObj<T_data> *C) {
    if (A->format != "CSR") {
      DEBUG_TRACE("carma_gemm needs a CSR matrix as input");
    }

    cusparseOperation_t transa = carma_char2cusparse_operation(op_A);
    cusparseStatus_t status;

    status = carma_check_cusparse_status(
        sparse_mm<T_data>(handle, transa, A->get_dims(1), B->get_dims(2), A->get_dims(2),
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


  template <class T_data>
  cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, char op_B,
                              CarmaSparseObj<T_data> *A,
                              CarmaSparseObj<T_data> *B,
                              CarmaSparseObj<T_data> *C) {
    if (A->format != "CSR") {
      DEBUG_TRACE("carma_gemm needs a CSR matrix as input");
    }

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
      status = carma_check_cusparse_status(sparse_gemm<T_data>(
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
#endif

template <class T_data>
cusparseStatus_t carma_csr2dense(CarmaSparseObj<T_data> *A, T_data *B) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_csr2dense needs a CSR matrix as input");
  }

  cusparseStatus_t status;
  cusparseHandle_t handle = A->current_context->get_cusparse_handle();

#if CUDA_VERSION < 11000
  status = sparse_dense<T_data>(handle, A->dims_data[1], A->dims_data[2], A->descr,
                     A->d_data, A->d_rowind, A->d_colind, B, A->dims_data[1]);
#else
  void* d_buffer = NULL;
  size_t bufferSize = 0;

  cusparseCreateDnMat(&(A->dn_descr), A->dims_data[1], A->dims_data[2], A->dims_data[2],
                      B, A->get_data_type(), CUSPARSE_ORDER_ROW);

  cusparseSparseToDense_bufferSize(handle, A->sp_descr, A->dn_descr,
                                        CUSPARSE_SPARSETODENSE_ALG_DEFAULT,
                                        &bufferSize);
  cudaMalloc(&d_buffer, bufferSize);
  status = cusparseSparseToDense(handle, A->sp_descr, A->dn_descr,
                                        CUSPARSE_SPARSETODENSE_ALG_DEFAULT,
                                        d_buffer);


  cusparseDestroyDnMat(A->dn_descr);
  cudaFree(d_buffer);

#endif

  if (status != CUSPARSE_STATUS_SUCCESS) {
    std::cerr << "Error | carma_csr2dense (sparse) | csr2dense failed"
              << std::endl;
    throw "Error | carma_csr2dense (sparse) | csr2dense failed";
  }
  return status;
}

template <class T_data>
cusparseStatus_t carma_csr2bsr(CarmaSparseObj<T_data> *A, int block_dim,
                                   CarmaSparseObj<T_data> *B) {
  if (A->format != "CSR") {
    DEBUG_TRACE("carma_csr2bsr needs a CSR matrix as input");
  }
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
  return sparse_bsr<T_data>(handle, dir, m, n, A->descr, A->d_data, A->d_rowind,
                 A->d_colind, B->block_dim, B->descr, B->d_data, B->d_rowind,
                 B->d_colind);
}

template <class T_data>
cusparseStatus_t carma_bsr2csr(CarmaSparseObj<T_data> *A,
                                   CarmaSparseObj<T_data> *B) {
  if (A->format != "BSR") {
    DEBUG_TRACE("carma_bsr2csr needs a BSR matrix as input");
  }
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
  return sparse_csr<T_data>(handle, dir, mb, nb, A->descr, A->d_data, A->d_rowind,
                 A->d_colind, A->block_dim, B->descr, B->d_data, B->d_rowind,
                 B->d_colind);
}

template cusparseStatus_t carma_gemv<float>(cusparseHandle_t handle,
                                            char op_A, float alpha,
                                            CarmaSparseObj<float> *A,
                                            float *x, float beta, float *y);

template cusparseStatus_t carma_gemv<double>(cusparseHandle_t handle,
                                             char op_A, double alpha,
                                             CarmaSparseObj<double> *A,
                                             double *x, double beta, double *y);


template cusparseStatus_t carma_gemm<float>(cusparseHandle_t handle,
                                            char op_A, float alpha,
                                            CarmaSparseObj<float> *A,
                                            CarmaObj<float> *B,
                                            float beta, CarmaObj<float> *C);

template cusparseStatus_t carma_gemm<double>(cusparseHandle_t handle,
                                             char op_A, double alpha,
                                             CarmaSparseObj<double> *A,
                                             CarmaObj<double> *B,
                                             double beta, CarmaObj<double> *C);


template cusparseStatus_t carma_gemm<float>(cusparseHandle_t handle,
                                            char op_A, char op_B,
                                            CarmaSparseObj<float> *A,
                                            CarmaSparseObj<float> *B,
                                            CarmaSparseObj<float> *C);

template cusparseStatus_t carma_gemm<double>(cusparseHandle_t handle,
                                             char op_A, char op_B,
                                             CarmaSparseObj<double> *A,
                                             CarmaSparseObj<double> *B,
                                             CarmaSparseObj<double> *C);

template cusparseStatus_t carma_csr2dense(CarmaSparseObj<float> *A, float *B);

template cusparseStatus_t carma_csr2dense(CarmaSparseObj<double> *A, double *B);


template cusparseStatus_t carma_csr2bsr(CarmaSparseObj<float> *src, int block_dim,
                                        CarmaSparseObj<float> *dest);

template cusparseStatus_t carma_csr2bsr(CarmaSparseObj<double> *src, int block_dim,
                                        CarmaSparseObj<double> *dest);



template cusparseStatus_t carma_bsr2csr(CarmaSparseObj<float> *src,
                                        CarmaSparseObj<float> *dest);

template cusparseStatus_t carma_bsr2csr(CarmaSparseObj<double> *src,
                                        CarmaSparseObj<double> *dest);
