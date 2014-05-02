#include <carma_cusparse.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>

#define carma_checkCusparseStatus(status) carma_checkCusparseStatus_v2(status, __LINE__, __FILE__)

cusparseStatus_t carma_checkCusparseStatus_v2(cusparseStatus_t status, int line,
    string file)
    /**< Generic CUSPARSE check status routine */
    {
  switch (status) {
  case CUSPARSE_STATUS_SUCCESS:
    return status;
  case CUSPARSE_STATUS_NOT_INITIALIZED:
    cerr << "!!!!! Cusparse error : CUBLAS_STATUS_NOT_INITIALIZED !!!!!\n";
    break;
  case CUSPARSE_STATUS_ALLOC_FAILED:
    cerr << "!!!!! Cusparse error : CUBLAS_STATUS_ALLOC_FAILED !!!!!\n";
    break;
  case CUSPARSE_STATUS_INVALID_VALUE:
    cerr << "!!!!! Cusparse error : CUBLAS_STATUS_ALLOC_FAILED !!!!!\n";
    break;
  case CUSPARSE_STATUS_ARCH_MISMATCH:
    cerr << "!!!!! Cusparse error : CUBLAS_STATUS_ARCH_MISMATCH !!!!!\n";
    break;
  case CUSPARSE_STATUS_MAPPING_ERROR:
    cerr << "!!!!! Cusparse error : CUBLAS_STATUS_MAPPING_ERROR !!!!!\n";
    break;
  case CUSPARSE_STATUS_EXECUTION_FAILED:
    cerr << "!!!!! Cusparse error : CUBLAS_STATUS_EXECUTION_FAILED !!!!!\n";
    break;
  case CUSPARSE_STATUS_INTERNAL_ERROR:
    cerr << "!!!!! Cusparse error : CUBLAS_STATUS_INTERNAL_ERROR !!!!!\n";
    break;
  case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
    cerr << "!!!!! Cusparse error : CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED !!!!!\n";
    break;
  }
  cerr << "!!!!! Cusparse error in " << file << "@" << line << " !!!!!\n";
  return status;
}

cusparseStatus_t carma_initCusparse(cusparseHandle_t *cusparse_handle)
/**< Generic CUSPARSE init routine */
{
  return carma_checkCusparseStatus(cusparseCreate(cusparse_handle));
}

cusparseStatus_t carma_shutdownCusparse(cusparseHandle_t cusparse_handle)
/**< Generic CUSPARSE shutdown routine */
{
  return carma_checkCusparseStatus(cusparseDestroy(cusparse_handle));
}

cusparseOperation_t carma_char2cusparseOperation(char operation) {
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


template<class T_data, cusparseStatus_t CUSPARSEAPI (*csrmv)(cusparseHandle_t handle,
    cusparseOperation_t transA,
    int m,
    int n,
    int nnz,
    const T_data *alpha,
    const cusparseMatDescr_t descrA,
    const T_data *csrValA,
    const int *csrRowPtrA,
    const int *csrColIndA,
    const T_data *x,
    const T_data *beta,
    T_data *y)>
cusparseStatus_t carma_gemv(cusparseHandle_t handle, char op_A,
    T_data alpha, carma_sparse_obj<T_data>* A, carma_obj<T_data>* x,
    T_data beta, carma_obj<T_data>* y){
  cusparseOperation_t trans=carma_char2cusparseOperation(op_A);

  return carma_checkCusparseStatus(csrmv(handle, trans, A->getDims(1), A->getDims(2), A->nz_elem, &alpha,
      A->descr, A->d_data, A->d_rowind, A->d_colind, *x, &beta, *y));

};

template<>
cusparseStatus_t carma_gemv<double>(cusparseHandle_t handle, char op_A, double alpha,
    carma_sparse_obj<double>* A, carma_obj<double>* x, double beta,
    carma_obj<double>* y){
  return carma_gemv<double,cusparseDcsrmv>(handle, op_A, alpha, A, x, beta, y);
}
template<>
cusparseStatus_t carma_gemv<float>(cusparseHandle_t handle, char op_A, float alpha,
    carma_sparse_obj<float>* A, carma_obj<float>* x, float beta,
    carma_obj<float>* y){
  return carma_gemv<float,cusparseScsrmv>(handle, op_A, alpha, A, x, beta, y);
}


template<class T_data, cusparseStatus_t CUSPARSEAPI (*csrmm)(cusparseHandle_t handle,
    cusparseOperation_t transA,
    int m,
    int n,
    int k,
    int nnz,
    const T_data *alpha,
    const cusparseMatDescr_t descrA,
    const T_data  *csrValA,
    const int *csrRowPtrA,
    const int *csrColIndA,
    const T_data *B,
    int ldb,
    const T_data *beta,
    T_data *C,
    int ldc)>
cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, T_data alpha,
    carma_sparse_obj<T_data>* A, carma_obj<T_data>* B, T_data beta,
    carma_obj<T_data>* C){
  //ofstream fichier;

  cusparseOperation_t transa = carma_char2cusparseOperation(op_A);
  cusparseStatus_t status;

  status = csrmm(handle, transa, A->getDims(1), B->getDims(2), A->getDims(2), A->nz_elem,
      &alpha, A->descr, A->d_data, A->d_rowind, A->d_colind, *B, B->getDims(1),
      &beta, *C, C->getDims(1));

  if (status != CUSPARSE_STATUS_SUCCESS) {
    cerr << "Error | carma_gemm (sparse) | Matrix-matrix multiplication failed"
        << endl;
    throw "Error | carma_gemm (sparse) | Matrix-matrix multiplication failed";
    //exit(EXIT_FAILURE);
  }
  return status;
}

template<>
cusparseStatus_t carma_gemm<float>(cusparseHandle_t handle, char op_A, float alpha,
    carma_sparse_obj<float>* A, carma_obj<float>* B, float beta,
    carma_obj<float>* C){
  return carma_gemm<float, cusparseScsrmm>(handle, op_A, alpha, A, B, beta, C);
}

template<>
cusparseStatus_t carma_gemm<double>(cusparseHandle_t handle, char op_A, double alpha,
    carma_sparse_obj<double>* A, carma_obj<double>* B, double beta,
    carma_obj<double>* C){
  return carma_gemm<double, cusparseDcsrmm>(handle, op_A, alpha, A, B, beta, C);
}
