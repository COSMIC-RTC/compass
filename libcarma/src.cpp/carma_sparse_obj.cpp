/*
 * carma_sparse_obj.cpp
 *
 *  Created on: Apr 8, 2014
 *      Author: ???
 */

#include "carma_sparse_obj.h"

template<class T_data>
carma_sparse_obj<T_data>::carma_sparse_obj() {
  _create(0, 0, 0);
}

template<class T_data>
carma_sparse_obj<T_data>::carma_sparse_obj(carma_obj<T_data>* M) {
  this->current_context = M->getContext();
  int *nnzPerRow;
  cudaMalloc((void**)&nnzPerRow, M->getDims(1)*sizeof(int));
  int nnzTotalDevHostPtr;

  cusparseMatDescr_t my_desc;
  cusparseCreateMatDescr(&my_desc);
  cusparseHandle_t handle = current_context->get_cusparseHandle();
  cusparseDnnz(handle, CUSPARSE_DIRECTION_ROW, M->getDims(1), M->getDims(2),
      my_desc, M->getData(), M->getDims(1), nnzPerRow, &nnzTotalDevHostPtr);

  _create(nnzTotalDevHostPtr, M->getDims(1), M->getDims(2));

  cusparseDdense2csr(handle, M->getDims(1), M->getDims(2), my_desc, M->getData(),
      M->getDims(1), nnzPerRow, this->d_data, this->csrRowPtr, this->d_colind);

  cudaFree(nnzPerRow);
}


template<class T_data>
carma_sparse_obj<T_data>::carma_sparse_obj(
    carma_sparse_obj<T_data>* M) {
  this->current_context=M->current_context;

  _create(M->nz_elem, M->getDims(1), M->getDims(2));

  cudaMemcpy(d_data, M->d_data, nz_elem*sizeof(T_data), cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_rowind, M->d_rowind, nz_elem*sizeof(int), cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_colind, M->d_colind, nz_elem*sizeof(int), cudaMemcpyDeviceToDevice);

  majorDim = M->majorDim;

  cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M->descr));
  cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M->descr));
  cusparseSetMatIndexBase(descr, cusparseGetMatIndexBase(M->descr));
  cusparseSetMatType(descr, cusparseGetMatType(M->descr));

}

template<class T_data>
carma_sparse_obj<T_data>::carma_sparse_obj(carma_context *current_context,
     carma_sparse_host_obj<T_data>* M) {

  this->device=current_context->get_activeDevice();
  this->current_context=current_context;

  _create(M->nz_elem, M->getDims(1), M->getDims(2));

  cudaMemcpy(d_data, M->h_data, nz_elem * sizeof(T_data),
      cudaMemcpyHostToDevice);
  cudaMemcpy(d_rowind, M->rowind, nz_elem * sizeof(int),
      cudaMemcpyHostToDevice);
  cudaMemcpy(d_colind, M->colind, nz_elem * sizeof(int),
      cudaMemcpyHostToDevice);

  majorDim = M->get_majorDim();

  cusparseHandle_t cusparseHandle=current_context->get_cusparseHandle();
  cusparseXcoo2csr(cusparseHandle, d_rowind,
      nz_elem, getDims(1), csrRowPtr, CUSPARSE_INDEX_BASE_ZERO);

}

template<class T_data>
void carma_sparse_obj<T_data>::operator=( carma_sparse_obj<T_data> &M) {
  resize(M.nz_elem, M.getDims(1), M.getDims(2));
  cudaMemcpy(d_data, M.d_data, nz_elem*sizeof(T_data), cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_rowind, M.d_rowind, nz_elem*sizeof(int), cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_colind, M.d_colind, nz_elem*sizeof(int), cudaMemcpyDeviceToDevice);

  majorDim = M.majorDim;

  cusparseSetMatDiagType(descr, cusparseGetMatDiagType(M.descr));
  cusparseSetMatFillMode(descr, cusparseGetMatFillMode(M.descr));
  cusparseSetMatIndexBase(descr, cusparseGetMatIndexBase(M.descr));
  cusparseSetMatType(descr, cusparseGetMatType(M.descr));
}

template<class T_data>
void carma_sparse_obj<T_data>::operator=(
    carma_sparse_host_obj<T_data> &M) {
  resize(M.nz_elem, M.getDims(1), M.getDims(2));
  cudaMemcpy(d_data, M.h_data, nz_elem * sizeof(T_data),
      cudaMemcpyHostToDevice);
  cudaMemcpy(d_rowind, M.rowind, nz_elem * sizeof(int),
      cudaMemcpyHostToDevice);
  cudaMemcpy(d_colind, M.colind, nz_elem * sizeof(int),
      cudaMemcpyHostToDevice);

  majorDim = M.get_majorDim();
}

template<class T_data>
void carma_sparse_obj<T_data>::resize(int nz_elem_, int dim1_, int dim2_) {
  if (nz_elem != nz_elem_) {
    _clear();
    _create(nz_elem_, dim1_, dim2_);
  } else {
    dims_data[0] = 2;
    dims_data[1] = dim1_;
    dims_data[2] = dim2_;
    majorDim = 'U';
  }

}

template<class T_data>
void carma_sparse_obj<T_data>::_create(int nz_elem_, int dim1_, int dim2_) {
  cusparseStatus_t status;
  nz_elem = nz_elem_;
  dims_data[0] = 2;
  dims_data[1] = dim1_;
  dims_data[2] = dim2_;

  if (nz_elem > 0) {
    cudaMalloc((void**) &d_data, nz_elem * sizeof(T_data));
    cudaMalloc((void**) &d_rowind, nz_elem * sizeof(int));
    cudaMalloc((void**) &d_colind, nz_elem * sizeof(int));
    cudaMalloc((void**) &csrRowPtr, (dim1_ + 1) * sizeof(int));
  }

  majorDim = 'U';

  status = cusparseCreateMatDescr(&descr);
  if (status != CUSPARSE_STATUS_SUCCESS) {
    cerr
        << "Error | carma_sparse_obj<T_data>::_create | Matrix descriptor initialization failed"
        << endl;
    throw "Error | carma_sparse_obj<T_data>::_create | Matrix descriptor initialization failed";
    //exit(EXIT_FAILURE);

  }

}

template<class T_data>
void carma_sparse_obj<T_data>::_clear() {
  cusparseStatus_t status;

  if (nz_elem > 0) {
    if (csrRowPtr == NULL || d_data == NULL
        || d_rowind == NULL || d_colind == NULL) {
      cerr << "Error | carma_sparse_obj<T_data>::_clear | double clear" << endl;
      throw "Error | carma_sparse_obj<T_data>::_clear | double clear";
    }
    cudaFree(d_data);
    cudaFree (d_rowind);
    cudaFree (d_colind);
    cudaFree(csrRowPtr);
  }
  d_data = NULL;
  d_rowind = NULL;
  d_colind = NULL;
  csrRowPtr = NULL;
  nz_elem = 0;
  dims_data[0] = 2;
  dims_data[1] = 0;
  dims_data[2] = 0;
  majorDim = 'U';

  status = cusparseDestroyMatDescr(descr);

  descr = 0;

  if (status != CUSPARSE_STATUS_SUCCESS) {
    cerr
        << "Error | carma_sparse_obj<T_data>::_clear | Matrix descriptor destruction failed"
        << endl;
    throw "Error | carma_sparse_obj<T_data>::_clear | Matrix descriptor destruction failed";
    //exit(EXIT_FAILURE);
  }
}

template<class T_data>
void carma_sparse_obj<T_data>::init_from_transpose(
    carma_sparse_obj<T_data>* M) {
  resize(M->nz_elem, M->getDims(1), M->getDims(2));

  cudaMemcpy(d_data, M->d_data, nz_elem*sizeof(T_data), cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_rowind, M->d_rowind, nz_elem*sizeof(int), cudaMemcpyDeviceToDevice);
  cudaMemcpy(d_colind, M->d_colind, nz_elem*sizeof(int), cudaMemcpyDeviceToDevice);

  if (M->majorDim == 'C')
    majorDim = 'R';
  else if (M->majorDim == 'R')
    majorDim = 'C';
  else
    majorDim = 'U';

}

template<class T_data>
bool carma_sparse_obj<T_data>::isColumnMajor() {
  bool colMajor = true;
/* TODO: demerdé ça
  int colm1 = 0;
  int rowm1 = 0;

  carma_sparse_host_obj<T_data> A_tmp;
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

template<class T_data>
carma_sparse_obj<T_data>::~carma_sparse_obj<T_data>() {
  _clear();
}

#define EXPLICITE_TEMPLATE(T_data) template carma_sparse_obj<T_data>::carma_sparse_obj(); \
    template carma_sparse_obj<T_data>::carma_sparse_obj( \
        carma_sparse_obj<T_data>* M); \
    template carma_sparse_obj<T_data>::carma_sparse_obj( \
        carma_obj<T_data>* M); \
    template carma_sparse_obj<T_data>::carma_sparse_obj(carma_context *current_context, \
        carma_sparse_host_obj<T_data>* M); \
    template carma_sparse_obj<T_data>::~carma_sparse_obj<T_data>(); \
    template void carma_sparse_obj<T_data>::operator=( carma_sparse_obj<T_data> &M); \
    template void carma_sparse_obj<T_data>::operator=( \
        carma_sparse_host_obj<T_data> &M); \
    template void carma_sparse_obj<T_data>::resize(int nz_elem_, int dim1_, int dim2_); \
    template void carma_sparse_obj<T_data>::_create(int nz_elem_, int dim1_, int dim2_); \
    template void carma_sparse_obj<T_data>::_clear(); \
    template void carma_sparse_obj<T_data>::init_from_transpose(carma_sparse_obj<T_data>* M); \
    template bool carma_sparse_obj<T_data>::isColumnMajor();

EXPLICITE_TEMPLATE(double);

#undef EXPLICITE_TEMPLATE
