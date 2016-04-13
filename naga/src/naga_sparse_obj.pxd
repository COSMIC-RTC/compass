cimport numpy as np
from naga_context cimport *
from naga_context import *
from naga_obj cimport *
from naga_obj import *


cdef extern from "cusparse.h":
    ctypedef struct cusparseMatDescr_t:
        pass
    ctypedef struct cusparseStatus_t:
        pass
    ctypedef struct cusparseHandle_t:
        pass




#################################################
# C-Class carma_sparse_obj
#################################################
cdef extern from "carma_sparse_obj.h":

    cdef cppclass carma_sparse_obj[T_data]:
      long dims_data[3] #< dimensions of the array
      int nz_elem #< number of elements in the array
      int device #< device where the carma_obj is allocate
      carma_context *current_context

      # ZERO-BASED INDEXING CSR-FORMAT
      T_data* d_data #nz_elem elements
      int* d_rowind #dim1+1  elements
      int* d_colind #nz_elem elements
      cusparseMatDescr_t descr

      char majorDim
      char* format;
      int blockDim #blockdim for BSR format

      # Magma stuff
      #union {
      #    magma_d_sparse_matrix d_spMat;
      #    magma_s_sparse_matrix s_spMat;
      #}

      carma_sparse_obj(carma_context *current_context)
      carma_sparse_obj(carma_obj[T_data]* M)
      carma_sparse_obj(carma_sparse_obj[T_data]* M)
      #carma_sparse_obj(carma_context *current_context, carma_sparse_host_obj[T_data]* M)
      carma_sparse_obj(carma_context *current_context, const long *dims, T_data * M, bool loadFromHost)
      carma_sparse_obj(carma_context *current_context, const long *dims,
          T_data *values, int *colind, int *rowind, int nz, bool loadFromHost)

      #void operator=(carma_sparse_obj[T_data]& M)
      #void operator=(carma_sparse_host_obj[T_data]& M)

      void resize(int nnz_, int dim1_, int dim2_)
      void init_from_transpose(carma_sparse_obj[T_data]* M)
      bool isColumnMajor()
      char get_majorDim() 
      void set_majorDim(char c)

      #  /**< General Utilities */
      #  operator T_data*() {
      #    return d_data;
      #  }
      #  T_data* operator[](int index) {
      #    return &d_data[index];
      #      }
      T_data* getData()
      T_data* getData(int index) 
      const long *getDims() 
      long getDims(int i) 
      int getNzElem() 
      carma_context* getContext() 

      int getDevice()

      void sparse_to_host(int *h_rowInd, int *h_colInd, T_data *h_data)

      void _create(int nnz_, int dim1_, int dim2_)
      void _clear()

      #template<cusparseStatus_t CUSPARSEAPI (*ptr_nnz)(cusparseHandle_t handle,
      #    cusparseDirection_t dirA, int m, int n, const cusparseMatDescr_t descrA,
      #    const T_data *A, int lda, int *nnzPerRowCol, int *nnzTotalDevHostPtr),
      #    cusparseStatus_t CUSPARSEAPI (*ptr_dense2csr)(cusparseHandle_t handle,
      #        int m, int n, const cusparseMatDescr_t descrA, const T_data *A,
      #        int lda, const int *nnzPerRow, T_data *csrValA, int *csrRowPtrA,
      #        int *csrColIndA)>
      #void init_carma_sparse_obj(carma_context *current_context, const long *dims,
      #    T_data * M, bool loadFromHost);
      #};

      cusparseStatus_t carma_gemv(cusparseHandle_t handle, char op_A, T_data alpha,
          carma_sparse_obj[T_data]* A, carma_obj[T_data]* x, T_data beta,
          carma_obj[T_data]* y);
  
      cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, T_data alpha,
          carma_sparse_obj[T_data]* A, carma_obj[T_data]* B, T_data beta,
          carma_obj[T_data]* C);
  
      cusparseStatus_t carma_gemm(cusparseHandle_t handle, char op_A, char op_B,
          carma_sparse_obj[T_data]* A, carma_sparse_obj[T_data]* B,
          carma_sparse_obj[T_data]* C);
  
      cusparseStatus_t carma_csr2dense(carma_sparse_obj[T_data]* src, T_data *dest);
  
      cusparseStatus_t carma_csr2bsr(carma_sparse_obj[T_data]* src, int blockDim,
          carma_sparse_obj[T_data] *dest);
  
      cusparseStatus_t carma_bsr2csr(carma_sparse_obj[T_data]* src,
          carma_sparse_obj[T_data] *dest);
  
      int carma_magma_csr2ell(carma_sparse_obj[T_data] *dA);
  
      int carma_magma_spmv(T_data alpha, carma_sparse_obj[T_data] *dA, carma_obj[T_data] *dx, T_data beta, carma_obj[T_data] *dy);
  
      int carma_sparse_magma_free(carma_sparse_obj[T_data] *dA);
  
      #    int carma_kgemv(carma_sparse_obj[T_data]* A, T_data alpha,
      #    const T_data* __restrict x, T_data beta, T_data *y);

#################################################
# P-Classes naga_sparse_obj
#################################################
cdef class naga_sparse_obj_Float:
    cdef carma_sparse_obj[float] *c_sparse_obj
    cdef copy(self, carma_sparse_obj[float] *c_sparse_obj)

cdef class naga_sparse_obj_Double:
    cdef carma_sparse_obj[double] *c_sparse_obj
    cdef copy(self, carma_sparse_obj[double] *c_sparse_obj)



