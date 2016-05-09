cimport numpy as np

from naga_obj cimport * 
from naga_obj import *
from naga_host_obj cimport * 
from naga_host_obj import *

cdef extern from "cublas_api.h":

    ctypedef enum cublasStatus_t:
        CUBLAS_STATUS_SUCCESS = 0
        CUBLAS_STATUS_NOT_INITIALIZED = 1
        CUBLAS_STATUS_ALLOC_FAILED = 3
        CUBLAS_STATUS_INVALID_VALUE = 7
        CUBLAS_STATUS_ARCH_MISMATCH = 8
        CUBLAS_STATUS_MAPPING_ERROR = 11
        CUBLAS_STATUS_EXECUTION_FAILED = 13
        CUBLAS_STATUS_INTERNAL_ERROR = 14
        CUBLAS_STATUS_NOT_SUPPORTED = 15
        CUBLAS_STATUS_LICENSE_ERROR = 16


#################################################
# C-  MAGMA SIGNATURES
#################################################
cdef extern from "carma_obj.h":
    int carma_svd[T](carma_obj[T] * imat, carma_obj[T] * eigenvals,
        carma_obj[T] * mod2act, carma_obj[T] * mes2mod)
    int carma_getri[T](carma_obj[T] * d_iA)
    int carma_potri[T](carma_obj[T] * d_iA)

    int carma_syevd[T](char jobz, long N, T * mat, T * eigenvals)

cdef extern from "carma_host_obj.h":
    int carma_svd_cpu[T](carma_host_obj[T] * imat,
        carma_host_obj[T] * eigenvals, carma_host_obj[T] * mod2act,
        carma_host_obj[T] * mes2mod)
    int carma_getri_cpu[T](long N, T * h_A)
    int carma_potri_cpu[T](long N, T * h_A)
    int carma_syevd_cpu[T](char jobz, int N, T * h_A, T * eigenvals)



cdef extern from "carma_cublas.h":
    cublasStatus_t carma_axpy[T](cublasHandle_t cublas_handle, int n, T alpha, T * vectx,
        int incx, T * vecty, int incy)

