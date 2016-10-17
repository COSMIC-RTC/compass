cimport numpy as np

from libc.stdint cimport uintptr_t
from libcpp cimport bool

from naga_context cimport *
from naga_context import *
from naga_streams cimport *
from naga_streams import *

cdef extern from "vector_types.h":
    struct float2:
        pass
    struct double2:
        pass

cdef extern from "cuComplex.h":
    ctypedef float2 cuFloatComplex
    ctypedef double2 cuDoubleComplex

ctypedef extern int cufftHandle


cdef extern from "cublas_v2.h":
    pass

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


    ctypedef enum cublasFillMode_t:
            CUBLAS_FILL_MODE_LOWER = 0
            CUBLAS_FILL_MODE_UPPER = 1
    ctypedef enum cublasSideMode_t:
            CUBLAS_SIDE_LEFT = 0
            CUBLAS_SIDE_RIGHT = 1

    struct cublasContext:
        pass
    ctypedef cublasContext * cublasHandle_t



# cuda plan
cdef extern from "cufft.h":
    ctypedef enum cufftType:
        CUFFT_R2C = 42
        CUFFT_C2R = 44
        CUFFT_C2C = 41
        CUFFT_D2Z = 106
        CUFFT_Z2D = 108
        CUFFT_Z2Z = 105


#################################################
# C-Class carma_obj
#################################################
cdef extern from "carma_obj.h":

    cdef cppclass tuple_t[T]:
        pass

    cdef cppclass carma_obj[T]:
        T * d_data
        int ndim
        long * dims_data
        int nb_elem
        int device
        carma_streams * stream

        carma_obj(carma_context * context, carma_obj[T] * obj)
        carma_obj(carma_obj[T] * obj)
        carma_obj(carma_context * context, const long * dims)
        carma_obj(carma_context * context, long * dims, T * data)

        T * getData()
        T * getData(int index)
        long * getDims()  #
        long getDims (int i)
        int getNbElem()  #
        carma_context * getContext()
        int getDevice()  #
        bool is_rng_init()  #

        int host2device(T * data)  #
        int device2host(T * data)  #
        int device2hostOpt(T * data)
        int host2deviceVect(T * data, int incx, int incy)
        int device2hostVect(T * data, int incx, int incy)
        int host2deviceMat(T * data, int lda, int ldb)
        int device2hostMat(T * data, int lda, int ldb)


        int copyInto(T * data, int nb_elem)  #
        int copyFrom(T * data, int nb_elem)  #

        int reset()

        cufftHandle * getPlan()
        unsigned int * getValues()
        T sum()  # (float and double only)

        int transpose(carma_obj[T] * source)  #

  # /**< Cublas V2 */
        int imax(int incx)  #
        int imin(int incx)  #
        T asum(int incx)  #
        T nrm2(int incx)  #
        T dot(carma_obj[T] * source, int incx, int incy)  #
        void scale(T alpha, int incx)  #
        void swap(carma_obj[T] * source, int incx, int incy)  #
        void copy(carma_obj[T] * source, int incx, int incy)  #
        void axpy(T alpha, carma_obj[T] * source, int incx, int incy)  #

# int carma_axpy_cpu<float>(long N, float alpha, float *h_X, long incX,
#    float *h_Y, long incY) {


#        void rot(carma_obj[T] *source, int incx, int incy, T sc, T ss)

        void gemv(char trans, T alpha, carma_obj[T] * matA, int lda,
            carma_obj[T] * vectx, int incx, T beta, int incy)  #
        void ger(T alpha, carma_obj[T] * vectx, int incx,
            carma_obj[T] * vecty, int incy, int lda)  #
        void symv(cublasFillMode_t uplo, T alpha, carma_obj[T] * matA,
            int lda, carma_obj[T] * vectx, int incx, T beta, int incy)  #

        void gemm(char transa, char transb, T alpha, carma_obj[T] * matA,  #
            int lda, carma_obj[T] * matB, int ldb, T beta, int ldc)
        void symm(cublasSideMode_t side, cublasFillMode_t uplo, T alpha,
            carma_obj[T] * matA, int lda, carma_obj[T] * matB, int ldb,
            T beta, int ldc)  #
        void syrk(cublasFillMode_t uplo, char transa, T alpha,
            carma_obj[T] * matA, int lda, T beta, int ldc)  #
        void syrkx(cublasFillMode_t uplo, char transa, T alpha,
            carma_obj[T] * matA, int lda, carma_obj[T] * matB, int ldb,
            T beta, int ldc)  #
        void geam(char transa, char transb, T alpha, carma_obj[T] * matA,
            int lda, T beta, carma_obj[T] * matB, int ldb, int ldc)  #
        void dgmm(cublasSideMode_t side, carma_obj[T] * matA, int lda,
            carma_obj[T] * vectx, int incx, int ldc)  #

        # /**< Curand */
        int init_prng(int device)
#        int destroy_prng()
#        int prng(T *output, char gtype, float alpha, float beta)
#        int prng(T *output, char gtype, float alpha)
#        int prng(char gtype, float alpha, float beta)
#        int prng(char gtype, float alpha)
        int prng(char gtype)
        int init_prng_host(int seed)
        int prng_host(char gtype)
#        int prng_host(char gtype, T alpha)
#        int destroy_prng_host()

        int prng_montagn( float init_montagn )



    cdef void  carma_initfft[T_in, T_out](const long * dims_data, cufftHandle * plan, cufftType tPlan)
    cdef int carma_fft[T_in, T_out](T_in * input, T_out * output, int dir, cufftHandle plan)

    cdef int snapTransformSize(unsigned int dataSize)
    cdef int carma_initfftconv(carma_obj[float] * data_in, carma_obj[float] * kernel_in,
        carma_obj[float] * padded_data, carma_obj[cuFloatComplex] * padded_spectrum, int kernelY, int kernelX)
    cdef int carma_fftconv(carma_obj[float] * data_out, carma_obj[float] * padded_data,
        carma_obj[cuFloatComplex] * padded_spectrum, int kernelY, int kernelX)


#################################################
# P-Classes naga_obj
#################################################
cdef class naga_obj_Int1D:
    cdef carma_obj[int] * c_o
    cdef naga_context context
cdef class naga_obj_Int2D:
    cdef carma_obj[int] * c_o
    cdef naga_context context
cdef class naga_obj_Int3D:
    cdef carma_obj[int] * c_o
    cdef naga_context context
cdef class naga_obj_Int4D:
    cdef carma_obj[int] * c_o
    cdef naga_context context

cdef class naga_obj_UInt1D:
    cdef carma_obj[unsigned int] * c_o
    cdef naga_context context
cdef class naga_obj_UInt2D:
    cdef carma_obj[unsigned int] * c_o
    cdef naga_context context
cdef class naga_obj_UInt3D:
    cdef carma_obj[unsigned int] * c_o
    cdef naga_context context
cdef class naga_obj_UInt4D:
    cdef carma_obj[unsigned int] * c_o
    cdef naga_context context

cdef class naga_obj_Float1D:
    cdef carma_obj[float] * c_o
    cdef naga_context context
cdef class naga_obj_Float2D:
    cdef carma_obj[float] * c_o
    cdef naga_context context
cdef class naga_obj_Float3D:
    cdef carma_obj[float] * c_o
    cdef naga_context context
cdef class naga_obj_Float4D:
    cdef carma_obj[float] * c_o
    cdef naga_context context

cdef class naga_obj_Double1D:
    cdef carma_obj[double] * c_o
    cdef naga_context context
cdef class naga_obj_Double2D:
    cdef carma_obj[double] * c_o
    cdef naga_context context
cdef class naga_obj_Double3D:
    cdef carma_obj[double] * c_o
    cdef naga_context context
cdef class naga_obj_Double4D:
    cdef carma_obj[double] * c_o
    cdef naga_context context

# use float2 and double2
# cuFloatComplex (cuDoubleComplex) dont go through compilation
cdef class naga_obj_ComplexS1D:
    cdef carma_obj[float2] * c_o
    cdef naga_context context
cdef class naga_obj_ComplexS2D:
    cdef carma_obj[float2] * c_o
    cdef naga_context context
cdef class naga_obj_ComplexS3D:
    cdef carma_obj[float2] * c_o
    cdef naga_context context
cdef class naga_obj_ComplexS4D:
    cdef carma_obj[float2] * c_o
    cdef naga_context context

cdef class naga_obj_ComplexD1D:
    cdef carma_obj[double2] * c_o
    cdef naga_context context
cdef class naga_obj_ComplexD2D:
    cdef carma_obj[double2] * c_o
    cdef naga_context context
cdef class naga_obj_ComplexD3D:
    cdef carma_obj[double2] * c_o
    cdef naga_context context
cdef class naga_obj_ComplexD4D:
    cdef carma_obj[double2] * c_o
    cdef naga_context context
