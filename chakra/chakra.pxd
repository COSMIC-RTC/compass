
cimport numpy as np
from libc.stdint cimport uintptr_t
from libcpp cimport bool

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
    ctypedef enum cublasFillMode_t:
            CUBLAS_FILL_MODE_LOWER = 0
            CUBLAS_FILL_MODE_UPPER = 1
    ctypedef enum cublasSideMode_t:
            CUBLAS_SIDE_LEFT = 0
            CUBLAS_SIDE_RIGHT = 1

#cuda plan
cdef extern from "cufft.h":
    ctypedef enum cufftType:
        CUFFT_R2C=42
        CUFFT_C2R=44
        CUFFT_C2C=41
        CUFFT_D2Z=106
        CUFFT_Z2D=108
        CUFFT_Z2Z=105


cdef extern from "driver_types.h":
    ctypedef cudaStream_t
    ctypedef cudaEvent_t

    ctypedef enum cudaMemcpyKind:
            cudaMemcpyHostToHost = 0
            cudaMemcpyHostToDevice =1
            cudaMemcpyDeviceToHost = 2
            cudaMemcpyDeviceToDevice = 3
            cudaMemcpyDefault = 4



#################################################
# C-Class carma_context
#################################################
cdef extern from "carma_context.h":
    cdef cppclass carma_device:
        pass

    cdef cppclass carma_context :
        #carma_context()
        int get_ndevice()
        int set_activeDevice(int newDevice, int silent)
        int set_activeDeviceForce(int newDevice, int silent)
        int set_activeDeviceForCpy(int newDevice, int silent)
        int get_activeDevice()
        @staticmethod
        carma_context* instance()

#################################################
# P-Class chakra_context
#################################################
cdef class chakra_context:
    cdef carma_context* c


#################################################
# C-Class carma_streams
#################################################


cdef extern from "carma_streams.h":
    cdef cppclass carma_streams:
        int eventflags
        carma_streams()
        carma_streams(unsigned int nbStreams)
        #~carma_streams()

        cudaStream_t get_stream(int stream)
        cudaEvent_t get_event(int stream)
        cudaStream_t operator[](int idx)

        int get_nbStreams()
        int add_stream()
        int add_stream(int nb)
        int del_stream()
        int del_stream(int nb)
        int del_all_streams()
        int wait_event(int stream)
        int wait_stream(int stream)
        int wait_all_streams()



#################################################
# C-Class carma_obj
#################################################
cdef extern from "carma_obj.h":
    cdef cppclass carma_obj[T]:
        T *d_data
        int ndim
        long *dims_data
        int nb_elem
        int device
        carma_streams *stream

        carma_obj(carma_context *context,carma_obj[T] *obj)
        carma_obj(carma_obj[T] *obj)
        carma_obj(carma_context *context, const long *dims)
        carma_obj(carma_context *context,long *dims, T *data)

        T* getData()                        
        T* getData(int index)                        
        long *getDims()                      #
        long getDimsi "getDims"(int i)                  
        int getNbElem()                     #
        carma_context *getContext()         
        int getDevice()                     #
        bool is_rng_init()                  #

        int host2device(T *data)            #
        int device2host(T *data)            #
        int device2hostOpt(T *data)
        int host2deviceVect(T *data, int incx, int incy)
        int device2hostVect(T *data, int incx, int incy)
        int host2deviceMat(T *data, int lda, int ldb)
        int device2hostMat(T *data, int lda, int ldb)


        int copyInto(T *data, int nb_elem)  #
        int copyFrom(T *data, int nb_elem)  #

        cufftHandle* getPlan()
        unsigned int * getValues()
        T sum()                   # (float and double only)        

        int transpose(carma_obj[T] *source) #

  #/**< Cublas V2 */
        int imax(int incx)  #
        int imin(int incx)  #
        T asum(int incx)    #
        T nrm2(int incx)    #
        T dot(carma_obj[T] *source, int incx, int incy)             #
        void scale(T alpha, int incx)                               #
        void swap(carma_obj[T] *source, int incx, int incy)         #
        void copy(carma_obj[T] *source, int incx, int incy)         #
        void axpy(T alpha, carma_obj[T] *source, int incx, int incy)#

#int carma_axpy_cpu<float>(long N, float alpha, float *h_X, long incX,
#    float *h_Y, long incY) {


#        void rot(carma_obj[T] *source, int incx, int incy, T sc, T ss)

        void gemv(char trans, T alpha, carma_obj[T] *matA, int lda,
            carma_obj[T] *vectx, int incx, T beta, int incy)                #
        void ger(T alpha, carma_obj[T] *vectx, int incx,
            carma_obj[T] *vecty, int incy, int lda)                         #
        void symv(cublasFillMode_t uplo, T alpha, carma_obj[T] *matA,
            int lda, carma_obj[T] *vectx, int incx, T beta, int incy)       #

        void gemm(char transa, char transb, T alpha, carma_obj[T] *matA,    #
            int lda, carma_obj[T] *matB, int ldb, T beta, int ldc)
        void symm(cublasSideMode_t side, cublasFillMode_t uplo, T alpha,
            carma_obj[T] *matA, int lda, carma_obj[T] *matB, int ldb,
            T beta, int ldc)                                                #
        void syrk(cublasFillMode_t uplo, char transa, T alpha,
            carma_obj[T] *matA, int lda, T beta, int ldc)                   #
        void syrkx(cublasFillMode_t uplo, char transa, T alpha,
            carma_obj[T] *matA, int lda, carma_obj[T] *matB, int ldb,
            T beta, int ldc)                                                #
        void geam(char transa, char transb, T alpha, carma_obj[T] *matA,
            int lda, T beta, carma_obj[T] *matB, int ldb, int ldc)          #
        void dgmm(cublasSideMode_t side, carma_obj[T] *matA, int lda,
            carma_obj[T] *vectx, int incx, int ldc)                         #

        #/**< Curand */
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



    cdef void  carma_initfft[T_in, T_out](const long *dims_data, cufftHandle *plan, cufftType tPlan)
    cdef int carma_fft[T_in, T_out](T_in *input, T_out *output, int dir, cufftHandle plan)

    cdef int snapTransformSize(unsigned int dataSize)
    cdef int carma_initfftconv(carma_obj[float] *data_in, carma_obj[float] *kernel_in, 
        carma_obj[float] *padded_data, carma_obj[cuFloatComplex] *padded_spectrum, int kernelY, int kernelX)
    cdef int carma_fftconv(carma_obj[float] *data_out, carma_obj[float] *padded_data,
        carma_obj[cuFloatComplex] *padded_spectrum, int kernelY, int kernelX)


#################################################
# P-Classes chakra_obj
#################################################
cdef class chakra_obj_Int1D:
    cdef carma_obj[int] *c_o
    cdef chakra_context context
cdef class chakra_obj_Int2D:
    cdef carma_obj[int] *c_o
    cdef chakra_context context
cdef class chakra_obj_Int3D:
    cdef carma_obj[int] *c_o
    cdef chakra_context context
cdef class chakra_obj_Int4D:
    cdef carma_obj[int] *c_o
    cdef chakra_context context

cdef class chakra_obj_UInt1D:
    cdef carma_obj[unsigned int] *c_o
    cdef chakra_context context
cdef class chakra_obj_UInt2D:
    cdef carma_obj[unsigned int] *c_o
    cdef chakra_context context
cdef class chakra_obj_UInt3D:
    cdef carma_obj[unsigned int] *c_o
    cdef chakra_context context
cdef class chakra_obj_UInt4D:
    cdef carma_obj[unsigned int] *c_o
    cdef chakra_context context

cdef class chakra_obj_Float1D:
    cdef carma_obj[float] *c_o
    cdef chakra_context context
cdef class chakra_obj_Float2D:
    cdef carma_obj[float] *c_o
    cdef chakra_context context
cdef class chakra_obj_Float3D:
    cdef carma_obj[float] *c_o
    cdef chakra_context context
cdef class chakra_obj_Float4D:
    cdef carma_obj[float] *c_o
    cdef chakra_context context

cdef class chakra_obj_Double1D:
    cdef carma_obj[double] *c_o
    cdef chakra_context context
cdef class chakra_obj_Double2D:
    cdef carma_obj[double] *c_o
    cdef chakra_context context
cdef class chakra_obj_Double3D:
    cdef carma_obj[double] *c_o
    cdef chakra_context context
cdef class chakra_obj_Double4D:
    cdef carma_obj[double] *c_o
    cdef chakra_context context

#use float2 and double2
# cuFloatComplex (cuDoubleComplex) dont go through compilation
cdef class chakra_obj_ComplexS1D:
    cdef carma_obj[float2] *c_o
    cdef chakra_context context
cdef class chakra_obj_ComplexS2D:
    cdef carma_obj[float2] *c_o
    cdef chakra_context context
cdef class chakra_obj_ComplexS3D:
    cdef carma_obj[float2] *c_o
    cdef chakra_context context
cdef class chakra_obj_ComplexS4D:
    cdef carma_obj[float2] *c_o
    cdef chakra_context context

cdef class chakra_obj_ComplexD1D:
    cdef carma_obj[double2] *c_o
    cdef chakra_context context
cdef class chakra_obj_ComplexD2D:
    cdef carma_obj[double2] *c_o
    cdef chakra_context context
cdef class chakra_obj_ComplexD3D:
    cdef carma_obj[double2] *c_o
    cdef chakra_context context
cdef class chakra_obj_ComplexD4D:
    cdef carma_obj[double2] *c_o
    cdef chakra_context context




#################################################
# C-Class carma_host_obj
#################################################
cdef extern from "carma_host_obj.h":
    ctypedef enum MemAlloc:
        MA_MALLOC
        MA_PAGELOCK
        MA_ZEROCPY
        MA_PORTABLE
        MA_WRICOMB
        MA_GENEPIN


    cdef cppclass carma_host_obj[T]:

        T *h_data #< Input data
        T *data_UA #< unpadded input dara for generic pinned mem
        long *dims_data #< dimensions of the array
        int nb_elem #< number of elments in the array
        MemAlloc mallocType #< type of host alloc
        carma_streams *streams


        carma_host_obj(const long *dims_data)
        carma_host_obj(const long *dims_data, MemAlloc mallocType)
        carma_host_obj(const carma_host_obj[T] *obj)
        carma_host_obj(const carma_host_obj[T] *obj, MemAlloc mallocType)
        carma_host_obj(const long *dims_data, T *data)
        carma_host_obj(const long *dims_data, T *data, MemAlloc mallocType)
        carma_host_obj(const long *dims_data, int nb_streams)
        carma_host_obj(const long *dims_data, MemAlloc mallocType, int nb_streams)
        carma_host_obj(const carma_host_obj[T] *obj, int nb_streams)
        carma_host_obj(const carma_host_obj[T] *obj, MemAlloc mallocType,
            int nb_streams)
        carma_host_obj(const long *dims_data, T *data, int nb_streams)
        carma_host_obj(const long *dims_data, T *data, MemAlloc mallocType,
            int nb_streams)
        #~carma_host_obj()

        void get_devpntr(void **pntr_dev)

        int get_nbStreams()
        int add_stream()
        int add_stream(int nb)
        int del_stream()
        int del_stream(int nb)
        cudaStream_t get_cudaStream_t(int stream)
        int wait_stream(int stream)
        int wait_all_streams()

        int cpy_obj(carma_obj[T]* caObj, cudaMemcpyKind flag)
        int cpy_obj(carma_obj[T]* caObj, cudaMemcpyKind flag, unsigned int stream)

        T* getData()
        T* getData(int index)
        const long *getDims()
        int getNbElem()

        int fill_from(const T *data)
        int fill_into(T *data)


#################################################
# P-Classes chakra_host_obj
#################################################
cdef class chakra_host_obj_Int1D:
    cdef carma_host_obj[int] *c_h
    cdef int* data_h
cdef class chakra_host_obj_Int2D:
    cdef carma_host_obj[int] *c_h
    cdef int* data_h
cdef class chakra_host_obj_Int3D:
    cdef carma_host_obj[int] *c_h
    cdef int* data_h
cdef class chakra_host_obj_Int4D:
    cdef carma_host_obj[int] *c_h
    cdef int* data_h


cdef class chakra_host_obj_UInt1D:
    cdef carma_host_obj[unsigned int] *c_h
    cdef unsigned int* data_h
cdef class chakra_host_obj_UInt2D:
    cdef carma_host_obj[unsigned int] *c_h
    cdef unsigned int* data_h
cdef class chakra_host_obj_UInt3D:
    cdef carma_host_obj[unsigned int] *c_h
    cdef unsigned int* data_h
cdef class chakra_host_obj_UInt4D:
    cdef carma_host_obj[unsigned int] *c_h
    cdef unsigned int* data_h


cdef class chakra_host_obj_Float1D:
    cdef carma_host_obj[float] *c_h
    cdef float* data_h
cdef class chakra_host_obj_Float2D:
    cdef carma_host_obj[float] *c_h
    cdef float* data_h
cdef class chakra_host_obj_Float3D:
    cdef carma_host_obj[float] *c_h
    cdef float* data_h
cdef class chakra_host_obj_Float4D:
    cdef carma_host_obj[float] *c_h
    cdef float* data_h


cdef class chakra_host_obj_Double1D:
    cdef carma_host_obj[double] *c_h
    cdef double* data_h
cdef class chakra_host_obj_Double2D:
    cdef carma_host_obj[double] *c_h
    cdef double* data_h
cdef class chakra_host_obj_Double3D:
    cdef carma_host_obj[double] *c_h
    cdef double* data_h
cdef class chakra_host_obj_Double4D:
    cdef carma_host_obj[double] *c_h
    cdef double* data_h


cdef class chakra_host_obj_ComplexS1D:
    cdef carma_host_obj[float2] *c_h
    cdef float2* data_h
cdef class chakra_host_obj_ComplexS2D:
    cdef carma_host_obj[float2] *c_h
    cdef float2* data_h
cdef class chakra_host_obj_ComplexS3D:
    cdef carma_host_obj[float2] *c_h
    cdef float2* data_h
cdef class chakra_host_obj_ComplexS4D:
    cdef carma_host_obj[float2] *c_h
    cdef float2* data_h


cdef class chakra_host_obj_ComplexD1D:
    cdef carma_host_obj[double2] *c_h
    cdef double2* data_h
cdef class chakra_host_obj_ComplexD2D:
    cdef carma_host_obj[double2] *c_h
    cdef double2* data_h
cdef class chakra_host_obj_ComplexD3D:
    cdef carma_host_obj[double2] *c_h
    cdef double2* data_h
cdef class chakra_host_obj_ComplexD4D:
    cdef carma_host_obj[double2] *c_h
    cdef double2* data_h


#################################################
# C-  MAGMA SIGNATURES
#################################################
cdef extern from "carma_obj.h":
    int carma_svd[T](carma_obj[T] *imat, carma_obj[T] *eigenvals,
        carma_obj[T] *mod2act, carma_obj[T] *mes2mod)
    int carma_getri[T](carma_obj[T] *d_iA)
    int carma_potri[T](carma_obj[T] *d_iA)

    int carma_syevd[T](char jobz, long N, T *mat, T *eigenvals)

cdef extern from "carma_host_obj.h":
    int carma_svd_cpu[T](carma_host_obj[T] *imat,
        carma_host_obj[T] *eigenvals, carma_host_obj[T] *mod2act,
        carma_host_obj[T] *mes2mod)
    int carma_getri_cpu[T](long N, T *h_A)
    int carma_potri_cpu[T](long N, T *h_A)
    int carma_syevd_cpu[T](char jobz, int N, T *h_A, T *eigenvals)
