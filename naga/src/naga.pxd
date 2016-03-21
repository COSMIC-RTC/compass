
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

    ctypedef enum cublasStatus_t:
        CUBLAS_STATUS_SUCCESS         =0
        CUBLAS_STATUS_NOT_INITIALIZED =1
        CUBLAS_STATUS_ALLOC_FAILED    =3
        CUBLAS_STATUS_INVALID_VALUE   =7
        CUBLAS_STATUS_ARCH_MISMATCH   =8
        CUBLAS_STATUS_MAPPING_ERROR   =11
        CUBLAS_STATUS_EXECUTION_FAILED=13
        CUBLAS_STATUS_INTERNAL_ERROR  =14
        CUBLAS_STATUS_NOT_SUPPORTED   =15
        CUBLAS_STATUS_LICENSE_ERROR   =16


    ctypedef enum cublasFillMode_t:
            CUBLAS_FILL_MODE_LOWER = 0
            CUBLAS_FILL_MODE_UPPER = 1
    ctypedef enum cublasSideMode_t:
            CUBLAS_SIDE_LEFT = 0
            CUBLAS_SIDE_RIGHT = 1

    struct cublasContext:
        pass
    ctypedef cublasContext *cublasHandle_t


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

cdef extern from "cusparse.h":
    ctypedef struct cusparseMatDescr_t:
        pass
    ctypedef struct cusparseStatus_t:
        pass
    ctypedef struct cusparseHandle_t:
        pass

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
        carma_device* get_device(int dev)
        cublasHandle_t get_cublasHandle()
        @staticmethod
        carma_context* instance()

#################################################
# P-Class naga_context
#################################################
cdef class naga_context:
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

    cdef cppclass tuple_t[T]:
        pass

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
        long getDims (int i)                  
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

        int reset()

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
# P-Classes naga_obj
#################################################
cdef class naga_obj_Int1D:
    cdef carma_obj[int] *c_o
    cdef naga_context context
cdef class naga_obj_Int2D:
    cdef carma_obj[int] *c_o
    cdef naga_context context
cdef class naga_obj_Int3D:
    cdef carma_obj[int] *c_o
    cdef naga_context context
cdef class naga_obj_Int4D:
    cdef carma_obj[int] *c_o
    cdef naga_context context

cdef class naga_obj_UInt1D:
    cdef carma_obj[unsigned int] *c_o
    cdef naga_context context
cdef class naga_obj_UInt2D:
    cdef carma_obj[unsigned int] *c_o
    cdef naga_context context
cdef class naga_obj_UInt3D:
    cdef carma_obj[unsigned int] *c_o
    cdef naga_context context
cdef class naga_obj_UInt4D:
    cdef carma_obj[unsigned int] *c_o
    cdef naga_context context

cdef class naga_obj_Float1D:
    cdef carma_obj[float] *c_o
    cdef naga_context context
cdef class naga_obj_Float2D:
    cdef carma_obj[float] *c_o
    cdef naga_context context
cdef class naga_obj_Float3D:
    cdef carma_obj[float] *c_o
    cdef naga_context context
cdef class naga_obj_Float4D:
    cdef carma_obj[float] *c_o
    cdef naga_context context

cdef class naga_obj_Double1D:
    cdef carma_obj[double] *c_o
    cdef naga_context context
cdef class naga_obj_Double2D:
    cdef carma_obj[double] *c_o
    cdef naga_context context
cdef class naga_obj_Double3D:
    cdef carma_obj[double] *c_o
    cdef naga_context context
cdef class naga_obj_Double4D:
    cdef carma_obj[double] *c_o
    cdef naga_context context

#use float2 and double2
# cuFloatComplex (cuDoubleComplex) dont go through compilation
cdef class naga_obj_ComplexS1D:
    cdef carma_obj[float2] *c_o
    cdef naga_context context
cdef class naga_obj_ComplexS2D:
    cdef carma_obj[float2] *c_o
    cdef naga_context context
cdef class naga_obj_ComplexS3D:
    cdef carma_obj[float2] *c_o
    cdef naga_context context
cdef class naga_obj_ComplexS4D:
    cdef carma_obj[float2] *c_o
    cdef naga_context context

cdef class naga_obj_ComplexD1D:
    cdef carma_obj[double2] *c_o
    cdef naga_context context
cdef class naga_obj_ComplexD2D:
    cdef carma_obj[double2] *c_o
    cdef naga_context context
cdef class naga_obj_ComplexD3D:
    cdef carma_obj[double2] *c_o
    cdef naga_context context
cdef class naga_obj_ComplexD4D:
    cdef carma_obj[double2] *c_o
    cdef naga_context context




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
# P-Classes naga_host_obj
#################################################
cdef class naga_host_obj_Int1D:
    cdef carma_host_obj[int] *c_h
    cdef int* data_h
cdef class naga_host_obj_Int2D:
    cdef carma_host_obj[int] *c_h
    cdef int* data_h
cdef class naga_host_obj_Int3D:
    cdef carma_host_obj[int] *c_h
    cdef int* data_h
cdef class naga_host_obj_Int4D:
    cdef carma_host_obj[int] *c_h
    cdef int* data_h


cdef class naga_host_obj_UInt1D:
    cdef carma_host_obj[unsigned int] *c_h
    cdef unsigned int* data_h
cdef class naga_host_obj_UInt2D:
    cdef carma_host_obj[unsigned int] *c_h
    cdef unsigned int* data_h
cdef class naga_host_obj_UInt3D:
    cdef carma_host_obj[unsigned int] *c_h
    cdef unsigned int* data_h
cdef class naga_host_obj_UInt4D:
    cdef carma_host_obj[unsigned int] *c_h
    cdef unsigned int* data_h


cdef class naga_host_obj_Float1D:
    cdef carma_host_obj[float] *c_h
    cdef float* data_h
cdef class naga_host_obj_Float2D:
    cdef carma_host_obj[float] *c_h
    cdef float* data_h
cdef class naga_host_obj_Float3D:
    cdef carma_host_obj[float] *c_h
    cdef float* data_h
cdef class naga_host_obj_Float4D:
    cdef carma_host_obj[float] *c_h
    cdef float* data_h


cdef class naga_host_obj_Double1D:
    cdef carma_host_obj[double] *c_h
    cdef double* data_h
cdef class naga_host_obj_Double2D:
    cdef carma_host_obj[double] *c_h
    cdef double* data_h
cdef class naga_host_obj_Double3D:
    cdef carma_host_obj[double] *c_h
    cdef double* data_h
cdef class naga_host_obj_Double4D:
    cdef carma_host_obj[double] *c_h
    cdef double* data_h


cdef class naga_host_obj_ComplexS1D:
    cdef carma_host_obj[float2] *c_h
    cdef float2* data_h
cdef class naga_host_obj_ComplexS2D:
    cdef carma_host_obj[float2] *c_h
    cdef float2* data_h
cdef class naga_host_obj_ComplexS3D:
    cdef carma_host_obj[float2] *c_h
    cdef float2* data_h
cdef class naga_host_obj_ComplexS4D:
    cdef carma_host_obj[float2] *c_h
    cdef float2* data_h


cdef class naga_host_obj_ComplexD1D:
    cdef carma_host_obj[double2] *c_h
    cdef double2* data_h
cdef class naga_host_obj_ComplexD2D:
    cdef carma_host_obj[double2] *c_h
    cdef double2* data_h
cdef class naga_host_obj_ComplexD3D:
    cdef carma_host_obj[double2] *c_h
    cdef double2* data_h
cdef class naga_host_obj_ComplexD4D:
    cdef carma_host_obj[double2] *c_h
    cdef double2* data_h


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



cdef extern from "carma_cublas.h":
    cublasStatus_t carma_axpy[T](cublasHandle_t cublas_handle, int n, T alpha, T *vectx,
        int incx, T *vecty, int incy)



#################################################
# C-Class carma_timer
#################################################
cdef extern from "carma_timer.h":
    cdef cppclass carma_timer:
        pass
        void start()
        void reset()
        void stop()
        double elapsed()

#################################################
# P-Class naga_timer
#################################################
cdef class naga_timer:
    cdef carma_timer *timer



cdef extern from "carma_utils.h":
    void carmaSafeDeviceSynchronize()
