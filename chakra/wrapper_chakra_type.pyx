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
