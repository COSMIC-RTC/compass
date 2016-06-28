cimport numpy as np

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

#################################################
# C-Class carma_context
#################################################
cdef extern from "carma_context.h":
    cdef cppclass carma_device:
        pass

    cdef cppclass carma_context :
        # carma_context()
        int get_ndevice()
        int set_activeDevice(int newDevice, int silent)
        int set_activeDeviceForce(int newDevice, int silent)
        int set_activeDeviceForCpy(int newDevice, int silent)
        int get_activeDevice()
        carma_device * get_device(int dev)
        cublasHandle_t get_cublasHandle()
        @staticmethod
        carma_context& instance()
        @staticmethod
        carma_context& instance_1gpu(int num_device)
        @staticmethod
        carma_context& instance_ngpu(int nb_devices, np.int32_t * devices_id)

#################################################
# P-Class naga_context
#################################################
cdef class naga_context:
    cdef carma_context *c

