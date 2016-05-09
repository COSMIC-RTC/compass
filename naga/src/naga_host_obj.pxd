cimport numpy as np

from libc.stdint cimport uintptr_t
cdef extern from "driver_types.h":
    ctypedef cudaStream_t
    ctypedef cudaEvent_t

    ctypedef enum cudaMemcpyKind:
            cudaMemcpyHostToHost = 0
            cudaMemcpyHostToDevice = 1
            cudaMemcpyDeviceToHost = 2
            cudaMemcpyDeviceToDevice = 3
            cudaMemcpyDefault = 4


from naga_context cimport * 
from naga_context import *
from naga_obj cimport * 
from naga_obj import *
from naga_streams cimport * 
from naga_streams import *

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

        T * h_data  # < Input data
        T * data_UA  # < unpadded input dara for generic pinned mem
        long * dims_data  # < dimensions of the array
        int nb_elem  # < number of elments in the array
        MemAlloc mallocType  # < type of host alloc
        carma_streams * streams


        carma_host_obj(const long * dims_data)
        carma_host_obj(const long * dims_data, MemAlloc mallocType)
        carma_host_obj(const carma_host_obj[T] * obj)
        carma_host_obj(const carma_host_obj[T] * obj, MemAlloc mallocType)
        carma_host_obj(const long * dims_data, T * data)
        carma_host_obj(const long * dims_data, T * data, MemAlloc mallocType)
        carma_host_obj(const long * dims_data, int nb_streams)
        carma_host_obj(const long * dims_data, MemAlloc mallocType, int nb_streams)
        carma_host_obj(const carma_host_obj[T] * obj, int nb_streams)
        carma_host_obj(const carma_host_obj[T] * obj, MemAlloc mallocType,
            int nb_streams)
        carma_host_obj(const long * dims_data, T * data, int nb_streams)
        carma_host_obj(const long * dims_data, T * data, MemAlloc mallocType,
            int nb_streams)
        # ~carma_host_obj()

        void get_devpntr(void ** pntr_dev)

        int get_nbStreams()
        int add_stream()
        int add_stream(int nb)
        int del_stream()
        int del_stream(int nb)
        cudaStream_t get_cudaStream_t(int stream)
        int wait_stream(int stream)
        int wait_all_streams()

        int cpy_obj(carma_obj[T] * caObj, cudaMemcpyKind flag)
        int cpy_obj(carma_obj[T] * caObj, cudaMemcpyKind flag, unsigned int stream)

        T * getData()
        T * getData(int index)
        const long * getDims()
        int getNbElem()

        int fill_from(const T * data)
        int fill_into(T * data)


#################################################
# P-Classes naga_host_obj
#################################################
cdef class naga_host_obj_Int1D:
    cdef carma_host_obj[int] * c_h
    cdef int * data_h
cdef class naga_host_obj_Int2D:
    cdef carma_host_obj[int] * c_h
    cdef int * data_h
cdef class naga_host_obj_Int3D:
    cdef carma_host_obj[int] * c_h
    cdef int * data_h
cdef class naga_host_obj_Int4D:
    cdef carma_host_obj[int] * c_h
    cdef int * data_h


cdef class naga_host_obj_UInt1D:
    cdef carma_host_obj[unsigned int] * c_h
    cdef unsigned int * data_h
cdef class naga_host_obj_UInt2D:
    cdef carma_host_obj[unsigned int] * c_h
    cdef unsigned int * data_h
cdef class naga_host_obj_UInt3D:
    cdef carma_host_obj[unsigned int] * c_h
    cdef unsigned int * data_h
cdef class naga_host_obj_UInt4D:
    cdef carma_host_obj[unsigned int] * c_h
    cdef unsigned int * data_h


cdef class naga_host_obj_Float1D:
    cdef carma_host_obj[float] * c_h
    cdef float * data_h
cdef class naga_host_obj_Float2D:
    cdef carma_host_obj[float] * c_h
    cdef float * data_h
cdef class naga_host_obj_Float3D:
    cdef carma_host_obj[float] * c_h
    cdef float * data_h
cdef class naga_host_obj_Float4D:
    cdef carma_host_obj[float] * c_h
    cdef float * data_h


cdef class naga_host_obj_Double1D:
    cdef carma_host_obj[double] * c_h
    cdef double * data_h
cdef class naga_host_obj_Double2D:
    cdef carma_host_obj[double] * c_h
    cdef double * data_h
cdef class naga_host_obj_Double3D:
    cdef carma_host_obj[double] * c_h
    cdef double * data_h
cdef class naga_host_obj_Double4D:
    cdef carma_host_obj[double] * c_h
    cdef double * data_h


cdef class naga_host_obj_ComplexS1D:
    cdef carma_host_obj[float2] * c_h
    cdef float2 * data_h
cdef class naga_host_obj_ComplexS2D:
    cdef carma_host_obj[float2] * c_h
    cdef float2 * data_h
cdef class naga_host_obj_ComplexS3D:
    cdef carma_host_obj[float2] * c_h
    cdef float2 * data_h
cdef class naga_host_obj_ComplexS4D:
    cdef carma_host_obj[float2] * c_h
    cdef float2 * data_h


cdef class naga_host_obj_ComplexD1D:
    cdef carma_host_obj[double2] * c_h
    cdef double2 * data_h
cdef class naga_host_obj_ComplexD2D:
    cdef carma_host_obj[double2] * c_h
    cdef double2 * data_h
cdef class naga_host_obj_ComplexD3D:
    cdef carma_host_obj[double2] * c_h
    cdef double2 * data_h
cdef class naga_host_obj_ComplexD4D:
    cdef carma_host_obj[double2] * c_h
    cdef double2 * data_h

