

import numpy as np
cimport numpy as np
#np.import_array()

cdef dict_MemAlloc={"malloc":(MA_MALLOC),
                    "pagelock":(MA_PAGELOCK),
                    "zerocpy":(MA_ZEROCPY),
                    "portable":(MA_PORTABLE),
                    "wricomb":(MA_WRICOMB),
                    "genepin":(MA_GENEPIN),
                    "none":(-1)}
#none : add into C code


cdef class naga_host_obj_Int1D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Int1D obj=None,
            np.ndarray[ndim=1,dtype=np.int32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Int1D constructor.

    constructors available:
        naga_host_obj_Int1D(np.ndarray dims)
        naga_host_obj_Int1D(np.ndarray dims, int nbStreams)
        naga_host_obj_Int1D(np.ndarray dims, str mallocType)
        naga_host_obj_Int1D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Int1D(np.ndarray data)
        naga_host_obj_Int1D(np.ndarray data, int nbStream)
        naga_host_obj_Int1D(np.ndarray data, str mallocType)
        naga_host_obj_Int1D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Int1D(naga_host_obj_Int1D obj)
        naga_host_obj_Int1D(naga_host_obj_Int1D obj, int nbStream)
        naga_host_obj_Int1D(naga_host_obj_Int1D obj, str mallocType)
        naga_host_obj_Int1D(naga_host_obj_Int1D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=1
        dims must be a np array of shape (1)
        """

        cdef long sh[1+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(1):
                        sh[i+1]=dims[i]
                    sh[0]=1
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[int](sh)
                    else:
                        self.c_h=new carma_host_obj[int](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[int](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[int](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=1
                for i in range(1):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data)
                    else:
                        self.c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[int](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[int](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[int](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[int](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((1),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(1):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Int1D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Int1D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[int]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Int1D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Int1D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[int]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((1),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=1,dtype=np.int32_t] data

        cdims=self.c_h.getDims()

        for i in range(1):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.int32,order="F")

        self.c_h.fill_into(<int*>data.data)
        return data


    def setData(self,np.ndarray[ndim=1,dtype=np.int32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<int*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Int2D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Int2D obj=None,
            np.ndarray[ndim=2,dtype=np.int32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Int2D constructor.

    constructors available:
        naga_host_obj_Int2D(np.ndarray dims)
        naga_host_obj_Int2D(np.ndarray dims, int nbStreams)
        naga_host_obj_Int2D(np.ndarray dims, str mallocType)
        naga_host_obj_Int2D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Int2D(np.ndarray data)
        naga_host_obj_Int2D(np.ndarray data, int nbStream)
        naga_host_obj_Int2D(np.ndarray data, str mallocType)
        naga_host_obj_Int2D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Int2D(naga_host_obj_Int2D obj)
        naga_host_obj_Int2D(naga_host_obj_Int2D obj, int nbStream)
        naga_host_obj_Int2D(naga_host_obj_Int2D obj, str mallocType)
        naga_host_obj_Int2D(naga_host_obj_Int2D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=2
        dims must be a np array of shape (2)
        """

        cdef long sh[2+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(2):
                        sh[i+1]=dims[i]
                    sh[0]=2
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[int](sh)
                    else:
                        self.c_h=new carma_host_obj[int](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[int](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[int](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=2
                for i in range(2):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data)
                    else:
                        self.c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[int](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[int](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[int](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[int](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((2),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(2):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Int2D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Int2D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[int]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Int2D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Int2D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[int]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((2),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=2,dtype=np.int32_t] data

        cdims=self.c_h.getDims()

        for i in range(2):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.int32,order="F")

        self.c_h.fill_into(<int*>data.data)
        return data


    def setData(self,np.ndarray[ndim=2,dtype=np.int32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<int*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Int3D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Int3D obj=None,
            np.ndarray[ndim=3,dtype=np.int32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Int3D constructor.

    constructors available:
        naga_host_obj_Int3D(np.ndarray dims)
        naga_host_obj_Int3D(np.ndarray dims, int nbStreams)
        naga_host_obj_Int3D(np.ndarray dims, str mallocType)
        naga_host_obj_Int3D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Int3D(np.ndarray data)
        naga_host_obj_Int3D(np.ndarray data, int nbStream)
        naga_host_obj_Int3D(np.ndarray data, str mallocType)
        naga_host_obj_Int3D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Int3D(naga_host_obj_Int3D obj)
        naga_host_obj_Int3D(naga_host_obj_Int3D obj, int nbStream)
        naga_host_obj_Int3D(naga_host_obj_Int3D obj, str mallocType)
        naga_host_obj_Int3D(naga_host_obj_Int3D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=3
        dims must be a np array of shape (3)
        """

        cdef long sh[3+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(3):
                        sh[i+1]=dims[i]
                    sh[0]=3
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[int](sh)
                    else:
                        self.c_h=new carma_host_obj[int](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[int](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[int](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=3
                for i in range(3):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data)
                    else:
                        self.c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[int](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[int](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[int](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[int](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((3),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(3):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Int3D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Int3D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[int]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Int3D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Int3D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[int]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((3),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=3,dtype=np.int32_t] data

        cdims=self.c_h.getDims()

        for i in range(3):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.int32,order="F")

        self.c_h.fill_into(<int*>data.data)
        return data


    def setData(self,np.ndarray[ndim=3,dtype=np.int32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<int*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Int4D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Int4D obj=None,
            np.ndarray[ndim=4,dtype=np.int32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Int4D constructor.

    constructors available:
        naga_host_obj_Int4D(np.ndarray dims)
        naga_host_obj_Int4D(np.ndarray dims, int nbStreams)
        naga_host_obj_Int4D(np.ndarray dims, str mallocType)
        naga_host_obj_Int4D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Int4D(np.ndarray data)
        naga_host_obj_Int4D(np.ndarray data, int nbStream)
        naga_host_obj_Int4D(np.ndarray data, str mallocType)
        naga_host_obj_Int4D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Int4D(naga_host_obj_Int4D obj)
        naga_host_obj_Int4D(naga_host_obj_Int4D obj, int nbStream)
        naga_host_obj_Int4D(naga_host_obj_Int4D obj, str mallocType)
        naga_host_obj_Int4D(naga_host_obj_Int4D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=4
        dims must be a np array of shape (4)
        """

        cdef long sh[4+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(4):
                        sh[i+1]=dims[i]
                    sh[0]=4
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[int](sh)
                    else:
                        self.c_h=new carma_host_obj[int](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[int](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[int](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=4
                for i in range(4):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data)
                    else:
                        self.c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[int](sh,
                                                            <int*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[int](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[int](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[int](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[int](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((4),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(4):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Int4D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Int4D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[int]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Int4D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Int4D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[int]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((4),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=4,dtype=np.int32_t] data

        cdims=self.c_h.getDims()

        for i in range(4):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.int32,order="F")

        self.c_h.fill_into(<int*>data.data)
        return data


    def setData(self,np.ndarray[ndim=4,dtype=np.int32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<int*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_UInt1D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_UInt1D obj=None,
            np.ndarray[ndim=1,dtype=np.uint32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_UInt1D constructor.

    constructors available:
        naga_host_obj_UInt1D(np.ndarray dims)
        naga_host_obj_UInt1D(np.ndarray dims, int nbStreams)
        naga_host_obj_UInt1D(np.ndarray dims, str mallocType)
        naga_host_obj_UInt1D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_UInt1D(np.ndarray data)
        naga_host_obj_UInt1D(np.ndarray data, int nbStream)
        naga_host_obj_UInt1D(np.ndarray data, str mallocType)
        naga_host_obj_UInt1D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_UInt1D(naga_host_obj_UInt1D obj)
        naga_host_obj_UInt1D(naga_host_obj_UInt1D obj, int nbStream)
        naga_host_obj_UInt1D(naga_host_obj_UInt1D obj, str mallocType)
        naga_host_obj_UInt1D(naga_host_obj_UInt1D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=1
        dims must be a np array of shape (1)
        """

        cdef long sh[1+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(1):
                        sh[i+1]=dims[i]
                    sh[0]=1
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[unsigned int](sh)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[unsigned int](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=1
                for i in range(1):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((1),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(1):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_UInt1D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_UInt1D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[unsigned int]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_UInt1D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_UInt1D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[unsigned int]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((1),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=1,dtype=np.uint32_t] data

        cdims=self.c_h.getDims()

        for i in range(1):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.uint32,order="F")

        self.c_h.fill_into(<unsigned int*>data.data)
        return data


    def setData(self,np.ndarray[ndim=1,dtype=np.uint32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<unsigned int*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_UInt2D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_UInt2D obj=None,
            np.ndarray[ndim=2,dtype=np.uint32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_UInt2D constructor.

    constructors available:
        naga_host_obj_UInt2D(np.ndarray dims)
        naga_host_obj_UInt2D(np.ndarray dims, int nbStreams)
        naga_host_obj_UInt2D(np.ndarray dims, str mallocType)
        naga_host_obj_UInt2D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_UInt2D(np.ndarray data)
        naga_host_obj_UInt2D(np.ndarray data, int nbStream)
        naga_host_obj_UInt2D(np.ndarray data, str mallocType)
        naga_host_obj_UInt2D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_UInt2D(naga_host_obj_UInt2D obj)
        naga_host_obj_UInt2D(naga_host_obj_UInt2D obj, int nbStream)
        naga_host_obj_UInt2D(naga_host_obj_UInt2D obj, str mallocType)
        naga_host_obj_UInt2D(naga_host_obj_UInt2D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=2
        dims must be a np array of shape (2)
        """

        cdef long sh[2+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(2):
                        sh[i+1]=dims[i]
                    sh[0]=2
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[unsigned int](sh)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[unsigned int](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=2
                for i in range(2):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((2),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(2):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_UInt2D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_UInt2D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[unsigned int]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_UInt2D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_UInt2D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[unsigned int]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((2),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=2,dtype=np.uint32_t] data

        cdims=self.c_h.getDims()

        for i in range(2):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.uint32,order="F")

        self.c_h.fill_into(<unsigned int*>data.data)
        return data


    def setData(self,np.ndarray[ndim=2,dtype=np.uint32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<unsigned int*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_UInt3D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_UInt3D obj=None,
            np.ndarray[ndim=3,dtype=np.uint32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_UInt3D constructor.

    constructors available:
        naga_host_obj_UInt3D(np.ndarray dims)
        naga_host_obj_UInt3D(np.ndarray dims, int nbStreams)
        naga_host_obj_UInt3D(np.ndarray dims, str mallocType)
        naga_host_obj_UInt3D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_UInt3D(np.ndarray data)
        naga_host_obj_UInt3D(np.ndarray data, int nbStream)
        naga_host_obj_UInt3D(np.ndarray data, str mallocType)
        naga_host_obj_UInt3D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_UInt3D(naga_host_obj_UInt3D obj)
        naga_host_obj_UInt3D(naga_host_obj_UInt3D obj, int nbStream)
        naga_host_obj_UInt3D(naga_host_obj_UInt3D obj, str mallocType)
        naga_host_obj_UInt3D(naga_host_obj_UInt3D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=3
        dims must be a np array of shape (3)
        """

        cdef long sh[3+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(3):
                        sh[i+1]=dims[i]
                    sh[0]=3
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[unsigned int](sh)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[unsigned int](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=3
                for i in range(3):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((3),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(3):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_UInt3D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_UInt3D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[unsigned int]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_UInt3D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_UInt3D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[unsigned int]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((3),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=3,dtype=np.uint32_t] data

        cdims=self.c_h.getDims()

        for i in range(3):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.uint32,order="F")

        self.c_h.fill_into(<unsigned int*>data.data)
        return data


    def setData(self,np.ndarray[ndim=3,dtype=np.uint32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<unsigned int*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_UInt4D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_UInt4D obj=None,
            np.ndarray[ndim=4,dtype=np.uint32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_UInt4D constructor.

    constructors available:
        naga_host_obj_UInt4D(np.ndarray dims)
        naga_host_obj_UInt4D(np.ndarray dims, int nbStreams)
        naga_host_obj_UInt4D(np.ndarray dims, str mallocType)
        naga_host_obj_UInt4D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_UInt4D(np.ndarray data)
        naga_host_obj_UInt4D(np.ndarray data, int nbStream)
        naga_host_obj_UInt4D(np.ndarray data, str mallocType)
        naga_host_obj_UInt4D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_UInt4D(naga_host_obj_UInt4D obj)
        naga_host_obj_UInt4D(naga_host_obj_UInt4D obj, int nbStream)
        naga_host_obj_UInt4D(naga_host_obj_UInt4D obj, str mallocType)
        naga_host_obj_UInt4D(naga_host_obj_UInt4D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=4
        dims must be a np array of shape (4)
        """

        cdef long sh[4+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(4):
                        sh[i+1]=dims[i]
                    sh[0]=4
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[unsigned int](sh)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[unsigned int](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=4
                for i in range(4):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data)
                    else:
                        self.c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[unsigned int](sh,
                                                            <unsigned int*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[unsigned int](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((4),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(4):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_UInt4D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_UInt4D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[unsigned int]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_UInt4D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_UInt4D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[unsigned int]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((4),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=4,dtype=np.uint32_t] data

        cdims=self.c_h.getDims()

        for i in range(4):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.uint32,order="F")

        self.c_h.fill_into(<unsigned int*>data.data)
        return data


    def setData(self,np.ndarray[ndim=4,dtype=np.uint32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<unsigned int*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Float1D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Float1D obj=None,
            np.ndarray[ndim=1,dtype=np.float32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Float1D constructor.

    constructors available:
        naga_host_obj_Float1D(np.ndarray dims)
        naga_host_obj_Float1D(np.ndarray dims, int nbStreams)
        naga_host_obj_Float1D(np.ndarray dims, str mallocType)
        naga_host_obj_Float1D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Float1D(np.ndarray data)
        naga_host_obj_Float1D(np.ndarray data, int nbStream)
        naga_host_obj_Float1D(np.ndarray data, str mallocType)
        naga_host_obj_Float1D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Float1D(naga_host_obj_Float1D obj)
        naga_host_obj_Float1D(naga_host_obj_Float1D obj, int nbStream)
        naga_host_obj_Float1D(naga_host_obj_Float1D obj, str mallocType)
        naga_host_obj_Float1D(naga_host_obj_Float1D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=1
        dims must be a np array of shape (1)
        """

        cdef long sh[1+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(1):
                        sh[i+1]=dims[i]
                    sh[0]=1
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[float](sh)
                    else:
                        self.c_h=new carma_host_obj[float](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[float](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[float](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=1
                for i in range(1):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data)
                    else:
                        self.c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[float](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[float](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[float](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[float](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((1),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(1):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Float1D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Float1D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[float]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Float1D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Float1D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[float]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((1),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=1,dtype=np.float32_t] data

        cdims=self.c_h.getDims()

        for i in range(1):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.float32,order="F")

        self.c_h.fill_into(<float*>data.data)
        return data


    def setData(self,np.ndarray[ndim=1,dtype=np.float32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<float*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Float2D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Float2D obj=None,
            np.ndarray[ndim=2,dtype=np.float32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Float2D constructor.

    constructors available:
        naga_host_obj_Float2D(np.ndarray dims)
        naga_host_obj_Float2D(np.ndarray dims, int nbStreams)
        naga_host_obj_Float2D(np.ndarray dims, str mallocType)
        naga_host_obj_Float2D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Float2D(np.ndarray data)
        naga_host_obj_Float2D(np.ndarray data, int nbStream)
        naga_host_obj_Float2D(np.ndarray data, str mallocType)
        naga_host_obj_Float2D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Float2D(naga_host_obj_Float2D obj)
        naga_host_obj_Float2D(naga_host_obj_Float2D obj, int nbStream)
        naga_host_obj_Float2D(naga_host_obj_Float2D obj, str mallocType)
        naga_host_obj_Float2D(naga_host_obj_Float2D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=2
        dims must be a np array of shape (2)
        """

        cdef long sh[2+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(2):
                        sh[i+1]=dims[i]
                    sh[0]=2
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[float](sh)
                    else:
                        self.c_h=new carma_host_obj[float](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[float](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[float](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=2
                for i in range(2):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data)
                    else:
                        self.c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[float](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[float](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[float](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[float](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((2),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(2):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Float2D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Float2D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[float]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Float2D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Float2D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[float]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((2),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data

        cdims=self.c_h.getDims()

        for i in range(2):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.float32,order="F")

        self.c_h.fill_into(<float*>data.data)
        return data


    def setData(self,np.ndarray[ndim=2,dtype=np.float32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<float*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Float3D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Float3D obj=None,
            np.ndarray[ndim=3,dtype=np.float32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Float3D constructor.

    constructors available:
        naga_host_obj_Float3D(np.ndarray dims)
        naga_host_obj_Float3D(np.ndarray dims, int nbStreams)
        naga_host_obj_Float3D(np.ndarray dims, str mallocType)
        naga_host_obj_Float3D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Float3D(np.ndarray data)
        naga_host_obj_Float3D(np.ndarray data, int nbStream)
        naga_host_obj_Float3D(np.ndarray data, str mallocType)
        naga_host_obj_Float3D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Float3D(naga_host_obj_Float3D obj)
        naga_host_obj_Float3D(naga_host_obj_Float3D obj, int nbStream)
        naga_host_obj_Float3D(naga_host_obj_Float3D obj, str mallocType)
        naga_host_obj_Float3D(naga_host_obj_Float3D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=3
        dims must be a np array of shape (3)
        """

        cdef long sh[3+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(3):
                        sh[i+1]=dims[i]
                    sh[0]=3
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[float](sh)
                    else:
                        self.c_h=new carma_host_obj[float](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[float](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[float](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=3
                for i in range(3):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data)
                    else:
                        self.c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[float](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[float](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[float](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[float](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((3),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(3):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Float3D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Float3D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[float]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Float3D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Float3D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[float]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((3),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=3,dtype=np.float32_t] data

        cdims=self.c_h.getDims()

        for i in range(3):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.float32,order="F")

        self.c_h.fill_into(<float*>data.data)
        return data


    def setData(self,np.ndarray[ndim=3,dtype=np.float32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<float*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Float4D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Float4D obj=None,
            np.ndarray[ndim=4,dtype=np.float32_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Float4D constructor.

    constructors available:
        naga_host_obj_Float4D(np.ndarray dims)
        naga_host_obj_Float4D(np.ndarray dims, int nbStreams)
        naga_host_obj_Float4D(np.ndarray dims, str mallocType)
        naga_host_obj_Float4D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Float4D(np.ndarray data)
        naga_host_obj_Float4D(np.ndarray data, int nbStream)
        naga_host_obj_Float4D(np.ndarray data, str mallocType)
        naga_host_obj_Float4D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Float4D(naga_host_obj_Float4D obj)
        naga_host_obj_Float4D(naga_host_obj_Float4D obj, int nbStream)
        naga_host_obj_Float4D(naga_host_obj_Float4D obj, str mallocType)
        naga_host_obj_Float4D(naga_host_obj_Float4D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=4
        dims must be a np array of shape (4)
        """

        cdef long sh[4+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(4):
                        sh[i+1]=dims[i]
                    sh[0]=4
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[float](sh)
                    else:
                        self.c_h=new carma_host_obj[float](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[float](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[float](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=4
                for i in range(4):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data)
                    else:
                        self.c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[float](sh,
                                                            <float*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[float](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[float](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[float](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[float](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((4),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(4):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Float4D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Float4D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[float]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Float4D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Float4D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[float]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((4),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=4,dtype=np.float32_t] data

        cdims=self.c_h.getDims()

        for i in range(4):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.float32,order="F")

        self.c_h.fill_into(<float*>data.data)
        return data


    def setData(self,np.ndarray[ndim=4,dtype=np.float32_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<float*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_ComplexS1D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_ComplexS1D obj=None,
            np.ndarray[ndim=1,dtype=np.complex64_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_ComplexS1D constructor.

    constructors available:
        naga_host_obj_ComplexS1D(np.ndarray dims)
        naga_host_obj_ComplexS1D(np.ndarray dims, int nbStreams)
        naga_host_obj_ComplexS1D(np.ndarray dims, str mallocType)
        naga_host_obj_ComplexS1D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_ComplexS1D(np.ndarray data)
        naga_host_obj_ComplexS1D(np.ndarray data, int nbStream)
        naga_host_obj_ComplexS1D(np.ndarray data, str mallocType)
        naga_host_obj_ComplexS1D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_ComplexS1D(naga_host_obj_ComplexS1D obj)
        naga_host_obj_ComplexS1D(naga_host_obj_ComplexS1D obj, int nbStream)
        naga_host_obj_ComplexS1D(naga_host_obj_ComplexS1D obj, str mallocType)
        naga_host_obj_ComplexS1D(naga_host_obj_ComplexS1D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=1
        dims must be a np array of shape (1)
        """

        cdef long sh[1+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(1):
                        sh[i+1]=dims[i]
                    sh[0]=1
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=1
                for i in range(1):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((1),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(1):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_ComplexS1D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_ComplexS1D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuFloatComplex]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_ComplexS1D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_ComplexS1D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuFloatComplex]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((1),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=1,dtype=np.complex64_t] data

        cdims=self.c_h.getDims()

        for i in range(1):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.complex64,order="F")

        self.c_h.fill_into(<cuFloatComplex*>data.data)
        return data


    def setData(self,np.ndarray[ndim=1,dtype=np.complex64_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<cuFloatComplex*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_ComplexS2D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_ComplexS2D obj=None,
            np.ndarray[ndim=2,dtype=np.complex64_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_ComplexS2D constructor.

    constructors available:
        naga_host_obj_ComplexS2D(np.ndarray dims)
        naga_host_obj_ComplexS2D(np.ndarray dims, int nbStreams)
        naga_host_obj_ComplexS2D(np.ndarray dims, str mallocType)
        naga_host_obj_ComplexS2D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_ComplexS2D(np.ndarray data)
        naga_host_obj_ComplexS2D(np.ndarray data, int nbStream)
        naga_host_obj_ComplexS2D(np.ndarray data, str mallocType)
        naga_host_obj_ComplexS2D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_ComplexS2D(naga_host_obj_ComplexS2D obj)
        naga_host_obj_ComplexS2D(naga_host_obj_ComplexS2D obj, int nbStream)
        naga_host_obj_ComplexS2D(naga_host_obj_ComplexS2D obj, str mallocType)
        naga_host_obj_ComplexS2D(naga_host_obj_ComplexS2D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=2
        dims must be a np array of shape (2)
        """

        cdef long sh[2+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(2):
                        sh[i+1]=dims[i]
                    sh[0]=2
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=2
                for i in range(2):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((2),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(2):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_ComplexS2D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_ComplexS2D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuFloatComplex]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_ComplexS2D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_ComplexS2D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuFloatComplex]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((2),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=2,dtype=np.complex64_t] data

        cdims=self.c_h.getDims()

        for i in range(2):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.complex64,order="F")

        self.c_h.fill_into(<cuFloatComplex*>data.data)
        return data


    def setData(self,np.ndarray[ndim=2,dtype=np.complex64_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<cuFloatComplex*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_ComplexS3D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_ComplexS3D obj=None,
            np.ndarray[ndim=3,dtype=np.complex64_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_ComplexS3D constructor.

    constructors available:
        naga_host_obj_ComplexS3D(np.ndarray dims)
        naga_host_obj_ComplexS3D(np.ndarray dims, int nbStreams)
        naga_host_obj_ComplexS3D(np.ndarray dims, str mallocType)
        naga_host_obj_ComplexS3D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_ComplexS3D(np.ndarray data)
        naga_host_obj_ComplexS3D(np.ndarray data, int nbStream)
        naga_host_obj_ComplexS3D(np.ndarray data, str mallocType)
        naga_host_obj_ComplexS3D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_ComplexS3D(naga_host_obj_ComplexS3D obj)
        naga_host_obj_ComplexS3D(naga_host_obj_ComplexS3D obj, int nbStream)
        naga_host_obj_ComplexS3D(naga_host_obj_ComplexS3D obj, str mallocType)
        naga_host_obj_ComplexS3D(naga_host_obj_ComplexS3D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=3
        dims must be a np array of shape (3)
        """

        cdef long sh[3+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(3):
                        sh[i+1]=dims[i]
                    sh[0]=3
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=3
                for i in range(3):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((3),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(3):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_ComplexS3D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_ComplexS3D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuFloatComplex]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_ComplexS3D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_ComplexS3D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuFloatComplex]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((3),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=3,dtype=np.complex64_t] data

        cdims=self.c_h.getDims()

        for i in range(3):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.complex64,order="F")

        self.c_h.fill_into(<cuFloatComplex*>data.data)
        return data


    def setData(self,np.ndarray[ndim=3,dtype=np.complex64_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<cuFloatComplex*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_ComplexS4D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_ComplexS4D obj=None,
            np.ndarray[ndim=4,dtype=np.complex64_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_ComplexS4D constructor.

    constructors available:
        naga_host_obj_ComplexS4D(np.ndarray dims)
        naga_host_obj_ComplexS4D(np.ndarray dims, int nbStreams)
        naga_host_obj_ComplexS4D(np.ndarray dims, str mallocType)
        naga_host_obj_ComplexS4D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_ComplexS4D(np.ndarray data)
        naga_host_obj_ComplexS4D(np.ndarray data, int nbStream)
        naga_host_obj_ComplexS4D(np.ndarray data, str mallocType)
        naga_host_obj_ComplexS4D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_ComplexS4D(naga_host_obj_ComplexS4D obj)
        naga_host_obj_ComplexS4D(naga_host_obj_ComplexS4D obj, int nbStream)
        naga_host_obj_ComplexS4D(naga_host_obj_ComplexS4D obj, str mallocType)
        naga_host_obj_ComplexS4D(naga_host_obj_ComplexS4D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=4
        dims must be a np array of shape (4)
        """

        cdef long sh[4+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(4):
                        sh[i+1]=dims[i]
                    sh[0]=4
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=4
                for i in range(4):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data)
                    else:
                        self.c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[cuFloatComplex](sh,
                                                            <cuFloatComplex*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[cuFloatComplex](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((4),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(4):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_ComplexS4D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_ComplexS4D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuFloatComplex]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_ComplexS4D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_ComplexS4D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuFloatComplex]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((4),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=4,dtype=np.complex64_t] data

        cdims=self.c_h.getDims()

        for i in range(4):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.complex64,order="F")

        self.c_h.fill_into(<cuFloatComplex*>data.data)
        return data


    def setData(self,np.ndarray[ndim=4,dtype=np.complex64_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<cuFloatComplex*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Double1D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Double1D obj=None,
            np.ndarray[ndim=1,dtype=np.float64_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Double1D constructor.

    constructors available:
        naga_host_obj_Double1D(np.ndarray dims)
        naga_host_obj_Double1D(np.ndarray dims, int nbStreams)
        naga_host_obj_Double1D(np.ndarray dims, str mallocType)
        naga_host_obj_Double1D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Double1D(np.ndarray data)
        naga_host_obj_Double1D(np.ndarray data, int nbStream)
        naga_host_obj_Double1D(np.ndarray data, str mallocType)
        naga_host_obj_Double1D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Double1D(naga_host_obj_Double1D obj)
        naga_host_obj_Double1D(naga_host_obj_Double1D obj, int nbStream)
        naga_host_obj_Double1D(naga_host_obj_Double1D obj, str mallocType)
        naga_host_obj_Double1D(naga_host_obj_Double1D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=1
        dims must be a np array of shape (1)
        """

        cdef long sh[1+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(1):
                        sh[i+1]=dims[i]
                    sh[0]=1
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[double](sh)
                    else:
                        self.c_h=new carma_host_obj[double](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[double](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[double](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=1
                for i in range(1):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data)
                    else:
                        self.c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[double](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[double](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[double](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[double](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((1),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(1):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Double1D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Double1D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[double]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Double1D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Double1D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[double]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((1),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=1,dtype=np.float64_t] data

        cdims=self.c_h.getDims()

        for i in range(1):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.float64,order="F")

        self.c_h.fill_into(<double*>data.data)
        return data


    def setData(self,np.ndarray[ndim=1,dtype=np.float64_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<double*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Double2D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Double2D obj=None,
            np.ndarray[ndim=2,dtype=np.float64_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Double2D constructor.

    constructors available:
        naga_host_obj_Double2D(np.ndarray dims)
        naga_host_obj_Double2D(np.ndarray dims, int nbStreams)
        naga_host_obj_Double2D(np.ndarray dims, str mallocType)
        naga_host_obj_Double2D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Double2D(np.ndarray data)
        naga_host_obj_Double2D(np.ndarray data, int nbStream)
        naga_host_obj_Double2D(np.ndarray data, str mallocType)
        naga_host_obj_Double2D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Double2D(naga_host_obj_Double2D obj)
        naga_host_obj_Double2D(naga_host_obj_Double2D obj, int nbStream)
        naga_host_obj_Double2D(naga_host_obj_Double2D obj, str mallocType)
        naga_host_obj_Double2D(naga_host_obj_Double2D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=2
        dims must be a np array of shape (2)
        """

        cdef long sh[2+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(2):
                        sh[i+1]=dims[i]
                    sh[0]=2
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[double](sh)
                    else:
                        self.c_h=new carma_host_obj[double](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[double](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[double](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=2
                for i in range(2):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data)
                    else:
                        self.c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[double](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[double](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[double](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[double](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((2),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(2):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Double2D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Double2D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[double]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Double2D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Double2D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[double]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((2),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=2,dtype=np.float64_t] data

        cdims=self.c_h.getDims()

        for i in range(2):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.float64,order="F")

        self.c_h.fill_into(<double*>data.data)
        return data


    def setData(self,np.ndarray[ndim=2,dtype=np.float64_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<double*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Double3D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Double3D obj=None,
            np.ndarray[ndim=3,dtype=np.float64_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Double3D constructor.

    constructors available:
        naga_host_obj_Double3D(np.ndarray dims)
        naga_host_obj_Double3D(np.ndarray dims, int nbStreams)
        naga_host_obj_Double3D(np.ndarray dims, str mallocType)
        naga_host_obj_Double3D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Double3D(np.ndarray data)
        naga_host_obj_Double3D(np.ndarray data, int nbStream)
        naga_host_obj_Double3D(np.ndarray data, str mallocType)
        naga_host_obj_Double3D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Double3D(naga_host_obj_Double3D obj)
        naga_host_obj_Double3D(naga_host_obj_Double3D obj, int nbStream)
        naga_host_obj_Double3D(naga_host_obj_Double3D obj, str mallocType)
        naga_host_obj_Double3D(naga_host_obj_Double3D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=3
        dims must be a np array of shape (3)
        """

        cdef long sh[3+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(3):
                        sh[i+1]=dims[i]
                    sh[0]=3
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[double](sh)
                    else:
                        self.c_h=new carma_host_obj[double](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[double](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[double](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=3
                for i in range(3):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data)
                    else:
                        self.c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[double](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[double](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[double](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[double](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((3),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(3):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Double3D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Double3D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[double]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Double3D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Double3D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[double]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((3),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=3,dtype=np.float64_t] data

        cdims=self.c_h.getDims()

        for i in range(3):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.float64,order="F")

        self.c_h.fill_into(<double*>data.data)
        return data


    def setData(self,np.ndarray[ndim=3,dtype=np.float64_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<double*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_Double4D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_Double4D obj=None,
            np.ndarray[ndim=4,dtype=np.float64_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_Double4D constructor.

    constructors available:
        naga_host_obj_Double4D(np.ndarray dims)
        naga_host_obj_Double4D(np.ndarray dims, int nbStreams)
        naga_host_obj_Double4D(np.ndarray dims, str mallocType)
        naga_host_obj_Double4D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_Double4D(np.ndarray data)
        naga_host_obj_Double4D(np.ndarray data, int nbStream)
        naga_host_obj_Double4D(np.ndarray data, str mallocType)
        naga_host_obj_Double4D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_Double4D(naga_host_obj_Double4D obj)
        naga_host_obj_Double4D(naga_host_obj_Double4D obj, int nbStream)
        naga_host_obj_Double4D(naga_host_obj_Double4D obj, str mallocType)
        naga_host_obj_Double4D(naga_host_obj_Double4D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=4
        dims must be a np array of shape (4)
        """

        cdef long sh[4+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(4):
                        sh[i+1]=dims[i]
                    sh[0]=4
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[double](sh)
                    else:
                        self.c_h=new carma_host_obj[double](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[double](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[double](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=4
                for i in range(4):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data)
                    else:
                        self.c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[double](sh,
                                                            <double*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[double](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[double](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[double](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[double](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((4),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(4):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_Double4D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_Double4D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[double]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_Double4D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_Double4D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[double]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((4),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=4,dtype=np.float64_t] data

        cdims=self.c_h.getDims()

        for i in range(4):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.float64,order="F")

        self.c_h.fill_into(<double*>data.data)
        return data


    def setData(self,np.ndarray[ndim=4,dtype=np.float64_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<double*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_ComplexD1D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_ComplexD1D obj=None,
            np.ndarray[ndim=1,dtype=np.complex128_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_ComplexD1D constructor.

    constructors available:
        naga_host_obj_ComplexD1D(np.ndarray dims)
        naga_host_obj_ComplexD1D(np.ndarray dims, int nbStreams)
        naga_host_obj_ComplexD1D(np.ndarray dims, str mallocType)
        naga_host_obj_ComplexD1D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_ComplexD1D(np.ndarray data)
        naga_host_obj_ComplexD1D(np.ndarray data, int nbStream)
        naga_host_obj_ComplexD1D(np.ndarray data, str mallocType)
        naga_host_obj_ComplexD1D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_ComplexD1D(naga_host_obj_ComplexD1D obj)
        naga_host_obj_ComplexD1D(naga_host_obj_ComplexD1D obj, int nbStream)
        naga_host_obj_ComplexD1D(naga_host_obj_ComplexD1D obj, str mallocType)
        naga_host_obj_ComplexD1D(naga_host_obj_ComplexD1D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=1
        dims must be a np array of shape (1)
        """

        cdef long sh[1+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(1):
                        sh[i+1]=dims[i]
                    sh[0]=1
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=1
                for i in range(1):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((1),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(1):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_ComplexD1D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_ComplexD1D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuDoubleComplex]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_ComplexD1D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_ComplexD1D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuDoubleComplex]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((1),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=1,dtype=np.complex128_t] data

        cdims=self.c_h.getDims()

        for i in range(1):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.complex128,order="F")

        self.c_h.fill_into(<cuDoubleComplex*>data.data)
        return data


    def setData(self,np.ndarray[ndim=1,dtype=np.complex128_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<cuDoubleComplex*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_ComplexD2D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_ComplexD2D obj=None,
            np.ndarray[ndim=2,dtype=np.complex128_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_ComplexD2D constructor.

    constructors available:
        naga_host_obj_ComplexD2D(np.ndarray dims)
        naga_host_obj_ComplexD2D(np.ndarray dims, int nbStreams)
        naga_host_obj_ComplexD2D(np.ndarray dims, str mallocType)
        naga_host_obj_ComplexD2D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_ComplexD2D(np.ndarray data)
        naga_host_obj_ComplexD2D(np.ndarray data, int nbStream)
        naga_host_obj_ComplexD2D(np.ndarray data, str mallocType)
        naga_host_obj_ComplexD2D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_ComplexD2D(naga_host_obj_ComplexD2D obj)
        naga_host_obj_ComplexD2D(naga_host_obj_ComplexD2D obj, int nbStream)
        naga_host_obj_ComplexD2D(naga_host_obj_ComplexD2D obj, str mallocType)
        naga_host_obj_ComplexD2D(naga_host_obj_ComplexD2D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=2
        dims must be a np array of shape (2)
        """

        cdef long sh[2+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(2):
                        sh[i+1]=dims[i]
                    sh[0]=2
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=2
                for i in range(2):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((2),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(2):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_ComplexD2D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_ComplexD2D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuDoubleComplex]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_ComplexD2D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_ComplexD2D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuDoubleComplex]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((2),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=2,dtype=np.complex128_t] data

        cdims=self.c_h.getDims()

        for i in range(2):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.complex128,order="F")

        self.c_h.fill_into(<cuDoubleComplex*>data.data)
        return data


    def setData(self,np.ndarray[ndim=2,dtype=np.complex128_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<cuDoubleComplex*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_ComplexD3D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_ComplexD3D obj=None,
            np.ndarray[ndim=3,dtype=np.complex128_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_ComplexD3D constructor.

    constructors available:
        naga_host_obj_ComplexD3D(np.ndarray dims)
        naga_host_obj_ComplexD3D(np.ndarray dims, int nbStreams)
        naga_host_obj_ComplexD3D(np.ndarray dims, str mallocType)
        naga_host_obj_ComplexD3D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_ComplexD3D(np.ndarray data)
        naga_host_obj_ComplexD3D(np.ndarray data, int nbStream)
        naga_host_obj_ComplexD3D(np.ndarray data, str mallocType)
        naga_host_obj_ComplexD3D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_ComplexD3D(naga_host_obj_ComplexD3D obj)
        naga_host_obj_ComplexD3D(naga_host_obj_ComplexD3D obj, int nbStream)
        naga_host_obj_ComplexD3D(naga_host_obj_ComplexD3D obj, str mallocType)
        naga_host_obj_ComplexD3D(naga_host_obj_ComplexD3D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=3
        dims must be a np array of shape (3)
        """

        cdef long sh[3+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(3):
                        sh[i+1]=dims[i]
                    sh[0]=3
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=3
                for i in range(3):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((3),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(3):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_ComplexD3D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_ComplexD3D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuDoubleComplex]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_ComplexD3D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_ComplexD3D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuDoubleComplex]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((3),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=3,dtype=np.complex128_t] data

        cdims=self.c_h.getDims()

        for i in range(3):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.complex128,order="F")

        self.c_h.fill_into(<cuDoubleComplex*>data.data)
        return data


    def setData(self,np.ndarray[ndim=3,dtype=np.complex128_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<cuDoubleComplex*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

cdef class naga_host_obj_ComplexD4D:

    def __cinit__(self,
            np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
            naga_host_obj_ComplexD4D obj=None,
            np.ndarray[ndim=4,dtype=np.complex128_t] data=None,
            str mallocType= "none",
            int nbStreams=-1
            ):

        """naga_host_obj_ComplexD4D constructor.

    constructors available:
        naga_host_obj_ComplexD4D(np.ndarray dims)
        naga_host_obj_ComplexD4D(np.ndarray dims, int nbStreams)
        naga_host_obj_ComplexD4D(np.ndarray dims, str mallocType)
        naga_host_obj_ComplexD4D(np.ndarray dims, str mallocType, int nbStreams)
        naga_host_obj_ComplexD4D(np.ndarray data)
        naga_host_obj_ComplexD4D(np.ndarray data, int nbStream)
        naga_host_obj_ComplexD4D(np.ndarray data, str mallocType)
        naga_host_obj_ComplexD4D(np.ndarray data, str mallocType, int nbStreams)
        naga_host_obj_ComplexD4D(naga_host_obj_ComplexD4D obj)
        naga_host_obj_ComplexD4D(naga_host_obj_ComplexD4D obj, int nbStream)
        naga_host_obj_ComplexD4D(naga_host_obj_ComplexD4D obj, str mallocType)
        naga_host_obj_ComplexD4D(naga_host_obj_ComplexD4D obj, str mallocType, int nbStreams)

        input data must be a np.ndarray with ndim=4
        dims must be a np array of shape (4)
        """

        cdef long sh[4+1]
        cdef int i


        cdef MemAlloc MT = dict_MemAlloc[mallocType]

        if(obj is None):
            if(data is None):
                try:
                    for i in range(4):
                        sh[i+1]=dims[i]
                    sh[0]=4
                except ValueError:
                    print("At least: argument 'dims' is required")

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, nbStreams)
                else:
                    if(nbStreams==-1):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, MT)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh, MT, nbStreams)

            else:
                if(not data.flags.fortran):
                    data=np.asfortranarray(data)
                sh[0]=4
                for i in range(4):
                     sh[i+1]=data.shape[i]

                if(MT ==-1):
                    if(nbStreams is None):
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data)
                    else:
                        self.c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            nbStreams)
                else:
                    if(nbStreams==-1):
                        self. c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            MT)
                    else:
                        self. c_h=new carma_host_obj[cuDoubleComplex](sh,
                                                            <cuDoubleComplex*>data.data,
                                                            MT,nbStreams)

        elif(dims is None and data is None):
            if(MT==-1):
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h)
                else:
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,nbStreams)
            else:
                if(nbStreams==-1):
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,MT)
                else:
                    self.c_h=new carma_host_obj[cuDoubleComplex](obj.c_h,MT,nbStreams)


        else:
           raise ValueError("----")


        self.data_h=self.c_h.getData()


    def __dealloc__(self):
        del self.c_h

    def get_Dims(self):
        cdef dims= np.ndarray((4),dtype=np.int64)
        cdef const long *cdims
        cdef int i

        cdims=self.c_h.getDims()
        for i in range(4):
            dims[i]=cdims[i+1]
        return dims

    def get_host_obj_ptr(self):
        """Return the pointer to the naga_host_obj."""
        cdef uintptr_t host_obj_ptr=<uintptr_t>self.c_h.getData()
        return host_obj_ptr


    def cpy_obj_from(self, naga_obj_ComplexD4D src):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        src -- naga_obj_ComplexD4D: object to copy from
        """
        cdef uintptr_t src_ptr = <uintptr_t>src.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuDoubleComplex]*>src_ptr,cudaMemcpyDeviceToHost)

    def cpy_obj_into(self,naga_obj_ComplexD4D dest):
        """Copy into naga_host_obj (cpu storage) the data from a naga_obj (gpu storage).

        dest -- naga_obj_ComplexD4D: object to copy into
        """
        cdef uintptr_t dest_ptr = <uintptr_t>dest.getCarma_ptr()
        self.c_h.cpy_obj(<carma_obj[cuDoubleComplex]*>dest_ptr,cudaMemcpyHostToDevice)


    def getData(self):
        cdef int i
        cdef sh=np.zeros((4),dtype=np.int64)
        cdef const long *cdims
        cdef np.ndarray[ndim=4,dtype=np.complex128_t] data

        cdims=self.c_h.getDims()

        for i in range(4):
            sh[i]=cdims[i+1]

        data=np.zeros(sh,dtype=np.complex128,order="F")

        self.c_h.fill_into(<cuDoubleComplex*>data.data)
        return data


    def setData(self,np.ndarray[ndim=4,dtype=np.complex128_t] data):
        if(not data.flags.fortran):
            data=np.asfortranarray(data)
        self.c_h.fill_from(<cuDoubleComplex*>data.data)

    def getNbElem(self):
        return self.c_h.getNbElem()

