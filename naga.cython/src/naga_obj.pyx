

import numpy as np
cimport numpy as np
#np.import_array()

assert sizeof(int) == sizeof(np.int32_t)
assert sizeof(float) == sizeof(np.float32_t)
assert sizeof(double) == sizeof(np.float64_t)


"""
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

        carma_obj(carma_obj[T] *obj)
        carma_obj(carma_context *context, const long *dims)
        carma_obj(carma_context *context,carma_obj[T] *obj)
        carma_obj(carma_context *context,long *dims, T *data)

        T* getData()
        T* getDataAt(int index)
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
        int init_prng(int seed)
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
"""
#Dictionary fft available:
# Dictionary contain fft input type
# key:   fft output type
# value: cufftType
cdef dict_fft_Float1D={naga_obj_ComplexS1D:(CUFFT_R2C)}
cdef dict_fft_Float2D={naga_obj_ComplexS2D:(CUFFT_R2C)}
cdef dict_fft_Float3D={naga_obj_ComplexS3D:(CUFFT_R2C)}
cdef dict_fft_Float4D={naga_obj_ComplexS4D:(CUFFT_R2C)}
cdef dict_fft_Double1D={naga_obj_ComplexD1D:(CUFFT_D2Z)}
cdef dict_fft_Double2D={naga_obj_ComplexD2D:(CUFFT_D2Z)}
cdef dict_fft_Double3D={naga_obj_ComplexD3D:(CUFFT_D2Z)}
cdef dict_fft_Double4D={naga_obj_ComplexD4D:(CUFFT_D2Z)}
cdef dict_fft_ComplexS1D={
    naga_obj_Float1D:(CUFFT_C2R),
    naga_obj_ComplexS1D:(CUFFT_C2C)
}
cdef dict_fft_ComplexS2D={
    naga_obj_Float2D:(CUFFT_C2R),
    naga_obj_ComplexS2D:(CUFFT_C2C)
}
cdef dict_fft_ComplexS3D={
    naga_obj_Float3D:(CUFFT_C2R),
    naga_obj_ComplexS3D:(CUFFT_C2C)
}
cdef dict_fft_ComplexS4D={
    naga_obj_Float4D:(CUFFT_C2R),
    naga_obj_ComplexS4D:(CUFFT_C2C)
}
cdef dict_fft_ComplexD1D={
    naga_obj_Double1D:(CUFFT_Z2D),
    naga_obj_ComplexD1D:(CUFFT_Z2Z)
}
cdef dict_fft_ComplexD2D={
    naga_obj_Double2D:(CUFFT_Z2D),
    naga_obj_ComplexD2D:(CUFFT_Z2Z)
}
cdef dict_fft_ComplexD3D={
    naga_obj_Double3D:(CUFFT_Z2D),
    naga_obj_ComplexD3D:(CUFFT_Z2Z)
}
cdef dict_fft_ComplexD4D={
    naga_obj_Double4D:(CUFFT_Z2D),
    naga_obj_ComplexD4D:(CUFFT_Z2Z)
}




##########################################################
#  P-Class `Int1D`
#########################################################


cdef class naga_obj_Int1D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Int1D obj=None,
                  np.ndarray[ndim=1, dtype=np.int32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Int1D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=1
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[int](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[int](ctxt.c,obj.c_o)
            return

        cdef long cdims[1+1]
        cdims[0]=1
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=1]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=1 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 1)")
            for i in range(1):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[int](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[int](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(1):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[int](ctxt.c,cdims,<int*>data.data)
        else:
            self.c_o=new carma_obj[int](ctxt.c,cdims,<int*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((1+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(1+1):
            dims[i]=cdims[i]
        return dims[1:1+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=1, dtype=np.int32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=1, dtype=np.int32_t] data):
        data -- np.int32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<int*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<int*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.int32) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] res=np.zeros(sh,dtype=np.int32)
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.int32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.int32) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] res=np.zeros(sh, dtype=np.int32)
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.int32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Int1D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Int1D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Int1D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Int1D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



##########################################################
#  P-Class `Int2D`
#########################################################


cdef class naga_obj_Int2D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Int2D obj=None,
                  np.ndarray[ndim=2, dtype=np.int32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Int2D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=2
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[int](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[int](ctxt.c,obj.c_o)
            return

        cdef long cdims[2+1]
        cdims[0]=2
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=2]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=2 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 2)")
            for i in range(2):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[int](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[int](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(2):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[int](ctxt.c,cdims,<int*>data.data)
        else:
            self.c_o=new carma_obj[int](ctxt.c,cdims,<int*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((2+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(2+1):
            dims[i]=cdims[i]
        return dims[1:2+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=2, dtype=np.int32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=2, dtype=np.int32_t] data):
        data -- np.int32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<int*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<int*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.int32) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.int32_t, ndim=2] res=np.zeros(sh,dtype=np.int32)
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.int32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.int32) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.int32_t, ndim=2] res=np.zeros(sh, dtype=np.int32)
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.int32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Int2D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Int2D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Int2D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Int2D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



##########################################################
#  P-Class `Int3D`
#########################################################


cdef class naga_obj_Int3D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Int3D obj=None,
                  np.ndarray[ndim=3, dtype=np.int32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Int3D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=3
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[int](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[int](ctxt.c,obj.c_o)
            return

        cdef long cdims[3+1]
        cdims[0]=3
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=3]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=3 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 3)")
            for i in range(3):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[int](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[int](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(3):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[int](ctxt.c,cdims,<int*>data.data)
        else:
            self.c_o=new carma_obj[int](ctxt.c,cdims,<int*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((3+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(3+1):
            dims[i]=cdims[i]
        return dims[1:3+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=3, dtype=np.int32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=3, dtype=np.int32_t] data):
        data -- np.int32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<int*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<int*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.int32) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.int32_t, ndim=3] res=np.zeros(sh,dtype=np.int32)
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.int32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.int32) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.int32_t, ndim=3] res=np.zeros(sh, dtype=np.int32)
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.int32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Int3D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Int3D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Int3D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Int3D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



##########################################################
#  P-Class `Int4D`
#########################################################


cdef class naga_obj_Int4D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Int4D obj=None,
                  np.ndarray[ndim=4, dtype=np.int32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Int4D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=4
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[int](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[int](ctxt.c,obj.c_o)
            return

        cdef long cdims[4+1]
        cdims[0]=4
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=4]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=4 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 4)")
            for i in range(4):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[int](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[int](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(4):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[int](ctxt.c,cdims,<int*>data.data)
        else:
            self.c_o=new carma_obj[int](ctxt.c,cdims,<int*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((4+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(4+1):
            dims[i]=cdims[i]
        return dims[1:4+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=4, dtype=np.int32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=4, dtype=np.int32_t] data):
        data -- np.int32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<int*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<int*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.int32) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.int32_t, ndim=4] res=np.zeros(sh,dtype=np.int32)
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.int32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.int32) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.int32_t, ndim=4] res=np.zeros(sh, dtype=np.int32)
        cdef np.ndarray[ dtype=np.int32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.int32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Int4D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Int4D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Int4D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Int4D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



##########################################################
#  P-Class `UInt1D`
#########################################################


cdef class naga_obj_UInt1D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_UInt1D obj=None,
                  np.ndarray[ndim=1, dtype=np.uint32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_UInt1D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=1
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[unsigned int](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[unsigned int](ctxt.c,obj.c_o)
            return

        cdef long cdims[1+1]
        cdims[0]=1
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=1]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=1 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 1)")
            for i in range(1):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[unsigned int](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[unsigned int](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(1):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[unsigned int](ctxt.c,cdims,<unsigned int*>data.data)
        else:
            self.c_o=new carma_obj[unsigned int](ctxt.c,cdims,<unsigned int*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((1+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(1+1):
            dims[i]=cdims[i]
        return dims[1:1+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=1, dtype=np.uint32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=1, dtype=np.uint32_t] data):
        data -- np.uint32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<unsigned int*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<unsigned int*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.uint32) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] res=np.zeros(sh,dtype=np.uint32)
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.uint32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<unsigned int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.uint32) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] res=np.zeros(sh, dtype=np.uint32)
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.uint32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<unsigned int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_UInt1D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_UInt1D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_UInt1D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_UInt1D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



##########################################################
#  P-Class `UInt2D`
#########################################################


cdef class naga_obj_UInt2D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_UInt2D obj=None,
                  np.ndarray[ndim=2, dtype=np.uint32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_UInt2D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=2
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[unsigned int](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[unsigned int](ctxt.c,obj.c_o)
            return

        cdef long cdims[2+1]
        cdims[0]=2
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=2]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=2 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 2)")
            for i in range(2):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[unsigned int](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[unsigned int](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(2):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[unsigned int](ctxt.c,cdims,<unsigned int*>data.data)
        else:
            self.c_o=new carma_obj[unsigned int](ctxt.c,cdims,<unsigned int*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((2+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(2+1):
            dims[i]=cdims[i]
        return dims[1:2+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=2, dtype=np.uint32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=2, dtype=np.uint32_t] data):
        data -- np.uint32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<unsigned int*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<unsigned int*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.uint32) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.uint32_t, ndim=2] res=np.zeros(sh,dtype=np.uint32)
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.uint32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<unsigned int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.uint32) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.uint32_t, ndim=2] res=np.zeros(sh, dtype=np.uint32)
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.uint32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<unsigned int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_UInt2D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_UInt2D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_UInt2D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_UInt2D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



##########################################################
#  P-Class `UInt3D`
#########################################################


cdef class naga_obj_UInt3D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_UInt3D obj=None,
                  np.ndarray[ndim=3, dtype=np.uint32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_UInt3D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=3
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[unsigned int](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[unsigned int](ctxt.c,obj.c_o)
            return

        cdef long cdims[3+1]
        cdims[0]=3
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=3]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=3 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 3)")
            for i in range(3):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[unsigned int](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[unsigned int](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(3):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[unsigned int](ctxt.c,cdims,<unsigned int*>data.data)
        else:
            self.c_o=new carma_obj[unsigned int](ctxt.c,cdims,<unsigned int*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((3+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(3+1):
            dims[i]=cdims[i]
        return dims[1:3+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=3, dtype=np.uint32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=3, dtype=np.uint32_t] data):
        data -- np.uint32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<unsigned int*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<unsigned int*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.uint32) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.uint32_t, ndim=3] res=np.zeros(sh,dtype=np.uint32)
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.uint32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<unsigned int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.uint32) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.uint32_t, ndim=3] res=np.zeros(sh, dtype=np.uint32)
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.uint32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<unsigned int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_UInt3D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_UInt3D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_UInt3D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_UInt3D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



##########################################################
#  P-Class `UInt4D`
#########################################################


cdef class naga_obj_UInt4D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_UInt4D obj=None,
                  np.ndarray[ndim=4, dtype=np.uint32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_UInt4D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=4
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[unsigned int](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[unsigned int](ctxt.c,obj.c_o)
            return

        cdef long cdims[4+1]
        cdims[0]=4
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=4]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=4 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 4)")
            for i in range(4):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[unsigned int](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[unsigned int](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(4):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[unsigned int](ctxt.c,cdims,<unsigned int*>data.data)
        else:
            self.c_o=new carma_obj[unsigned int](ctxt.c,cdims,<unsigned int*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((4+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(4+1):
            dims[i]=cdims[i]
        return dims[1:4+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=4, dtype=np.uint32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=4, dtype=np.uint32_t] data):
        data -- np.uint32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<unsigned int*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<unsigned int*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.uint32) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.uint32_t, ndim=4] res=np.zeros(sh,dtype=np.uint32)
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.uint32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<unsigned int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.uint32) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.uint32_t, ndim=4] res=np.zeros(sh, dtype=np.uint32)
        cdef np.ndarray[ dtype=np.uint32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.uint32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<unsigned int*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_UInt4D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_UInt4D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_UInt4D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_UInt4D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



##########################################################
#  P-Class `Float1D`
#########################################################


cdef class naga_obj_Float1D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Float1D obj=None,
                  np.ndarray[ndim=1, dtype=np.float32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Float1D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=1
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[float](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[float](ctxt.c,obj.c_o)
            return

        cdef long cdims[1+1]
        cdims[0]=1
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=1]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=1 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 1)")
            for i in range(1):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[float](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[float](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(1):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[float](ctxt.c,cdims,<float*>data.data)
        else:
            self.c_o=new carma_obj[float](ctxt.c,cdims,<float*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((1+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(1+1):
            dims[i]=cdims[i]
        return dims[1:1+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=1, dtype=np.float32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=1, dtype=np.float32_t] data):
        data -- np.float32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<float*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<float*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.float32) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] res=np.zeros(sh,dtype=np.float32)
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<float*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.float32) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] res=np.zeros(sh, dtype=np.float32)
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<float*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Float1D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Float1D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Float1D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Float1D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_Float1D src):
        """Cublas dot

        src -- naga_obj_Float1D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.float32_t alpha):
        """Cublas scal

        alpha -- np.float32: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef float a
        a=alpha

        self.c_o.scale(a,1)


    def swap(self, naga_obj_Float1D src):
        """Cublas swap

        src -- naga_obj_Float1D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_Float1D src):
        """Cublas copy

        src -- naga_obj_Float1D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)


    def init_prng(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.init_prng(seed)


    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng_host(seed)
        self.c_o.prng_host(b'U')


    def montagn(self, float ta_jolie_maman,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.prng_montagn(ta_jolie_maman)



    def axpy( self, np.float32_t alpha, naga_obj_Float1D dest=None):
        """cublas axpy

        dest -- naga_obj_Float1D
        alpha-- np.float32
        beta -- np.float32
        Return dest=alpha*self +dest
        """
        if(dest==None):
            c=self.getContext()
            dest=naga_obj_Float1D(ctxt=self.getContext(),dims=self.get_Dims())
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("axpy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef float a
        a=alpha

        dest.c_o.axpy(a,self.c_o,1,1)
        return dest


    def ger(self, naga_obj_Float1D Y, np.float32_t alpha=1.0,
            naga_obj_Float2D A=None):
        """Cublas ger

        Y -- naga_obj_Float1D
        alpha -- np.float32 (default = 1)
        A -- naga_obj_Float2D (default = None)

        Return A=alpha*self*t(y)+A
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimX=self.get_Dims()
        cdef dimY=Y.get_Dims()
        if(A is None):
            dimA=np.ndarray(2,dtype=np.int64)
            dimA[0]=dimY[0]
            dimA[1]=dimX[0]
            A=naga_obj_Float2D(self.getContext(),dims=dimA)
        if(device !=Y.getDevice() or
           device !=A.getDevice()):
            raise ValueError(" ger only on the same device")
        cdef float a
        a=alpha

        A.c_o.ger(a,Y.c_o,1,self.c_o,1,dimY[0])

        return A



    def nrm2(self):
        """Cublas nrm2. Return the Euclidean norm """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.nrm2(1)


    def asum(self):
        """Cublas asum. Return the sum of the absolute values of the data's elements """
        self.context.set_activeDeviceForCpy(self.getDevice())
        return self.c_o.asum(1)


    def sum(self):
        """Return the sum of the data's elements """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.sum()



    def fft(self, naga_obj_ComplexS1D dest=None, int direction=1):
        cdef plan
        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd

        cdef long dims_src[1+1]
        dims_src[0]=1
        for i in range(1):
            dims_src[i+1]=ds[i]

        if(dest is None):
            dd=self.get_Dims()
            dd[1-1]=dd[1-1]/2+1
            dest=naga_obj_ComplexS1D(self.getContext(),dims=dd)

        try:
            plan=dict_fft_Float1D[dest.__class__]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return

        output = dest.getData_ptr()
        carma_initfft[float,cuFloatComplex](dims_src,self.c_o.getPlan(),plan)
        carma_fft[float,cuFloatComplex](self.c_o.getData(), <cuFloatComplex*><void*>output, direction,handle[0])

        return dest


##########################################################
#  P-Class `Float2D`
#########################################################


cdef class naga_obj_Float2D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Float2D obj=None,
                  np.ndarray[ndim=2, dtype=np.float32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Float2D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=2
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[float](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[float](ctxt.c,obj.c_o)
            return

        cdef long cdims[2+1]
        cdims[0]=2
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=2]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=2 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 2)")
            for i in range(2):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[float](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[float](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(2):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[float](ctxt.c,cdims,<float*>data.data)
        else:
            self.c_o=new carma_obj[float](ctxt.c,cdims,<float*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((2+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(2+1):
            dims[i]=cdims[i]
        return dims[1:2+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=2, dtype=np.float32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=2, dtype=np.float32_t] data):
        data -- np.float32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<float*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<float*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.float32) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float32_t, ndim=2] res=np.zeros(sh,dtype=np.float32)
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<float*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.float32) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float32_t, ndim=2] res=np.zeros(sh, dtype=np.float32)
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<float*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Float2D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Float2D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Float2D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Float2D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_Float2D src):
        """Cublas dot

        src -- naga_obj_Float2D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.float32_t alpha):
        """Cublas scal

        alpha -- np.float32: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef float a
        a=alpha

        self.c_o.scale(a,1)


    def swap(self, naga_obj_Float2D src):
        """Cublas swap

        src -- naga_obj_Float2D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_Float2D src):
        """Cublas copy

        src -- naga_obj_Float2D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)


    def init_prng(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.init_prng(seed)


    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng_host(seed)
        self.c_o.prng_host(b'U')


    def montagn(self, float ta_jolie_maman,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.prng_montagn(ta_jolie_maman)



    def gemv(self, naga_obj_Float1D Vx, np.float32_t alpha=1,
             naga_obj_Float1D Vy=None, np.float32_t beta=0 ):
        """Cublas gemv

        Vx -- naga_obj_Float1D
        alpha -- np.float32 (default = 1)
        Vy -- naga_obj_Float1D (default = None)
        beta -- np.float32 (default = 0)
        Return  Vy=alpha*self*Vx+beta*Vy
        """
        device=self.getDevice()
        self.context.set_activeDeviceForCpy(device)
        cdef dimM=self.get_Dims()
        if(Vy is None):
            Vy=naga_obj_Float1D(ctxt=self.getContext(),dims=dimM[0:1])

        if( (self.getDevice()!=Vx.getDevice()) or
             self.getDevice()!= Vy.getDevice()):
            raise ValueError("gemv only on the same device")
        cdef float a
        cdef float b
        a=alpha
        b=beta


        #Vy.c_o.gemv(b'N',a, self.c_o,dimM[1],Vx.c_o,1,b,1)
        Vy.c_o.gemv(b'N',a, self.c_o,dimM[0],Vx.c_o,1,b,1)

        return Vy


    def symv(self,naga_obj_Float1D Vx, np.float32_t alpha=1,
             naga_obj_Float1D Vy=None , np.float32_t beta=0):
        """Cublas symv

        Vx -- naga_obj_Float1D
        alpha -- np.float32 (default = 1)
        Vy -- naga_obj_Float1D (default = None)
        beta -- np.float32 (default = 0)
        Return Vy=alpha*self*Vx+beta*Vy
        """
        cdef device= self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimM=self.get_Dims()
        if(Vy is None):
            #dimY=np.array([1,dimM[1]],dtype=np.int64)
            #Vy=naga_obj_Float2D(ctxt=self.getContext(),dims=dimY)
            Vy=naga_obj_Float2D(ctxt=self.getContext(),dims=dimM)
        if(device !=Vx.getDevice() or
           device !=Vy.getDevice() ):
            raise ValueError("symv only on the same device")

        cdef float a
        cdef float b
        a=alpha
        b=beta

        Vy.c_o.symv(CUBLAS_FILL_MODE_LOWER, a, self.c_o, dimM[0],
                    Vx.c_o, 1, b, 1)

        return Vy


    def gemm(self,naga_obj_Float2D B, str opA='n', str opB='n',
             np.float32_t alpha=1, naga_obj_Float2D C=None, np.float32_t beta=0):
        """Cublas gemm

        B -- naga_obj_Float2D
        opA -- char (default = 'n')
        opB -- char (default = 'n')
        alpha -- np.float32 (default = 1)
        C -- naga_obj_Float2D (default = None)
        beta  -- np.float32 (default = 0)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        Return C=alpha opA(self)*opB(B)+beta*C
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            dimC[0]= dimA[0] if opA=='n' else dimA[1]
            dimC[1]= dimB[1] if opB=='n' else dimB[0]
            C=naga_obj_Float2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device!=B.getDevice() or
           device!=C.getDevice() ):
            raise ValueError("gemm only on the same device")

        cdef float a
        cdef float b
        a=alpha
        b=beta


        C.c_o.gemm(ord(opA),ord(opB),a,self.c_o,dimA[0],B.c_o,dimB[0],b,dimC[0])

        return C


    def symm(self, naga_obj_Float2D B, str side='l', np.float32_t alpha=1,
             naga_obj_Float2D C=None, np.float32_t beta=0):
        """Cublas symm

        B -- naga_obj_Float2D
        side -- char (default = 'l')
        alpha -- np.float32 (default =1)
        C -- naga_obj_Float2D (default = None)
        beta -- np.float32 (default =0)

        return alpha*A*B+beta*C     if side='l'
               alpha*B*A+beta*C     otherwise
        """
        prodSide=CUBLAS_SIDE_LEFT if side=='l' else CUBLAS_SIDE_RIGHT

        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(side=='l'):
                dimC[0]= dimA[0]
                dimC[1]= dimB[1]
            else:
                dimC[0]= dimA[1]
                dimC[1]= dimB[0]

            C=naga_obj_Float2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("symm only on the same device")

        cdef float a
        cdef float b
        a=alpha
        b=beta

        C.c_o.symm(prodSide,CUBLAS_FILL_MODE_LOWER, a,self.c_o,dimA[0],
                    B.c_o,dimB[0],b,dimC[0])
        return C


    def syrk(self, str opA='n', np.float32_t alpha=1, naga_obj_Float2D C=None,
             np.float32_t beta=0):
        """Cublas syrk

        opA -- char (default = 'n')
        alpha -- np.float32 (default = 1)
        C -- naga_obj_Float2D (default = None)
        beta -- np.float32 (default = 0)

        opA: transposition on matrix self
        'n': no transposition
        't':transpose matrix
        Return alpha*opA(self)*opA(self)T+beta*C
        """

        cdef device = self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[0]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[1]

            C=naga_obj_Float2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != C.getDevice()):
            raise ValueError("syrk only on the same device")

        cdef float a
        cdef float b
        a=alpha
        b=beta

        C.c_o.syrk(CUBLAS_FILL_MODE_UPPER, ord(opA),a, self.c_o, dimA[0], b,
                    dimC[0])
        return C


    def syrkx(self, naga_obj_Float2D B, str opA='n',
              np.float32_t alpha=1, naga_obj_Float2D C=None, np.float32_t beta=0):
        """Cublas syrkx

        B -- naga_obj_Float2D
        opA -- char (default = 'n')
        apha -- np.float32 (default = 1)
        C -- naga_obj_Float2D (default = None)
        beta -- np.float32 (default = 0)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        Return alpha*opA(self)*opB(B)T+beta*C
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimB[0]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimB[1]
            C=naga_obj_Float2D(self.getContext(),dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("syrkx only on the same device")

        cdef float a
        cdef float b
        a=alpha
        b=beta

        C.c_o.syrkx(CUBLAS_FILL_MODE_UPPER, ord(opA),a,self.c_o,dimA[0],
                    B.c_o, dimB[0],b,dimC[0])

        return C


    def geam(self, naga_obj_Float2D B, np.float32_t alpha=1, np.float32_t beta=0,
             str opA='n', str opB='n', naga_obj_Float2D C=None):
        """Cublas geam

        B -- naga_obj_Float2D
        alpha -- np.float32 (default = 1)
        beta -- np.float32 (default = 0)
        opA -- char (default = 'n')
        opB -- char (default = 'n')
        C -- naga_obj_Float2D (default = None)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        return C= alpha*opA(self)+beta*opB(B)
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if( C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[1]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[0]
            C=naga_obj_Float2D(self.getContext(),dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("geam only on the same device")


        cdef float a
        cdef float b
        a=alpha
        b=beta

        C.c_o.geam(ord(opA), ord(opB), a, self.c_o, dimA[0],
                   b, B.c_o, dimB[0], dimC[0])
        return C


    def dgmm(self, naga_obj_Float1D X, str sidec='l',
             naga_obj_Float2D C=None):
        """Cublas dgmm

        X -- naga_obj_Float1D
        side -- char (default = 'l')
        C -- naga_obj_Float2D (default = None)

        Return self*diag(X)     if sidec='l'
               diag(X)*self     otherwise
        """

        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimX=X.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(sidec=='l'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[1]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[0]
            C=naga_obj_Float2D(self.getContext(),dims=dimC)
        else:
            dimC=C.getDims()

        if(device != X.getDevice() or
           device != C.getDevice() ):
            raise ValueError("dgmm only in the same device")

        side=CUBLAS_SIDE_LEFT if sidec=='l' else CUBLAS_SIDE_RIGHT

        C.c_o.dgmm(side,self.c_o,dimA[1], X.c_o, 1, dimC[1])
        return C


    def transpose(self, naga_obj_Float2D dest=None):
        """
        """
        cdef device=self.getDevice()
        self.context.set_activeDeviceForCpy(device)
        cdef dimA=self.getDims()
        if(dest is None):
            dimD=np.ndarray(2,dtype=np.int64)
            dimD[0]=dimA[1]
            dimD[1]=dimA[0]
            dest= naga_obj_Float2D(self.getContext,dims=dimD)
        else:
            dimD=dest.getDims()

        if(device != dest.getDevice()):
            raise ValueError("transpose only on the same device")

        dest.c_o.transpose(self.c_o)
        return dest



    def nrm2(self):
        """Cublas nrm2. Return the Euclidean norm """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.nrm2(1)


    def asum(self):
        """Cublas asum. Return the sum of the absolute values of the data's elements """
        self.context.set_activeDeviceForCpy(self.getDevice())
        return self.c_o.asum(1)


    def sum(self):
        """Return the sum of the data's elements """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.sum()



    def fft(self, naga_obj_ComplexS2D dest=None, int direction=1):
        cdef plan
        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd

        cdef long dims_src[2+1]
        dims_src[0]=2
        for i in range(2):
            dims_src[i+1]=ds[i]

        if(dest is None):
            dd=self.get_Dims()
            dd[2-1]=dd[2-1]/2+1
            dest=naga_obj_ComplexS2D(self.getContext(),dims=dd)

        try:
            plan=dict_fft_Float2D[dest.__class__]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return

        output = dest.getData_ptr()
        carma_initfft[float,cuFloatComplex](dims_src,self.c_o.getPlan(),plan)
        carma_fft[float,cuFloatComplex](self.c_o.getData(), <cuFloatComplex*><void*>output, direction,handle[0])

        return dest


##########################################################
#  P-Class `Float3D`
#########################################################


cdef class naga_obj_Float3D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Float3D obj=None,
                  np.ndarray[ndim=3, dtype=np.float32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Float3D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=3
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[float](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[float](ctxt.c,obj.c_o)
            return

        cdef long cdims[3+1]
        cdims[0]=3
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=3]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=3 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 3)")
            for i in range(3):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[float](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[float](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(3):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[float](ctxt.c,cdims,<float*>data.data)
        else:
            self.c_o=new carma_obj[float](ctxt.c,cdims,<float*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((3+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(3+1):
            dims[i]=cdims[i]
        return dims[1:3+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=3, dtype=np.float32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=3, dtype=np.float32_t] data):
        data -- np.float32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<float*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<float*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.float32) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float32_t, ndim=3] res=np.zeros(sh,dtype=np.float32)
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<float*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.float32) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float32_t, ndim=3] res=np.zeros(sh, dtype=np.float32)
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<float*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Float3D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Float3D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Float3D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Float3D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_Float3D src):
        """Cublas dot

        src -- naga_obj_Float3D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.float32_t alpha):
        """Cublas scal

        alpha -- np.float32: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef float a
        a=alpha

        self.c_o.scale(a,1)


    def swap(self, naga_obj_Float3D src):
        """Cublas swap

        src -- naga_obj_Float3D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_Float3D src):
        """Cublas copy

        src -- naga_obj_Float3D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)


    def init_prng(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.init_prng(seed)


    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng_host(seed)
        self.c_o.prng_host(b'U')


    def montagn(self, float ta_jolie_maman,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.prng_montagn(ta_jolie_maman)



    def nrm2(self):
        """Cublas nrm2. Return the Euclidean norm """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.nrm2(1)


    def asum(self):
        """Cublas asum. Return the sum of the absolute values of the data's elements """
        self.context.set_activeDeviceForCpy(self.getDevice())
        return self.c_o.asum(1)


    def sum(self):
        """Return the sum of the data's elements """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.sum()



    def fft(self, naga_obj_ComplexS3D dest=None, int direction=1):
        cdef plan
        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd

        cdef long dims_src[3+1]
        dims_src[0]=3
        for i in range(3):
            dims_src[i+1]=ds[i]

        if(dest is None):
            dd=self.get_Dims()
            dd[3-1]=dd[3-1]/2+1
            dest=naga_obj_ComplexS3D(self.getContext(),dims=dd)

        try:
            plan=dict_fft_Float3D[dest.__class__]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return

        output = dest.getData_ptr()
        carma_initfft[float,cuFloatComplex](dims_src,self.c_o.getPlan(),plan)
        carma_fft[float,cuFloatComplex](self.c_o.getData(), <cuFloatComplex*><void*>output, direction,handle[0])

        return dest


##########################################################
#  P-Class `Float4D`
#########################################################


cdef class naga_obj_Float4D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Float4D obj=None,
                  np.ndarray[ndim=4, dtype=np.float32_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Float4D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=4
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[float](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[float](ctxt.c,obj.c_o)
            return

        cdef long cdims[4+1]
        cdims[0]=4
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=4]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=4 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 4)")
            for i in range(4):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[float](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[float](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(4):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[float](ctxt.c,cdims,<float*>data.data)
        else:
            self.c_o=new carma_obj[float](ctxt.c,cdims,<float*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((4+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(4+1):
            dims[i]=cdims[i]
        return dims[1:4+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=4, dtype=np.float32_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=4, dtype=np.float32_t] data):
        data -- np.float32: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<float*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<float*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.float32) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float32_t, ndim=4] res=np.zeros(sh,dtype=np.float32)
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<float*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.float32) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float32_t, ndim=4] res=np.zeros(sh, dtype=np.float32)
        cdef np.ndarray[ dtype=np.float32_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float32)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<float*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Float4D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Float4D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Float4D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Float4D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_Float4D src):
        """Cublas dot

        src -- naga_obj_Float4D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.float32_t alpha):
        """Cublas scal

        alpha -- np.float32: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef float a
        a=alpha

        self.c_o.scale(a,1)


    def swap(self, naga_obj_Float4D src):
        """Cublas swap

        src -- naga_obj_Float4D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_Float4D src):
        """Cublas copy

        src -- naga_obj_Float4D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)


    def init_prng(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.init_prng(seed)


    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng_host(seed)
        self.c_o.prng_host(b'U')


    def montagn(self, float ta_jolie_maman,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.prng_montagn(ta_jolie_maman)



    def nrm2(self):
        """Cublas nrm2. Return the Euclidean norm """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.nrm2(1)


    def asum(self):
        """Cublas asum. Return the sum of the absolute values of the data's elements """
        self.context.set_activeDeviceForCpy(self.getDevice())
        return self.c_o.asum(1)


    def sum(self):
        """Return the sum of the data's elements """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.sum()



    def fft(self, naga_obj_ComplexS4D dest=None, int direction=1):
        cdef plan
        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd

        cdef long dims_src[4+1]
        dims_src[0]=4
        for i in range(4):
            dims_src[i+1]=ds[i]

        if(dest is None):
            dd=self.get_Dims()
            dd[4-1]=dd[4-1]/2+1
            dest=naga_obj_ComplexS4D(self.getContext(),dims=dd)

        try:
            plan=dict_fft_Float4D[dest.__class__]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return

        output = dest.getData_ptr()
        carma_initfft[float,cuFloatComplex](dims_src,self.c_o.getPlan(),plan)
        carma_fft[float,cuFloatComplex](self.c_o.getData(), <cuFloatComplex*><void*>output, direction,handle[0])

        return dest


##########################################################
#  P-Class `ComplexS1D`
#########################################################


cdef class naga_obj_ComplexS1D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_ComplexS1D obj=None,
                  np.ndarray[ndim=1, dtype=np.complex64_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_ComplexS1D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=1
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[cuFloatComplex](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[cuFloatComplex](ctxt.c,obj.c_o)
            return

        cdef long cdims[1+1]
        cdims[0]=1
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=1]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=1 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 1)")
            for i in range(1):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[cuFloatComplex](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[cuFloatComplex](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(1):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[cuFloatComplex](ctxt.c,cdims,<cuFloatComplex*>data.data)
        else:
            self.c_o=new carma_obj[cuFloatComplex](ctxt.c,cdims,<cuFloatComplex*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((1+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(1+1):
            dims[i]=cdims[i]
        return dims[1:1+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=1, dtype=np.complex64_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=1, dtype=np.complex64_t] data):
        data -- np.complex64: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<cuFloatComplex*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<cuFloatComplex*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.complex64) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] res=np.zeros(sh,dtype=np.complex64)
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<cuFloatComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.complex64) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] res=np.zeros(sh, dtype=np.complex64)
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<cuFloatComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_ComplexS1D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_ComplexS1D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_ComplexS1D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_ComplexS1D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_ComplexS1D src):
        """Cublas dot

        src -- naga_obj_ComplexS1D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.complex64_t alpha):
        """Cublas scal

        alpha -- np.complex64: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuFloatComplex a
        a.x=alpha.x
        a.y=alpha.y

        self.c_o.scale(a,1)


    def swap(self, naga_obj_ComplexS1D src):
        """Cublas swap

        src -- naga_obj_ComplexS1D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_ComplexS1D src):
        """Cublas copy

        src -- naga_obj_ComplexS1D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)



    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng(seed)
        self.c_o.prng(b'U')




    def axpy( self, np.complex64_t alpha, naga_obj_ComplexS1D dest=None):
        """cublas axpy

        dest -- naga_obj_ComplexS1D
        alpha-- np.complex64
        beta -- np.complex64
        Return dest=alpha*self +dest
        """
        if(dest==None):
            c=self.getContext()
            dest=naga_obj_ComplexS1D(ctxt=self.getContext(),dims=self.get_Dims())
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("axpy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuFloatComplex a
        a.x=alpha.x
        a.y=alpha.y

        dest.c_o.axpy(a,self.c_o,1,1)
        return dest


    def ger(self, naga_obj_ComplexS1D Y, np.complex64_t alpha=1.0,
            naga_obj_ComplexS2D A=None):
        """Cublas ger

        Y -- naga_obj_ComplexS1D
        alpha -- np.complex64 (default = 1)
        A -- naga_obj_ComplexS2D (default = None)

        Return A=alpha*self*t(y)+A
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimX=self.get_Dims()
        cdef dimY=Y.get_Dims()
        if(A is None):
            dimA=np.ndarray(2,dtype=np.int64)
            dimA[0]=dimY[0]
            dimA[1]=dimX[0]
            A=naga_obj_ComplexS2D(self.getContext(),dims=dimA)
        if(device !=Y.getDevice() or
           device !=A.getDevice()):
            raise ValueError(" ger only on the same device")
        cdef cuFloatComplex a
        a.x=alpha.x
        a.y=alpha.y

        A.c_o.ger(a,Y.c_o,1,self.c_o,1,dimY[0])

        return A



    def fft(self, bool inplace=False, bool C2R=False, dest=None, int direction=1):
        """Compute fft, using "cufftExec"

        dest -- naga_obj (default = None)
        dir -- integer (default 1)

        dir: fft's direction
        if dest is None, inplace fft (only available for C2C fft)

        Return dest= fft(self,dir)
        """

        cdef plan

        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd=self.get_Dims()

        cdef int i

        cdef long dims_src[1+1]
        cdef long dims_dest[1+1]
        dims_src[0]=1
        dims_dest[0]=1
        for i in range(1):
            dims_src[i+1]=ds[i]
            dims_dest[i+1]=ds[i]

        if(C2R):
            dims_dest[1]=(dims_dest[1]-1)*2
            dd[1-1]=(dd[1-1]-1)*2

        if(dest is None):
            if(inplace):
                dest=self
            else:
                if(C2R):
                    dest=naga_obj_Float1D(self.getContext(),dims=dd)
                else:
                    dest=naga_obj_ComplexS1D(self.getContext(),dims=dd)

        cdef t_dest=dest.__class__
        try:
            plan=dict_fft_ComplexS1D[t_dest]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return
        output = dest.getData_ptr()

        if(t_dest==self.__class__):
            carma_initfft[cuFloatComplex,cuFloatComplex](dims_src,self.c_o.getPlan(),plan)
            carma_fft[cuFloatComplex,cuFloatComplex](self.c_o.getData(), <cuFloatComplex*><void*>output, direction,handle[0])
        else:
            carma_initfft[cuFloatComplex,float](dims_dest,self.c_o.getPlan(),plan)
            carma_fft[cuFloatComplex,float](self.c_o.getData(), <float*><void*>output, direction,handle[0])

        if(not inplace):
            return dest


##########################################################
#  P-Class `ComplexS2D`
#########################################################


cdef class naga_obj_ComplexS2D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_ComplexS2D obj=None,
                  np.ndarray[ndim=2, dtype=np.complex64_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_ComplexS2D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=2
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[cuFloatComplex](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[cuFloatComplex](ctxt.c,obj.c_o)
            return

        cdef long cdims[2+1]
        cdims[0]=2
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=2]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=2 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 2)")
            for i in range(2):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[cuFloatComplex](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[cuFloatComplex](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(2):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[cuFloatComplex](ctxt.c,cdims,<cuFloatComplex*>data.data)
        else:
            self.c_o=new carma_obj[cuFloatComplex](ctxt.c,cdims,<cuFloatComplex*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((2+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(2+1):
            dims[i]=cdims[i]
        return dims[1:2+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=2, dtype=np.complex64_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=2, dtype=np.complex64_t] data):
        data -- np.complex64: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<cuFloatComplex*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<cuFloatComplex*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.complex64) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex64_t, ndim=2] res=np.zeros(sh,dtype=np.complex64)
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<cuFloatComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.complex64) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex64_t, ndim=2] res=np.zeros(sh, dtype=np.complex64)
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<cuFloatComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_ComplexS2D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_ComplexS2D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_ComplexS2D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_ComplexS2D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_ComplexS2D src):
        """Cublas dot

        src -- naga_obj_ComplexS2D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.complex64_t alpha):
        """Cublas scal

        alpha -- np.complex64: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuFloatComplex a
        a.x=alpha.x
        a.y=alpha.y

        self.c_o.scale(a,1)


    def swap(self, naga_obj_ComplexS2D src):
        """Cublas swap

        src -- naga_obj_ComplexS2D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_ComplexS2D src):
        """Cublas copy

        src -- naga_obj_ComplexS2D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)



    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng(seed)
        self.c_o.prng(b'U')




    def gemv(self, naga_obj_ComplexS1D Vx, np.complex64_t alpha=1,
             naga_obj_ComplexS1D Vy=None, np.complex64_t beta=0 ):
        """Cublas gemv

        Vx -- naga_obj_ComplexS1D
        alpha -- np.complex64 (default = 1)
        Vy -- naga_obj_ComplexS1D (default = None)
        beta -- np.complex64 (default = 0)
        Return  Vy=alpha*self*Vx+beta*Vy
        """
        device=self.getDevice()
        self.context.set_activeDeviceForCpy(device)
        cdef dimM=self.get_Dims()
        if(Vy is None):
            Vy=naga_obj_ComplexS1D(ctxt=self.getContext(),dims=dimM[0:1])

        if( (self.getDevice()!=Vx.getDevice()) or
             self.getDevice()!= Vy.getDevice()):
            raise ValueError("gemv only on the same device")
        cdef cuFloatComplex a
        cdef cuFloatComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y


        #Vy.c_o.gemv(b'N',a, self.c_o,dimM[1],Vx.c_o,1,b,1)
        Vy.c_o.gemv(b'N',a, self.c_o,dimM[0],Vx.c_o,1,b,1)

        return Vy


    def symv(self,naga_obj_ComplexS1D Vx, np.complex64_t alpha=1,
             naga_obj_ComplexS1D Vy=None , np.complex64_t beta=0):
        """Cublas symv

        Vx -- naga_obj_ComplexS1D
        alpha -- np.complex64 (default = 1)
        Vy -- naga_obj_ComplexS1D (default = None)
        beta -- np.complex64 (default = 0)
        Return Vy=alpha*self*Vx+beta*Vy
        """
        cdef device= self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimM=self.get_Dims()
        if(Vy is None):
            #dimY=np.array([1,dimM[1]],dtype=np.int64)
            #Vy=naga_obj_ComplexS2D(ctxt=self.getContext(),dims=dimY)
            Vy=naga_obj_ComplexS2D(ctxt=self.getContext(),dims=dimM)
        if(device !=Vx.getDevice() or
           device !=Vy.getDevice() ):
            raise ValueError("symv only on the same device")

        cdef cuFloatComplex a
        cdef cuFloatComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        Vy.c_o.symv(CUBLAS_FILL_MODE_LOWER, a, self.c_o, dimM[0],
                    Vx.c_o, 1, b, 1)

        return Vy


    def gemm(self,naga_obj_ComplexS2D B, str opA='n', str opB='n',
             np.complex64_t alpha=1, naga_obj_ComplexS2D C=None, np.complex64_t beta=0):
        """Cublas gemm

        B -- naga_obj_ComplexS2D
        opA -- char (default = 'n')
        opB -- char (default = 'n')
        alpha -- np.complex64 (default = 1)
        C -- naga_obj_ComplexS2D (default = None)
        beta  -- np.complex64 (default = 0)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        Return C=alpha opA(self)*opB(B)+beta*C
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            dimC[0]= dimA[0] if opA=='n' else dimA[1]
            dimC[1]= dimB[1] if opB=='n' else dimB[0]
            C=naga_obj_ComplexS2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device!=B.getDevice() or
           device!=C.getDevice() ):
            raise ValueError("gemm only on the same device")

        cdef cuFloatComplex a
        cdef cuFloatComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y


        C.c_o.gemm(ord(opA),ord(opB),a,self.c_o,dimA[0],B.c_o,dimB[0],b,dimC[0])

        return C


    def symm(self, naga_obj_ComplexS2D B, str side='l', np.complex64_t alpha=1,
             naga_obj_ComplexS2D C=None, np.complex64_t beta=0):
        """Cublas symm

        B -- naga_obj_ComplexS2D
        side -- char (default = 'l')
        alpha -- np.complex64 (default =1)
        C -- naga_obj_ComplexS2D (default = None)
        beta -- np.complex64 (default =0)

        return alpha*A*B+beta*C     if side='l'
               alpha*B*A+beta*C     otherwise
        """
        prodSide=CUBLAS_SIDE_LEFT if side=='l' else CUBLAS_SIDE_RIGHT

        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(side=='l'):
                dimC[0]= dimA[0]
                dimC[1]= dimB[1]
            else:
                dimC[0]= dimA[1]
                dimC[1]= dimB[0]

            C=naga_obj_ComplexS2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("symm only on the same device")

        cdef cuFloatComplex a
        cdef cuFloatComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        C.c_o.symm(prodSide,CUBLAS_FILL_MODE_LOWER, a,self.c_o,dimA[0],
                    B.c_o,dimB[0],b,dimC[0])
        return C


    def syrk(self, str opA='n', np.complex64_t alpha=1, naga_obj_ComplexS2D C=None,
             np.complex64_t beta=0):
        """Cublas syrk

        opA -- char (default = 'n')
        alpha -- np.complex64 (default = 1)
        C -- naga_obj_ComplexS2D (default = None)
        beta -- np.complex64 (default = 0)

        opA: transposition on matrix self
        'n': no transposition
        't':transpose matrix
        Return alpha*opA(self)*opA(self)T+beta*C
        """

        cdef device = self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[0]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[1]

            C=naga_obj_ComplexS2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != C.getDevice()):
            raise ValueError("syrk only on the same device")

        cdef cuFloatComplex a
        cdef cuFloatComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        C.c_o.syrk(CUBLAS_FILL_MODE_UPPER, ord(opA),a, self.c_o, dimA[0], b,
                    dimC[0])
        return C


    def syrkx(self, naga_obj_ComplexS2D B, str opA='n',
              np.complex64_t alpha=1, naga_obj_ComplexS2D C=None, np.complex64_t beta=0):
        """Cublas syrkx

        B -- naga_obj_ComplexS2D
        opA -- char (default = 'n')
        apha -- np.complex64 (default = 1)
        C -- naga_obj_ComplexS2D (default = None)
        beta -- np.complex64 (default = 0)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        Return alpha*opA(self)*opB(B)T+beta*C
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimB[0]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimB[1]
            C=naga_obj_ComplexS2D(self.getContext(),dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("syrkx only on the same device")

        cdef cuFloatComplex a
        cdef cuFloatComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        C.c_o.syrkx(CUBLAS_FILL_MODE_UPPER, ord(opA),a,self.c_o,dimA[0],
                    B.c_o, dimB[0],b,dimC[0])

        return C


    def geam(self, naga_obj_ComplexS2D B, np.complex64_t alpha=1, np.complex64_t beta=0,
             str opA='n', str opB='n', naga_obj_ComplexS2D C=None):
        """Cublas geam

        B -- naga_obj_ComplexS2D
        alpha -- np.complex64 (default = 1)
        beta -- np.complex64 (default = 0)
        opA -- char (default = 'n')
        opB -- char (default = 'n')
        C -- naga_obj_ComplexS2D (default = None)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        return C= alpha*opA(self)+beta*opB(B)
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if( C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[1]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[0]
            C=naga_obj_ComplexS2D(self.getContext(),dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("geam only on the same device")


        cdef cuFloatComplex a
        cdef cuFloatComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        C.c_o.geam(ord(opA), ord(opB), a, self.c_o, dimA[0],
                   b, B.c_o, dimB[0], dimC[0])
        return C


    def dgmm(self, naga_obj_ComplexS1D X, str sidec='l',
             naga_obj_ComplexS2D C=None):
        """Cublas dgmm

        X -- naga_obj_ComplexS1D
        side -- char (default = 'l')
        C -- naga_obj_ComplexS2D (default = None)

        Return self*diag(X)     if sidec='l'
               diag(X)*self     otherwise
        """

        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimX=X.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(sidec=='l'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[1]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[0]
            C=naga_obj_ComplexS2D(self.getContext(),dims=dimC)
        else:
            dimC=C.getDims()

        if(device != X.getDevice() or
           device != C.getDevice() ):
            raise ValueError("dgmm only in the same device")

        side=CUBLAS_SIDE_LEFT if sidec=='l' else CUBLAS_SIDE_RIGHT

        C.c_o.dgmm(side,self.c_o,dimA[1], X.c_o, 1, dimC[1])
        return C


    def transpose(self, naga_obj_ComplexS2D dest=None):
        """
        """
        cdef device=self.getDevice()
        self.context.set_activeDeviceForCpy(device)
        cdef dimA=self.getDims()
        if(dest is None):
            dimD=np.ndarray(2,dtype=np.int64)
            dimD[0]=dimA[1]
            dimD[1]=dimA[0]
            dest= naga_obj_ComplexS2D(self.getContext,dims=dimD)
        else:
            dimD=dest.getDims()

        if(device != dest.getDevice()):
            raise ValueError("transpose only on the same device")

        dest.c_o.transpose(self.c_o)
        return dest



    def fft(self, bool inplace=False, bool C2R=False, dest=None, int direction=1):
        """Compute fft, using "cufftExec"

        dest -- naga_obj (default = None)
        dir -- integer (default 1)

        dir: fft's direction
        if dest is None, inplace fft (only available for C2C fft)

        Return dest= fft(self,dir)
        """

        cdef plan

        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd=self.get_Dims()

        cdef int i

        cdef long dims_src[2+1]
        cdef long dims_dest[2+1]
        dims_src[0]=2
        dims_dest[0]=2
        for i in range(2):
            dims_src[i+1]=ds[i]
            dims_dest[i+1]=ds[i]

        if(C2R):
            dims_dest[2]=(dims_dest[2]-1)*2
            dd[2-1]=(dd[2-1]-1)*2

        if(dest is None):
            if(inplace):
                dest=self
            else:
                if(C2R):
                    dest=naga_obj_Float2D(self.getContext(),dims=dd)
                else:
                    dest=naga_obj_ComplexS2D(self.getContext(),dims=dd)

        cdef t_dest=dest.__class__
        try:
            plan=dict_fft_ComplexS2D[t_dest]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return
        output = dest.getData_ptr()

        if(t_dest==self.__class__):
            carma_initfft[cuFloatComplex,cuFloatComplex](dims_src,self.c_o.getPlan(),plan)
            carma_fft[cuFloatComplex,cuFloatComplex](self.c_o.getData(), <cuFloatComplex*><void*>output, direction,handle[0])
        else:
            carma_initfft[cuFloatComplex,float](dims_dest,self.c_o.getPlan(),plan)
            carma_fft[cuFloatComplex,float](self.c_o.getData(), <float*><void*>output, direction,handle[0])

        if(not inplace):
            return dest


##########################################################
#  P-Class `ComplexS3D`
#########################################################


cdef class naga_obj_ComplexS3D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_ComplexS3D obj=None,
                  np.ndarray[ndim=3, dtype=np.complex64_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_ComplexS3D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=3
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[cuFloatComplex](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[cuFloatComplex](ctxt.c,obj.c_o)
            return

        cdef long cdims[3+1]
        cdims[0]=3
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=3]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=3 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 3)")
            for i in range(3):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[cuFloatComplex](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[cuFloatComplex](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(3):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[cuFloatComplex](ctxt.c,cdims,<cuFloatComplex*>data.data)
        else:
            self.c_o=new carma_obj[cuFloatComplex](ctxt.c,cdims,<cuFloatComplex*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((3+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(3+1):
            dims[i]=cdims[i]
        return dims[1:3+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=3, dtype=np.complex64_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=3, dtype=np.complex64_t] data):
        data -- np.complex64: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<cuFloatComplex*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<cuFloatComplex*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.complex64) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex64_t, ndim=3] res=np.zeros(sh,dtype=np.complex64)
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<cuFloatComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.complex64) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex64_t, ndim=3] res=np.zeros(sh, dtype=np.complex64)
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<cuFloatComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_ComplexS3D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_ComplexS3D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_ComplexS3D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_ComplexS3D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_ComplexS3D src):
        """Cublas dot

        src -- naga_obj_ComplexS3D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.complex64_t alpha):
        """Cublas scal

        alpha -- np.complex64: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuFloatComplex a
        a.x=alpha.x
        a.y=alpha.y

        self.c_o.scale(a,1)


    def swap(self, naga_obj_ComplexS3D src):
        """Cublas swap

        src -- naga_obj_ComplexS3D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_ComplexS3D src):
        """Cublas copy

        src -- naga_obj_ComplexS3D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)



    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng(seed)
        self.c_o.prng(b'U')




    def fft(self, bool inplace=False, bool C2R=False, dest=None, int direction=1):
        """Compute fft, using "cufftExec"

        dest -- naga_obj (default = None)
        dir -- integer (default 1)

        dir: fft's direction
        if dest is None, inplace fft (only available for C2C fft)

        Return dest= fft(self,dir)
        """

        cdef plan

        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd=self.get_Dims()

        cdef int i

        cdef long dims_src[3+1]
        cdef long dims_dest[3+1]
        dims_src[0]=3
        dims_dest[0]=3
        for i in range(3):
            dims_src[i+1]=ds[i]
            dims_dest[i+1]=ds[i]

        if(C2R):
            dims_dest[3]=(dims_dest[3]-1)*2
            dd[3-1]=(dd[3-1]-1)*2

        if(dest is None):
            if(inplace):
                dest=self
            else:
                if(C2R):
                    dest=naga_obj_Float3D(self.getContext(),dims=dd)
                else:
                    dest=naga_obj_ComplexS3D(self.getContext(),dims=dd)

        cdef t_dest=dest.__class__
        try:
            plan=dict_fft_ComplexS3D[t_dest]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return
        output = dest.getData_ptr()

        if(t_dest==self.__class__):
            carma_initfft[cuFloatComplex,cuFloatComplex](dims_src,self.c_o.getPlan(),plan)
            carma_fft[cuFloatComplex,cuFloatComplex](self.c_o.getData(), <cuFloatComplex*><void*>output, direction,handle[0])
        else:
            carma_initfft[cuFloatComplex,float](dims_dest,self.c_o.getPlan(),plan)
            carma_fft[cuFloatComplex,float](self.c_o.getData(), <float*><void*>output, direction,handle[0])

        if(not inplace):
            return dest


##########################################################
#  P-Class `ComplexS4D`
#########################################################


cdef class naga_obj_ComplexS4D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_ComplexS4D obj=None,
                  np.ndarray[ndim=4, dtype=np.complex64_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_ComplexS4D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=4
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[cuFloatComplex](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[cuFloatComplex](ctxt.c,obj.c_o)
            return

        cdef long cdims[4+1]
        cdims[0]=4
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=4]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=4 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 4)")
            for i in range(4):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[cuFloatComplex](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[cuFloatComplex](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(4):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[cuFloatComplex](ctxt.c,cdims,<cuFloatComplex*>data.data)
        else:
            self.c_o=new carma_obj[cuFloatComplex](ctxt.c,cdims,<cuFloatComplex*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((4+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(4+1):
            dims[i]=cdims[i]
        return dims[1:4+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=4, dtype=np.complex64_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=4, dtype=np.complex64_t] data):
        data -- np.complex64: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<cuFloatComplex*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<cuFloatComplex*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.complex64) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex64_t, ndim=4] res=np.zeros(sh,dtype=np.complex64)
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<cuFloatComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.complex64) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex64_t, ndim=4] res=np.zeros(sh, dtype=np.complex64)
        cdef np.ndarray[ dtype=np.complex64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<cuFloatComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_ComplexS4D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_ComplexS4D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_ComplexS4D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_ComplexS4D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_ComplexS4D src):
        """Cublas dot

        src -- naga_obj_ComplexS4D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.complex64_t alpha):
        """Cublas scal

        alpha -- np.complex64: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuFloatComplex a
        a.x=alpha.x
        a.y=alpha.y

        self.c_o.scale(a,1)


    def swap(self, naga_obj_ComplexS4D src):
        """Cublas swap

        src -- naga_obj_ComplexS4D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_ComplexS4D src):
        """Cublas copy

        src -- naga_obj_ComplexS4D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)



    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng(seed)
        self.c_o.prng(b'U')




    def fft(self, bool inplace=False, bool C2R=False, dest=None, int direction=1):
        """Compute fft, using "cufftExec"

        dest -- naga_obj (default = None)
        dir -- integer (default 1)

        dir: fft's direction
        if dest is None, inplace fft (only available for C2C fft)

        Return dest= fft(self,dir)
        """

        cdef plan

        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd=self.get_Dims()

        cdef int i

        cdef long dims_src[4+1]
        cdef long dims_dest[4+1]
        dims_src[0]=4
        dims_dest[0]=4
        for i in range(4):
            dims_src[i+1]=ds[i]
            dims_dest[i+1]=ds[i]

        if(C2R):
            dims_dest[4]=(dims_dest[4]-1)*2
            dd[4-1]=(dd[4-1]-1)*2

        if(dest is None):
            if(inplace):
                dest=self
            else:
                if(C2R):
                    dest=naga_obj_Float4D(self.getContext(),dims=dd)
                else:
                    dest=naga_obj_ComplexS4D(self.getContext(),dims=dd)

        cdef t_dest=dest.__class__
        try:
            plan=dict_fft_ComplexS4D[t_dest]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return
        output = dest.getData_ptr()

        if(t_dest==self.__class__):
            carma_initfft[cuFloatComplex,cuFloatComplex](dims_src,self.c_o.getPlan(),plan)
            carma_fft[cuFloatComplex,cuFloatComplex](self.c_o.getData(), <cuFloatComplex*><void*>output, direction,handle[0])
        else:
            carma_initfft[cuFloatComplex,float](dims_dest,self.c_o.getPlan(),plan)
            carma_fft[cuFloatComplex,float](self.c_o.getData(), <float*><void*>output, direction,handle[0])

        if(not inplace):
            return dest


##########################################################
#  P-Class `Double1D`
#########################################################


cdef class naga_obj_Double1D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Double1D obj=None,
                  np.ndarray[ndim=1, dtype=np.float64_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Double1D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=1
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[double](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[double](ctxt.c,obj.c_o)
            return

        cdef long cdims[1+1]
        cdims[0]=1
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=1]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=1 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 1)")
            for i in range(1):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[double](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[double](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(1):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[double](ctxt.c,cdims,<double*>data.data)
        else:
            self.c_o=new carma_obj[double](ctxt.c,cdims,<double*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((1+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(1+1):
            dims[i]=cdims[i]
        return dims[1:1+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=1, dtype=np.float64_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=1, dtype=np.float64_t] data):
        data -- np.float64: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<double*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<double*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.float64) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] res=np.zeros(sh,dtype=np.float64)
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<double*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.float64) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] res=np.zeros(sh, dtype=np.float64)
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<double*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Double1D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Double1D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Double1D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Double1D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_Double1D src):
        """Cublas dot

        src -- naga_obj_Double1D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.float64_t alpha):
        """Cublas scal

        alpha -- np.float64: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef double a
        a=alpha

        self.c_o.scale(a,1)


    def swap(self, naga_obj_Double1D src):
        """Cublas swap

        src -- naga_obj_Double1D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_Double1D src):
        """Cublas copy

        src -- naga_obj_Double1D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)


    def init_prng(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.init_prng(seed)


    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng_host(seed)
        self.c_o.prng_host(b'U')


    def montagn(self, float ta_jolie_maman,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.prng_montagn(ta_jolie_maman)



    def axpy( self, np.float64_t alpha, naga_obj_Double1D dest=None):
        """cublas axpy

        dest -- naga_obj_Double1D
        alpha-- np.float64
        beta -- np.float64
        Return dest=alpha*self +dest
        """
        if(dest==None):
            c=self.getContext()
            dest=naga_obj_Double1D(ctxt=self.getContext(),dims=self.get_Dims())
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("axpy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef double a
        a=alpha

        dest.c_o.axpy(a,self.c_o,1,1)
        return dest


    def ger(self, naga_obj_Double1D Y, np.float64_t alpha=1.0,
            naga_obj_Double2D A=None):
        """Cublas ger

        Y -- naga_obj_Double1D
        alpha -- np.float64 (default = 1)
        A -- naga_obj_Double2D (default = None)

        Return A=alpha*self*t(y)+A
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimX=self.get_Dims()
        cdef dimY=Y.get_Dims()
        if(A is None):
            dimA=np.ndarray(2,dtype=np.int64)
            dimA[0]=dimY[0]
            dimA[1]=dimX[0]
            A=naga_obj_Double2D(self.getContext(),dims=dimA)
        if(device !=Y.getDevice() or
           device !=A.getDevice()):
            raise ValueError(" ger only on the same device")
        cdef double a
        a=alpha

        A.c_o.ger(a,Y.c_o,1,self.c_o,1,dimY[0])

        return A



    def nrm2(self):
        """Cublas nrm2. Return the Euclidean norm """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.nrm2(1)


    def asum(self):
        """Cublas asum. Return the sum of the absolute values of the data's elements """
        self.context.set_activeDeviceForCpy(self.getDevice())
        return self.c_o.asum(1)


    def sum(self):
        """Return the sum of the data's elements """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.sum()



    def fft(self, naga_obj_ComplexD1D dest=None, int direction=1):
        cdef plan
        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd

        cdef long dims_src[1+1]
        dims_src[0]=1
        for i in range(1):
            dims_src[i+1]=ds[i]

        if(dest is None):
            dd=self.get_Dims()
            dd[1-1]=dd[1-1]/2+1
            dest=naga_obj_ComplexD1D(self.getContext(),dims=dd)

        try:
            plan=dict_fft_Double1D[dest.__class__]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return

        output = dest.getData_ptr()
        carma_initfft[double,cuDoubleComplex](dims_src,self.c_o.getPlan(),plan)
        carma_fft[double,cuDoubleComplex](self.c_o.getData(), <cuDoubleComplex*><void*>output, direction,handle[0])

        return dest


##########################################################
#  P-Class `Double2D`
#########################################################


cdef class naga_obj_Double2D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Double2D obj=None,
                  np.ndarray[ndim=2, dtype=np.float64_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Double2D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=2
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[double](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[double](ctxt.c,obj.c_o)
            return

        cdef long cdims[2+1]
        cdims[0]=2
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=2]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=2 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 2)")
            for i in range(2):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[double](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[double](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(2):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[double](ctxt.c,cdims,<double*>data.data)
        else:
            self.c_o=new carma_obj[double](ctxt.c,cdims,<double*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((2+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(2+1):
            dims[i]=cdims[i]
        return dims[1:2+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=2, dtype=np.float64_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=2, dtype=np.float64_t] data):
        data -- np.float64: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<double*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<double*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.float64) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float64_t, ndim=2] res=np.zeros(sh,dtype=np.float64)
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<double*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.float64) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float64_t, ndim=2] res=np.zeros(sh, dtype=np.float64)
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<double*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Double2D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Double2D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Double2D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Double2D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_Double2D src):
        """Cublas dot

        src -- naga_obj_Double2D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.float64_t alpha):
        """Cublas scal

        alpha -- np.float64: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef double a
        a=alpha

        self.c_o.scale(a,1)


    def swap(self, naga_obj_Double2D src):
        """Cublas swap

        src -- naga_obj_Double2D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_Double2D src):
        """Cublas copy

        src -- naga_obj_Double2D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)


    def init_prng(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.init_prng(seed)


    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng_host(seed)
        self.c_o.prng_host(b'U')


    def montagn(self, float ta_jolie_maman,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.prng_montagn(ta_jolie_maman)



    def gemv(self, naga_obj_Double1D Vx, np.float64_t alpha=1,
             naga_obj_Double1D Vy=None, np.float64_t beta=0 ):
        """Cublas gemv

        Vx -- naga_obj_Double1D
        alpha -- np.float64 (default = 1)
        Vy -- naga_obj_Double1D (default = None)
        beta -- np.float64 (default = 0)
        Return  Vy=alpha*self*Vx+beta*Vy
        """
        device=self.getDevice()
        self.context.set_activeDeviceForCpy(device)
        cdef dimM=self.get_Dims()
        if(Vy is None):
            Vy=naga_obj_Double1D(ctxt=self.getContext(),dims=dimM[0:1])

        if( (self.getDevice()!=Vx.getDevice()) or
             self.getDevice()!= Vy.getDevice()):
            raise ValueError("gemv only on the same device")
        cdef double a
        cdef double b
        a=alpha
        b=beta


        #Vy.c_o.gemv(b'N',a, self.c_o,dimM[1],Vx.c_o,1,b,1)
        Vy.c_o.gemv(b'N',a, self.c_o,dimM[0],Vx.c_o,1,b,1)

        return Vy


    def symv(self,naga_obj_Double1D Vx, np.float64_t alpha=1,
             naga_obj_Double1D Vy=None , np.float64_t beta=0):
        """Cublas symv

        Vx -- naga_obj_Double1D
        alpha -- np.float64 (default = 1)
        Vy -- naga_obj_Double1D (default = None)
        beta -- np.float64 (default = 0)
        Return Vy=alpha*self*Vx+beta*Vy
        """
        cdef device= self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimM=self.get_Dims()
        if(Vy is None):
            #dimY=np.array([1,dimM[1]],dtype=np.int64)
            #Vy=naga_obj_Double2D(ctxt=self.getContext(),dims=dimY)
            Vy=naga_obj_Double2D(ctxt=self.getContext(),dims=dimM)
        if(device !=Vx.getDevice() or
           device !=Vy.getDevice() ):
            raise ValueError("symv only on the same device")

        cdef double a
        cdef double b
        a=alpha
        b=beta

        Vy.c_o.symv(CUBLAS_FILL_MODE_LOWER, a, self.c_o, dimM[0],
                    Vx.c_o, 1, b, 1)

        return Vy


    def gemm(self,naga_obj_Double2D B, str opA='n', str opB='n',
             np.float64_t alpha=1, naga_obj_Double2D C=None, np.float64_t beta=0):
        """Cublas gemm

        B -- naga_obj_Double2D
        opA -- char (default = 'n')
        opB -- char (default = 'n')
        alpha -- np.float64 (default = 1)
        C -- naga_obj_Double2D (default = None)
        beta  -- np.float64 (default = 0)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        Return C=alpha opA(self)*opB(B)+beta*C
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            dimC[0]= dimA[0] if opA=='n' else dimA[1]
            dimC[1]= dimB[1] if opB=='n' else dimB[0]
            C=naga_obj_Double2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device!=B.getDevice() or
           device!=C.getDevice() ):
            raise ValueError("gemm only on the same device")

        cdef double a
        cdef double b
        a=alpha
        b=beta


        C.c_o.gemm(ord(opA),ord(opB),a,self.c_o,dimA[0],B.c_o,dimB[0],b,dimC[0])

        return C


    def symm(self, naga_obj_Double2D B, str side='l', np.float64_t alpha=1,
             naga_obj_Double2D C=None, np.float64_t beta=0):
        """Cublas symm

        B -- naga_obj_Double2D
        side -- char (default = 'l')
        alpha -- np.float64 (default =1)
        C -- naga_obj_Double2D (default = None)
        beta -- np.float64 (default =0)

        return alpha*A*B+beta*C     if side='l'
               alpha*B*A+beta*C     otherwise
        """
        prodSide=CUBLAS_SIDE_LEFT if side=='l' else CUBLAS_SIDE_RIGHT

        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(side=='l'):
                dimC[0]= dimA[0]
                dimC[1]= dimB[1]
            else:
                dimC[0]= dimA[1]
                dimC[1]= dimB[0]

            C=naga_obj_Double2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("symm only on the same device")

        cdef double a
        cdef double b
        a=alpha
        b=beta

        C.c_o.symm(prodSide,CUBLAS_FILL_MODE_LOWER, a,self.c_o,dimA[0],
                    B.c_o,dimB[0],b,dimC[0])
        return C


    def syrk(self, str opA='n', np.float64_t alpha=1, naga_obj_Double2D C=None,
             np.float64_t beta=0):
        """Cublas syrk

        opA -- char (default = 'n')
        alpha -- np.float64 (default = 1)
        C -- naga_obj_Double2D (default = None)
        beta -- np.float64 (default = 0)

        opA: transposition on matrix self
        'n': no transposition
        't':transpose matrix
        Return alpha*opA(self)*opA(self)T+beta*C
        """

        cdef device = self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[0]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[1]

            C=naga_obj_Double2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != C.getDevice()):
            raise ValueError("syrk only on the same device")

        cdef double a
        cdef double b
        a=alpha
        b=beta

        C.c_o.syrk(CUBLAS_FILL_MODE_UPPER, ord(opA),a, self.c_o, dimA[0], b,
                    dimC[0])
        return C


    def syrkx(self, naga_obj_Double2D B, str opA='n',
              np.float64_t alpha=1, naga_obj_Double2D C=None, np.float64_t beta=0):
        """Cublas syrkx

        B -- naga_obj_Double2D
        opA -- char (default = 'n')
        apha -- np.float64 (default = 1)
        C -- naga_obj_Double2D (default = None)
        beta -- np.float64 (default = 0)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        Return alpha*opA(self)*opB(B)T+beta*C
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimB[0]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimB[1]
            C=naga_obj_Double2D(self.getContext(),dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("syrkx only on the same device")

        cdef double a
        cdef double b
        a=alpha
        b=beta

        C.c_o.syrkx(CUBLAS_FILL_MODE_UPPER, ord(opA),a,self.c_o,dimA[0],
                    B.c_o, dimB[0],b,dimC[0])

        return C


    def geam(self, naga_obj_Double2D B, np.float64_t alpha=1, np.float64_t beta=0,
             str opA='n', str opB='n', naga_obj_Double2D C=None):
        """Cublas geam

        B -- naga_obj_Double2D
        alpha -- np.float64 (default = 1)
        beta -- np.float64 (default = 0)
        opA -- char (default = 'n')
        opB -- char (default = 'n')
        C -- naga_obj_Double2D (default = None)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        return C= alpha*opA(self)+beta*opB(B)
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if( C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[1]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[0]
            C=naga_obj_Double2D(self.getContext(),dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("geam only on the same device")


        cdef double a
        cdef double b
        a=alpha
        b=beta

        C.c_o.geam(ord(opA), ord(opB), a, self.c_o, dimA[0],
                   b, B.c_o, dimB[0], dimC[0])
        return C


    def dgmm(self, naga_obj_Double1D X, str sidec='l',
             naga_obj_Double2D C=None):
        """Cublas dgmm

        X -- naga_obj_Double1D
        side -- char (default = 'l')
        C -- naga_obj_Double2D (default = None)

        Return self*diag(X)     if sidec='l'
               diag(X)*self     otherwise
        """

        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimX=X.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(sidec=='l'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[1]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[0]
            C=naga_obj_Double2D(self.getContext(),dims=dimC)
        else:
            dimC=C.getDims()

        if(device != X.getDevice() or
           device != C.getDevice() ):
            raise ValueError("dgmm only in the same device")

        side=CUBLAS_SIDE_LEFT if sidec=='l' else CUBLAS_SIDE_RIGHT

        C.c_o.dgmm(side,self.c_o,dimA[1], X.c_o, 1, dimC[1])
        return C


    def transpose(self, naga_obj_Double2D dest=None):
        """
        """
        cdef device=self.getDevice()
        self.context.set_activeDeviceForCpy(device)
        cdef dimA=self.getDims()
        if(dest is None):
            dimD=np.ndarray(2,dtype=np.int64)
            dimD[0]=dimA[1]
            dimD[1]=dimA[0]
            dest= naga_obj_Double2D(self.getContext,dims=dimD)
        else:
            dimD=dest.getDims()

        if(device != dest.getDevice()):
            raise ValueError("transpose only on the same device")

        dest.c_o.transpose(self.c_o)
        return dest



    def nrm2(self):
        """Cublas nrm2. Return the Euclidean norm """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.nrm2(1)


    def asum(self):
        """Cublas asum. Return the sum of the absolute values of the data's elements """
        self.context.set_activeDeviceForCpy(self.getDevice())
        return self.c_o.asum(1)


    def sum(self):
        """Return the sum of the data's elements """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.sum()



    def fft(self, naga_obj_ComplexD2D dest=None, int direction=1):
        cdef plan
        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd

        cdef long dims_src[2+1]
        dims_src[0]=2
        for i in range(2):
            dims_src[i+1]=ds[i]

        if(dest is None):
            dd=self.get_Dims()
            dd[2-1]=dd[2-1]/2+1
            dest=naga_obj_ComplexD2D(self.getContext(),dims=dd)

        try:
            plan=dict_fft_Double2D[dest.__class__]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return

        output = dest.getData_ptr()
        carma_initfft[double,cuDoubleComplex](dims_src,self.c_o.getPlan(),plan)
        carma_fft[double,cuDoubleComplex](self.c_o.getData(), <cuDoubleComplex*><void*>output, direction,handle[0])

        return dest


##########################################################
#  P-Class `Double3D`
#########################################################


cdef class naga_obj_Double3D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Double3D obj=None,
                  np.ndarray[ndim=3, dtype=np.float64_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Double3D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=3
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[double](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[double](ctxt.c,obj.c_o)
            return

        cdef long cdims[3+1]
        cdims[0]=3
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=3]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=3 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 3)")
            for i in range(3):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[double](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[double](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(3):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[double](ctxt.c,cdims,<double*>data.data)
        else:
            self.c_o=new carma_obj[double](ctxt.c,cdims,<double*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((3+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(3+1):
            dims[i]=cdims[i]
        return dims[1:3+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=3, dtype=np.float64_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=3, dtype=np.float64_t] data):
        data -- np.float64: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<double*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<double*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.float64) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float64_t, ndim=3] res=np.zeros(sh,dtype=np.float64)
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<double*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.float64) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float64_t, ndim=3] res=np.zeros(sh, dtype=np.float64)
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<double*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Double3D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Double3D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Double3D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Double3D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_Double3D src):
        """Cublas dot

        src -- naga_obj_Double3D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.float64_t alpha):
        """Cublas scal

        alpha -- np.float64: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef double a
        a=alpha

        self.c_o.scale(a,1)


    def swap(self, naga_obj_Double3D src):
        """Cublas swap

        src -- naga_obj_Double3D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_Double3D src):
        """Cublas copy

        src -- naga_obj_Double3D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)


    def init_prng(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.init_prng(seed)


    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng_host(seed)
        self.c_o.prng_host(b'U')


    def montagn(self, float ta_jolie_maman,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.prng_montagn(ta_jolie_maman)



    def nrm2(self):
        """Cublas nrm2. Return the Euclidean norm """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.nrm2(1)


    def asum(self):
        """Cublas asum. Return the sum of the absolute values of the data's elements """
        self.context.set_activeDeviceForCpy(self.getDevice())
        return self.c_o.asum(1)


    def sum(self):
        """Return the sum of the data's elements """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.sum()



    def fft(self, naga_obj_ComplexD3D dest=None, int direction=1):
        cdef plan
        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd

        cdef long dims_src[3+1]
        dims_src[0]=3
        for i in range(3):
            dims_src[i+1]=ds[i]

        if(dest is None):
            dd=self.get_Dims()
            dd[3-1]=dd[3-1]/2+1
            dest=naga_obj_ComplexD3D(self.getContext(),dims=dd)

        try:
            plan=dict_fft_Double3D[dest.__class__]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return

        output = dest.getData_ptr()
        carma_initfft[double,cuDoubleComplex](dims_src,self.c_o.getPlan(),plan)
        carma_fft[double,cuDoubleComplex](self.c_o.getData(), <cuDoubleComplex*><void*>output, direction,handle[0])

        return dest


##########################################################
#  P-Class `Double4D`
#########################################################


cdef class naga_obj_Double4D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_Double4D obj=None,
                  np.ndarray[ndim=4, dtype=np.float64_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_Double4D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=4
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[double](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[double](ctxt.c,obj.c_o)
            return

        cdef long cdims[4+1]
        cdims[0]=4
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=4]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=4 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 4)")
            for i in range(4):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[double](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[double](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(4):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[double](ctxt.c,cdims,<double*>data.data)
        else:
            self.c_o=new carma_obj[double](ctxt.c,cdims,<double*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((4+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(4+1):
            dims[i]=cdims[i]
        return dims[1:4+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=4, dtype=np.float64_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=4, dtype=np.float64_t] data):
        data -- np.float64: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<double*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<double*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.float64) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float64_t, ndim=4] res=np.zeros(sh,dtype=np.float64)
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<double*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.float64) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.float64_t, ndim=4] res=np.zeros(sh, dtype=np.float64)
        cdef np.ndarray[ dtype=np.float64_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.float64)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<double*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_Double4D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_Double4D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_Double4D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_Double4D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_Double4D src):
        """Cublas dot

        src -- naga_obj_Double4D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.float64_t alpha):
        """Cublas scal

        alpha -- np.float64: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef double a
        a=alpha

        self.c_o.scale(a,1)


    def swap(self, naga_obj_Double4D src):
        """Cublas swap

        src -- naga_obj_Double4D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_Double4D src):
        """Cublas copy

        src -- naga_obj_Double4D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)


    def init_prng(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.init_prng(seed)


    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng_host(seed)
        self.c_o.prng_host(b'U')


    def montagn(self, float ta_jolie_maman,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        self.c_o.prng_montagn(ta_jolie_maman)



    def nrm2(self):
        """Cublas nrm2. Return the Euclidean norm """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.nrm2(1)


    def asum(self):
        """Cublas asum. Return the sum of the absolute values of the data's elements """
        self.context.set_activeDeviceForCpy(self.getDevice())
        return self.c_o.asum(1)


    def sum(self):
        """Return the sum of the data's elements """
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.sum()



    def fft(self, naga_obj_ComplexD4D dest=None, int direction=1):
        cdef plan
        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd

        cdef long dims_src[4+1]
        dims_src[0]=4
        for i in range(4):
            dims_src[i+1]=ds[i]

        if(dest is None):
            dd=self.get_Dims()
            dd[4-1]=dd[4-1]/2+1
            dest=naga_obj_ComplexD4D(self.getContext(),dims=dd)

        try:
            plan=dict_fft_Double4D[dest.__class__]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return

        output = dest.getData_ptr()
        carma_initfft[double,cuDoubleComplex](dims_src,self.c_o.getPlan(),plan)
        carma_fft[double,cuDoubleComplex](self.c_o.getData(), <cuDoubleComplex*><void*>output, direction,handle[0])

        return dest


##########################################################
#  P-Class `ComplexD1D`
#########################################################


cdef class naga_obj_ComplexD1D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_ComplexD1D obj=None,
                  np.ndarray[ndim=1, dtype=np.complex128_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_ComplexD1D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=1
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[cuDoubleComplex](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[cuDoubleComplex](ctxt.c,obj.c_o)
            return

        cdef long cdims[1+1]
        cdims[0]=1
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=1]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=1 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 1)")
            for i in range(1):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[cuDoubleComplex](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[cuDoubleComplex](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(1):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[cuDoubleComplex](ctxt.c,cdims,<cuDoubleComplex*>data.data)
        else:
            self.c_o=new carma_obj[cuDoubleComplex](ctxt.c,cdims,<cuDoubleComplex*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((1+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(1+1):
            dims[i]=cdims[i]
        return dims[1:1+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=1, dtype=np.complex128_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=1, dtype=np.complex128_t] data):
        data -- np.complex128: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<cuDoubleComplex*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<cuDoubleComplex*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.complex128) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] res=np.zeros(sh,dtype=np.complex128)
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex128)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<cuDoubleComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.complex128) of 1 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] res=np.zeros(sh, dtype=np.complex128)
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex128)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<cuDoubleComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_ComplexD1D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_ComplexD1D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_ComplexD1D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_ComplexD1D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_ComplexD1D src):
        """Cublas dot

        src -- naga_obj_ComplexD1D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.complex128_t alpha):
        """Cublas scal

        alpha -- np.complex128: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuDoubleComplex a
        a.x=alpha.x
        a.y=alpha.y

        self.c_o.scale(a,1)


    def swap(self, naga_obj_ComplexD1D src):
        """Cublas swap

        src -- naga_obj_ComplexD1D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_ComplexD1D src):
        """Cublas copy

        src -- naga_obj_ComplexD1D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)



    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng(seed)
        self.c_o.prng(b'U')




    def axpy( self, np.complex128_t alpha, naga_obj_ComplexD1D dest=None):
        """cublas axpy

        dest -- naga_obj_ComplexD1D
        alpha-- np.complex128
        beta -- np.complex128
        Return dest=alpha*self +dest
        """
        if(dest==None):
            c=self.getContext()
            dest=naga_obj_ComplexD1D(ctxt=self.getContext(),dims=self.get_Dims())
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("axpy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuDoubleComplex a
        a.x=alpha.x
        a.y=alpha.y

        dest.c_o.axpy(a,self.c_o,1,1)
        return dest


    def ger(self, naga_obj_ComplexD1D Y, np.complex128_t alpha=1.0,
            naga_obj_ComplexD2D A=None):
        """Cublas ger

        Y -- naga_obj_ComplexD1D
        alpha -- np.complex128 (default = 1)
        A -- naga_obj_ComplexD2D (default = None)

        Return A=alpha*self*t(y)+A
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimX=self.get_Dims()
        cdef dimY=Y.get_Dims()
        if(A is None):
            dimA=np.ndarray(2,dtype=np.int64)
            dimA[0]=dimY[0]
            dimA[1]=dimX[0]
            A=naga_obj_ComplexD2D(self.getContext(),dims=dimA)
        if(device !=Y.getDevice() or
           device !=A.getDevice()):
            raise ValueError(" ger only on the same device")
        cdef cuDoubleComplex a
        a.x=alpha.x
        a.y=alpha.y

        A.c_o.ger(a,Y.c_o,1,self.c_o,1,dimY[0])

        return A



    def fft(self, bool inplace=False, bool C2R=False, dest=None, int direction=1):
        """Compute fft, using "cufftExec"

        dest -- naga_obj (default = None)
        dir -- integer (default 1)

        dir: fft's direction
        if dest is None, inplace fft (only available for C2C fft)

        Return dest= fft(self,dir)
        """

        cdef plan

        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd=self.get_Dims()

        cdef int i

        cdef long dims_src[1+1]
        cdef long dims_dest[1+1]
        dims_src[0]=1
        dims_dest[0]=1
        for i in range(1):
            dims_src[i+1]=ds[i]
            dims_dest[i+1]=ds[i]

        if(C2R):
            dims_dest[1]=(dims_dest[1]-1)*2
            dd[1-1]=(dd[1-1]-1)*2

        if(dest is None):
            if(inplace):
                dest=self
            else:
                if(C2R):
                    dest=naga_obj_Double1D(self.getContext(),dims=dd)
                else:
                    dest=naga_obj_ComplexD1D(self.getContext(),dims=dd)

        cdef t_dest=dest.__class__
        try:
            plan=dict_fft_ComplexD1D[t_dest]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return
        output = dest.getData_ptr()

        if(t_dest==self.__class__):
            carma_initfft[cuDoubleComplex,cuDoubleComplex](dims_src,self.c_o.getPlan(),plan)
            carma_fft[cuDoubleComplex,cuDoubleComplex](self.c_o.getData(), <cuDoubleComplex*><void*>output, direction,handle[0])
        else:
            carma_initfft[cuDoubleComplex,double](dims_dest,self.c_o.getPlan(),plan)
            carma_fft[cuDoubleComplex,double](self.c_o.getData(), <double*><void*>output, direction,handle[0])

        if(not inplace):
            return dest


##########################################################
#  P-Class `ComplexD2D`
#########################################################


cdef class naga_obj_ComplexD2D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_ComplexD2D obj=None,
                  np.ndarray[ndim=2, dtype=np.complex128_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_ComplexD2D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=2
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[cuDoubleComplex](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[cuDoubleComplex](ctxt.c,obj.c_o)
            return

        cdef long cdims[2+1]
        cdims[0]=2
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=2]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=2 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 2)")
            for i in range(2):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[cuDoubleComplex](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[cuDoubleComplex](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(2):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[cuDoubleComplex](ctxt.c,cdims,<cuDoubleComplex*>data.data)
        else:
            self.c_o=new carma_obj[cuDoubleComplex](ctxt.c,cdims,<cuDoubleComplex*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((2+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(2+1):
            dims[i]=cdims[i]
        return dims[1:2+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=2, dtype=np.complex128_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=2, dtype=np.complex128_t] data):
        data -- np.complex128: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<cuDoubleComplex*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<cuDoubleComplex*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.complex128) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex128_t, ndim=2] res=np.zeros(sh,dtype=np.complex128)
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex128)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<cuDoubleComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.complex128) of 2 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex128_t, ndim=2] res=np.zeros(sh, dtype=np.complex128)
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex128)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<cuDoubleComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_ComplexD2D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_ComplexD2D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_ComplexD2D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_ComplexD2D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_ComplexD2D src):
        """Cublas dot

        src -- naga_obj_ComplexD2D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.complex128_t alpha):
        """Cublas scal

        alpha -- np.complex128: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuDoubleComplex a
        a.x=alpha.x
        a.y=alpha.y

        self.c_o.scale(a,1)


    def swap(self, naga_obj_ComplexD2D src):
        """Cublas swap

        src -- naga_obj_ComplexD2D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_ComplexD2D src):
        """Cublas copy

        src -- naga_obj_ComplexD2D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)



    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng(seed)
        self.c_o.prng(b'U')




    def gemv(self, naga_obj_ComplexD1D Vx, np.complex128_t alpha=1,
             naga_obj_ComplexD1D Vy=None, np.complex128_t beta=0 ):
        """Cublas gemv

        Vx -- naga_obj_ComplexD1D
        alpha -- np.complex128 (default = 1)
        Vy -- naga_obj_ComplexD1D (default = None)
        beta -- np.complex128 (default = 0)
        Return  Vy=alpha*self*Vx+beta*Vy
        """
        device=self.getDevice()
        self.context.set_activeDeviceForCpy(device)
        cdef dimM=self.get_Dims()
        if(Vy is None):
            Vy=naga_obj_ComplexD1D(ctxt=self.getContext(),dims=dimM[0:1])

        if( (self.getDevice()!=Vx.getDevice()) or
             self.getDevice()!= Vy.getDevice()):
            raise ValueError("gemv only on the same device")
        cdef cuDoubleComplex a
        cdef cuDoubleComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y


        #Vy.c_o.gemv(b'N',a, self.c_o,dimM[1],Vx.c_o,1,b,1)
        Vy.c_o.gemv(b'N',a, self.c_o,dimM[0],Vx.c_o,1,b,1)

        return Vy


    def symv(self,naga_obj_ComplexD1D Vx, np.complex128_t alpha=1,
             naga_obj_ComplexD1D Vy=None , np.complex128_t beta=0):
        """Cublas symv

        Vx -- naga_obj_ComplexD1D
        alpha -- np.complex128 (default = 1)
        Vy -- naga_obj_ComplexD1D (default = None)
        beta -- np.complex128 (default = 0)
        Return Vy=alpha*self*Vx+beta*Vy
        """
        cdef device= self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimM=self.get_Dims()
        if(Vy is None):
            #dimY=np.array([1,dimM[1]],dtype=np.int64)
            #Vy=naga_obj_ComplexD2D(ctxt=self.getContext(),dims=dimY)
            Vy=naga_obj_ComplexD2D(ctxt=self.getContext(),dims=dimM)
        if(device !=Vx.getDevice() or
           device !=Vy.getDevice() ):
            raise ValueError("symv only on the same device")

        cdef cuDoubleComplex a
        cdef cuDoubleComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        Vy.c_o.symv(CUBLAS_FILL_MODE_LOWER, a, self.c_o, dimM[0],
                    Vx.c_o, 1, b, 1)

        return Vy


    def gemm(self,naga_obj_ComplexD2D B, str opA='n', str opB='n',
             np.complex128_t alpha=1, naga_obj_ComplexD2D C=None, np.complex128_t beta=0):
        """Cublas gemm

        B -- naga_obj_ComplexD2D
        opA -- char (default = 'n')
        opB -- char (default = 'n')
        alpha -- np.complex128 (default = 1)
        C -- naga_obj_ComplexD2D (default = None)
        beta  -- np.complex128 (default = 0)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        Return C=alpha opA(self)*opB(B)+beta*C
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            dimC[0]= dimA[0] if opA=='n' else dimA[1]
            dimC[1]= dimB[1] if opB=='n' else dimB[0]
            C=naga_obj_ComplexD2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device!=B.getDevice() or
           device!=C.getDevice() ):
            raise ValueError("gemm only on the same device")

        cdef cuDoubleComplex a
        cdef cuDoubleComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y


        C.c_o.gemm(ord(opA),ord(opB),a,self.c_o,dimA[0],B.c_o,dimB[0],b,dimC[0])

        return C


    def symm(self, naga_obj_ComplexD2D B, str side='l', np.complex128_t alpha=1,
             naga_obj_ComplexD2D C=None, np.complex128_t beta=0):
        """Cublas symm

        B -- naga_obj_ComplexD2D
        side -- char (default = 'l')
        alpha -- np.complex128 (default =1)
        C -- naga_obj_ComplexD2D (default = None)
        beta -- np.complex128 (default =0)

        return alpha*A*B+beta*C     if side='l'
               alpha*B*A+beta*C     otherwise
        """
        prodSide=CUBLAS_SIDE_LEFT if side=='l' else CUBLAS_SIDE_RIGHT

        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(side=='l'):
                dimC[0]= dimA[0]
                dimC[1]= dimB[1]
            else:
                dimC[0]= dimA[1]
                dimC[1]= dimB[0]

            C=naga_obj_ComplexD2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("symm only on the same device")

        cdef cuDoubleComplex a
        cdef cuDoubleComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        C.c_o.symm(prodSide,CUBLAS_FILL_MODE_LOWER, a,self.c_o,dimA[0],
                    B.c_o,dimB[0],b,dimC[0])
        return C


    def syrk(self, str opA='n', np.complex128_t alpha=1, naga_obj_ComplexD2D C=None,
             np.complex128_t beta=0):
        """Cublas syrk

        opA -- char (default = 'n')
        alpha -- np.complex128 (default = 1)
        C -- naga_obj_ComplexD2D (default = None)
        beta -- np.complex128 (default = 0)

        opA: transposition on matrix self
        'n': no transposition
        't':transpose matrix
        Return alpha*opA(self)*opA(self)T+beta*C
        """

        cdef device = self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[0]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[1]

            C=naga_obj_ComplexD2D(ctxt=self.getContext(), dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != C.getDevice()):
            raise ValueError("syrk only on the same device")

        cdef cuDoubleComplex a
        cdef cuDoubleComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        C.c_o.syrk(CUBLAS_FILL_MODE_UPPER, ord(opA),a, self.c_o, dimA[0], b,
                    dimC[0])
        return C


    def syrkx(self, naga_obj_ComplexD2D B, str opA='n',
              np.complex128_t alpha=1, naga_obj_ComplexD2D C=None, np.complex128_t beta=0):
        """Cublas syrkx

        B -- naga_obj_ComplexD2D
        opA -- char (default = 'n')
        apha -- np.complex128 (default = 1)
        C -- naga_obj_ComplexD2D (default = None)
        beta -- np.complex128 (default = 0)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        Return alpha*opA(self)*opB(B)T+beta*C
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimB[0]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimB[1]
            C=naga_obj_ComplexD2D(self.getContext(),dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("syrkx only on the same device")

        cdef cuDoubleComplex a
        cdef cuDoubleComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        C.c_o.syrkx(CUBLAS_FILL_MODE_UPPER, ord(opA),a,self.c_o,dimA[0],
                    B.c_o, dimB[0],b,dimC[0])

        return C


    def geam(self, naga_obj_ComplexD2D B, np.complex128_t alpha=1, np.complex128_t beta=0,
             str opA='n', str opB='n', naga_obj_ComplexD2D C=None):
        """Cublas geam

        B -- naga_obj_ComplexD2D
        alpha -- np.complex128 (default = 1)
        beta -- np.complex128 (default = 0)
        opA -- char (default = 'n')
        opB -- char (default = 'n')
        C -- naga_obj_ComplexD2D (default = None)

        opA (opB): transposition on matrix self (B),
        'n': no transposition
        't':transpose matrix
        return C= alpha*opA(self)+beta*opB(B)
        """
        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimB=B.get_Dims()
        if( C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(opA=='n'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[1]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[0]
            C=naga_obj_ComplexD2D(self.getContext(),dims=dimC)
        else:
            dimC=C.get_Dims()

        if(device != B.getDevice() or
           device != C.getDevice() ):
            raise ValueError("geam only on the same device")


        cdef cuDoubleComplex a
        cdef cuDoubleComplex b
        a.x=alpha.x
        a.y=alpha.y
        b.x=beta.x
        b.y=beta.y

        C.c_o.geam(ord(opA), ord(opB), a, self.c_o, dimA[0],
                   b, B.c_o, dimB[0], dimC[0])
        return C


    def dgmm(self, naga_obj_ComplexD1D X, str sidec='l',
             naga_obj_ComplexD2D C=None):
        """Cublas dgmm

        X -- naga_obj_ComplexD1D
        side -- char (default = 'l')
        C -- naga_obj_ComplexD2D (default = None)

        Return self*diag(X)     if sidec='l'
               diag(X)*self     otherwise
        """

        cdef device=self.getDevice()
        self.context.set_activeDevice(device)
        cdef dimA=self.get_Dims()
        cdef dimX=X.get_Dims()
        if(C is None):
            dimC=np.ndarray(2,dtype=np.int64)
            if(sidec=='l'):
                dimC[0]=dimA[0]
                dimC[1]=dimA[1]
            else:
                dimC[0]=dimA[1]
                dimC[1]=dimA[0]
            C=naga_obj_ComplexD2D(self.getContext(),dims=dimC)
        else:
            dimC=C.getDims()

        if(device != X.getDevice() or
           device != C.getDevice() ):
            raise ValueError("dgmm only in the same device")

        side=CUBLAS_SIDE_LEFT if sidec=='l' else CUBLAS_SIDE_RIGHT

        C.c_o.dgmm(side,self.c_o,dimA[1], X.c_o, 1, dimC[1])
        return C


    def transpose(self, naga_obj_ComplexD2D dest=None):
        """
        """
        cdef device=self.getDevice()
        self.context.set_activeDeviceForCpy(device)
        cdef dimA=self.getDims()
        if(dest is None):
            dimD=np.ndarray(2,dtype=np.int64)
            dimD[0]=dimA[1]
            dimD[1]=dimA[0]
            dest= naga_obj_ComplexD2D(self.getContext,dims=dimD)
        else:
            dimD=dest.getDims()

        if(device != dest.getDevice()):
            raise ValueError("transpose only on the same device")

        dest.c_o.transpose(self.c_o)
        return dest



    def fft(self, bool inplace=False, bool C2R=False, dest=None, int direction=1):
        """Compute fft, using "cufftExec"

        dest -- naga_obj (default = None)
        dir -- integer (default 1)

        dir: fft's direction
        if dest is None, inplace fft (only available for C2C fft)

        Return dest= fft(self,dir)
        """

        cdef plan

        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd=self.get_Dims()

        cdef int i

        cdef long dims_src[2+1]
        cdef long dims_dest[2+1]
        dims_src[0]=2
        dims_dest[0]=2
        for i in range(2):
            dims_src[i+1]=ds[i]
            dims_dest[i+1]=ds[i]

        if(C2R):
            dims_dest[2]=(dims_dest[2]-1)*2
            dd[2-1]=(dd[2-1]-1)*2

        if(dest is None):
            if(inplace):
                dest=self
            else:
                if(C2R):
                    dest=naga_obj_Double2D(self.getContext(),dims=dd)
                else:
                    dest=naga_obj_ComplexD2D(self.getContext(),dims=dd)

        cdef t_dest=dest.__class__
        try:
            plan=dict_fft_ComplexD2D[t_dest]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return
        output = dest.getData_ptr()

        if(t_dest==self.__class__):
            carma_initfft[cuDoubleComplex,cuDoubleComplex](dims_src,self.c_o.getPlan(),plan)
            carma_fft[cuDoubleComplex,cuDoubleComplex](self.c_o.getData(), <cuDoubleComplex*><void*>output, direction,handle[0])
        else:
            carma_initfft[cuDoubleComplex,double](dims_dest,self.c_o.getPlan(),plan)
            carma_fft[cuDoubleComplex,double](self.c_o.getData(), <double*><void*>output, direction,handle[0])

        if(not inplace):
            return dest


##########################################################
#  P-Class `ComplexD3D`
#########################################################


cdef class naga_obj_ComplexD3D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_ComplexD3D obj=None,
                  np.ndarray[ndim=3, dtype=np.complex128_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_ComplexD3D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=3
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[cuDoubleComplex](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[cuDoubleComplex](ctxt.c,obj.c_o)
            return

        cdef long cdims[3+1]
        cdims[0]=3
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=3]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=3 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 3)")
            for i in range(3):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[cuDoubleComplex](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[cuDoubleComplex](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(3):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[cuDoubleComplex](ctxt.c,cdims,<cuDoubleComplex*>data.data)
        else:
            self.c_o=new carma_obj[cuDoubleComplex](ctxt.c,cdims,<cuDoubleComplex*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((3+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(3+1):
            dims[i]=cdims[i]
        return dims[1:3+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=3, dtype=np.complex128_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=3, dtype=np.complex128_t] data):
        data -- np.complex128: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<cuDoubleComplex*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<cuDoubleComplex*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.complex128) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex128_t, ndim=3] res=np.zeros(sh,dtype=np.complex128)
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex128)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<cuDoubleComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.complex128) of 3 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex128_t, ndim=3] res=np.zeros(sh, dtype=np.complex128)
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex128)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<cuDoubleComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_ComplexD3D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_ComplexD3D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_ComplexD3D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_ComplexD3D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_ComplexD3D src):
        """Cublas dot

        src -- naga_obj_ComplexD3D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.complex128_t alpha):
        """Cublas scal

        alpha -- np.complex128: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuDoubleComplex a
        a.x=alpha.x
        a.y=alpha.y

        self.c_o.scale(a,1)


    def swap(self, naga_obj_ComplexD3D src):
        """Cublas swap

        src -- naga_obj_ComplexD3D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_ComplexD3D src):
        """Cublas copy

        src -- naga_obj_ComplexD3D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)



    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng(seed)
        self.c_o.prng(b'U')




    def fft(self, bool inplace=False, bool C2R=False, dest=None, int direction=1):
        """Compute fft, using "cufftExec"

        dest -- naga_obj (default = None)
        dir -- integer (default 1)

        dir: fft's direction
        if dest is None, inplace fft (only available for C2C fft)

        Return dest= fft(self,dir)
        """

        cdef plan

        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd=self.get_Dims()

        cdef int i

        cdef long dims_src[3+1]
        cdef long dims_dest[3+1]
        dims_src[0]=3
        dims_dest[0]=3
        for i in range(3):
            dims_src[i+1]=ds[i]
            dims_dest[i+1]=ds[i]

        if(C2R):
            dims_dest[3]=(dims_dest[3]-1)*2
            dd[3-1]=(dd[3-1]-1)*2

        if(dest is None):
            if(inplace):
                dest=self
            else:
                if(C2R):
                    dest=naga_obj_Double3D(self.getContext(),dims=dd)
                else:
                    dest=naga_obj_ComplexD3D(self.getContext(),dims=dd)

        cdef t_dest=dest.__class__
        try:
            plan=dict_fft_ComplexD3D[t_dest]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return
        output = dest.getData_ptr()

        if(t_dest==self.__class__):
            carma_initfft[cuDoubleComplex,cuDoubleComplex](dims_src,self.c_o.getPlan(),plan)
            carma_fft[cuDoubleComplex,cuDoubleComplex](self.c_o.getData(), <cuDoubleComplex*><void*>output, direction,handle[0])
        else:
            carma_initfft[cuDoubleComplex,double](dims_dest,self.c_o.getPlan(),plan)
            carma_fft[cuDoubleComplex,double](self.c_o.getData(), <double*><void*>output, direction,handle[0])

        if(not inplace):
            return dest


##########################################################
#  P-Class `ComplexD4D`
#########################################################


cdef class naga_obj_ComplexD4D:

    def __cinit__(self,naga_context ctxt=None ,
                  naga_obj_ComplexD4D obj=None,
                  np.ndarray[ndim=4, dtype=np.complex128_t] data=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] dims=None,
                  ):
        """naga_obj_ComplexD4D constructor.

        naga_obj(naga_obj obj)
        naga_obj(naga_context ctxt, naga_obj obj)
        naga_obj(naga_context ctxt, np.ndarray dims)
        naga_obj(naga_context ctxt, np.ndarray data)

        input data must be a np.ndarray with ndim=4
        """

        if (ctxt is None):
            #if no context: copy constructor
            #using the context associated to the naga_obj argument
            if(obj is None):
                raise ValueError("Missing argument naga_context or naga_obj")
            else:
                obj.activateDevice()
                self.c_o=new carma_obj[cuDoubleComplex](obj.c_o)
                self.context=obj.getContext()
            return

        self.context=ctxt

        if(obj is not None):
            #copy constructor
            #using context given in argument
            ctxt.set_activeDevice(obj.getDevice())
            self.c_o = new carma_obj[cuDoubleComplex](ctxt.c,obj.c_o)
            return

        cdef long cdims[4+1]
        cdims[0]=4
        cdef int i

        if(data is None):
            if(dims is None):
                raise ValueError("Missing argument data (np.array[ndim=4]) or dims (np.ndarray[ndim=1])")

            if( dims.shape[0]!=4 ):
                raise ValueError("Wrong number of dimension: got ",dims.shape[0]," (expected: 4)")
            for i in range(4):
                cdims[i+1]=dims[i]
            #self.c_o= new carma_obj[cuDoubleComplex](ctxt.c,<long*>dims.data)
            self.c_o= new carma_obj[cuDoubleComplex](ctxt.c,cdims)
            return

        cdef np.ndarray data_F=np.asfortranarray(data)
        for i in range(4):
            cdims[i+1]=data.shape[i]
        if(data.flags.fortran):
            self.c_o=new carma_obj[cuDoubleComplex](ctxt.c,cdims,<cuDoubleComplex*>data.data)
        else:
            self.c_o=new carma_obj[cuDoubleComplex](ctxt.c,cdims,<cuDoubleComplex*>data_F.data)

    def __dealloc__(self):
        del self.c_o

    def getData_ptr(self):
        """Return the pointer value to the naga data (as an integer, type:uintptr_t)."""
        cdef uintptr_t data_ptr=<uintptr_t>self.c_o.getData()
        return data_ptr


    def getCarma_ptr(self):
        """Return the pointer value to the carma object of the naga (as an integer, type:uintptr_t)."""
        cdef uintptr_t carma_ptr=<uintptr_t>self.c_o
        return carma_ptr


    def get_Dims(self):
        """Return the dimensions of the naga_obj."""

        cdef np.ndarray dims = np.ndarray((4+1),dtype=np.int64)
        cdef const long *cdims
        cdef int i
        cdims=self.c_o.getDims()
        for i in range(4+1):
            dims[i]=cdims[i]
        return dims[1:4+1]


    def getNbElem(self):
        """Return the number of elements of the naga object."""
        return self.c_o.getNbElem()


    def getDevice(self):
        """Return the device used by the current naga_obj."""
        return self.c_o.getDevice()


    def activateDevice(self, int silent=1):
        """Activate the device used by the current naga_obj."""
        cdef int dev = self.c_o.getDevice()
        self.context.set_activeDevice(dev,silent)


    def getContext(self):
        """Return a pointer to the carma_context associated with the current naga_obj."""
        return self.context


    def is_rng_init(self):
        self.c_o.is_rng_init()


    def host2device(self,np.ndarray[ndim=4, dtype=np.complex128_t] data):
        """Copy data from host to device.

        host2device(np.ndarray[ndim=4, dtype=np.complex128_t] data):
        data -- np.complex128: data to copy from host to device
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        if(data.flags.fortran):
            self.c_o.host2device(<cuDoubleComplex*>data.data)
            return
        cdef np.ndarray data_F=data.flatten(order="F")
        self.c_o.host2device(<cuDoubleComplex*>data_F.data)


    def device2host(self):
        """Copy data from device to host.

        return np.ndarray(dtype=np.complex128) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex128_t, ndim=4] res=np.zeros(sh,dtype=np.complex128)
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex128)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2host(<cuDoubleComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res

    def device2hostOpt(self):
        """Copy data from device o_data to host.

        return np.ndarray(dtype=np.complex128) of 4 dimension(s)
        """
        sh=self.get_Dims()
        cdef np.ndarray[ dtype=np.complex128_t, ndim=4] res=np.zeros(sh, dtype=np.complex128)
        cdef np.ndarray[ dtype=np.complex128_t, ndim=1] tmp=np.zeros(self.c_o.getNbElem(),
                                                                        dtype=np.complex128)
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.device2hostOpt(<cuDoubleComplex*>tmp.data)
        res=tmp.reshape(sh,order="F")
        return res


    def copyFrom(self, naga_obj_ComplexD4D src):
        """Copy the data from src to the current naga_obj.

        src -- naga_obj_ComplexD4D: object to copy the data from.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyFrom( src.c_o.getData(), src.c_o.getNbElem())


    def copyInto(self, naga_obj_ComplexD4D dest):
        """Copy data from current naga_obj to dest.

        dest -- naga_obj_ComplexD4D: object to copy the data into.
        """
        if(self.getDevice()!=dest.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copyInto(dest.c_o.getData(), dest.c_o.getNbElem())

    def reset(self):
        """Set naga_obj to zero.

        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.reset()


    def getValues(self):
        return <np.int64_t> self.c_o.getValues()



    def imax (self):
        """Cublas amax.

        Return the smallest index of the maximum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imax(1)


    def imin (self):
        """Cublas amin

        Return the smallest index of the minimum absolute magnitude element."""
        self.context.set_activeDevice(self.getDevice())
        return self.c_o.imin(1)


    def dot(self,naga_obj_ComplexD4D src):
        """Cublas dot

        src -- naga_obj_ComplexD4D
        return the dot product of src and self.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("dot only on the same device")
        return self.c_o.dot(src.c_o,1,1)


    def scale(self, np.complex128_t alpha):
        """Cublas scal

        alpha -- np.complex128: caling factor
        self = alpha.self
        """
        self.context.set_activeDeviceForCpy(self.getDevice())
        cdef cuDoubleComplex a
        a.x=alpha.x
        a.y=alpha.y

        self.c_o.scale(a,1)


    def swap(self, naga_obj_ComplexD4D src):
        """Cublas swap

        src -- naga_obj_ComplexD4D
        Swap data contents of naga objects self and src.
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("swap only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.swap(src.c_o,1,1)


    def copy(self, naga_obj_ComplexD4D src):
        """Cublas copy

        src -- naga_obj_ComplexD4D
        Copy data from src into self
        """
        if(self.getDevice()!=src.getDevice()):
            raise ValueError("copy only on the same device")
        self.context.set_activeDeviceForCpy(self.getDevice())
        self.c_o.copy(src.c_o,1,1)



    def random(self,int seed=1234):
        """Generate random values for this naga datas.

        seed -- integer: seed for random function (default:1234)
        """
        if(not self.c_o.is_rng_init()):
            self.c_o.init_prng(seed)
        self.c_o.prng(b'U')




    def fft(self, bool inplace=False, bool C2R=False, dest=None, int direction=1):
        """Compute fft, using "cufftExec"

        dest -- naga_obj (default = None)
        dir -- integer (default 1)

        dir: fft's direction
        if dest is None, inplace fft (only available for C2C fft)

        Return dest= fft(self,dir)
        """

        cdef plan

        cdef cufftHandle *handle = self.c_o.getPlan()
        cdef uintptr_t output

        cdef np.ndarray ds=self.get_Dims()
        cdef np.ndarray dd=self.get_Dims()

        cdef int i

        cdef long dims_src[4+1]
        cdef long dims_dest[4+1]
        dims_src[0]=4
        dims_dest[0]=4
        for i in range(4):
            dims_src[i+1]=ds[i]
            dims_dest[i+1]=ds[i]

        if(C2R):
            dims_dest[4]=(dims_dest[4]-1)*2
            dd[4-1]=(dd[4-1]-1)*2

        if(dest is None):
            if(inplace):
                dest=self
            else:
                if(C2R):
                    dest=naga_obj_Double4D(self.getContext(),dims=dd)
                else:
                    dest=naga_obj_ComplexD4D(self.getContext(),dims=dd)

        cdef t_dest=dest.__class__
        try:
            plan=dict_fft_ComplexD4D[t_dest]
        except KeyError:
            print("fft not available from ", self.__class__.__name__, " to ",dest.__class__.__name__)
            return
        output = dest.getData_ptr()

        if(t_dest==self.__class__):
            carma_initfft[cuDoubleComplex,cuDoubleComplex](dims_src,self.c_o.getPlan(),plan)
            carma_fft[cuDoubleComplex,cuDoubleComplex](self.c_o.getData(), <cuDoubleComplex*><void*>output, direction,handle[0])
        else:
            carma_initfft[cuDoubleComplex,double](dims_dest,self.c_o.getPlan(),plan)
            carma_fft[cuDoubleComplex,double](self.c_o.getData(), <double*><void*>output, direction,handle[0])

        if(not inplace):
            return dest


