import numpy as np
cimport numpy as np
np.import_array()

def psfrecs_init( bytes type, int nactus, int niter,
                  np.ndarray[ndim=1, dtype=np.float32_t] IFvalue,
                  np.ndarray[ndim=1, dtype=np.int32_t] IFrowind,
                  np.ndarray[ndim=1, dtype=np.int32_t] IFcolind,
                  np.ndarray[ndim=2, dtype=np.float32_t] TT,
                  np.ndarray[ndim=2, dtype=np.float32_t] spupil,
                  float scale, int nmodes=2,
                  np.ndarray[ndim=2, dtype=np.float32_t] Btt=None,
                  np.ndarray[ndim=2, dtype=np.float32_t] covmodes=None):
    """ Initialize a shesha_psfrecs object for psf reconstructions
    :parameters:
        type : (str) : reconstruction method used ("roket" or "Vii")
        nactus : (int) : number of actuators
        niter : (int) : number of iterations performed with roket
        IFvalue : (np.ndarray[ndim=1,dtype=float32_t]) : Non zeros values of pzt influence function matrix
        IFrowind : (np.ndarray[ndim=1,dtype=int32_t]) : Row indices of nnz values (csr sparse format)
        IFcolind : (np.ndarray[ndim=1,dtype=int32_t]) : Column indices of nnz values (csr sparse format)
        TT : (np.ndarray[ndim=1,dtype=float32_t])np.ndarray[ndim=1,dtype=float32_t]) : Tip-tilt influence functions
        spupil : (np.ndarray[ndim=2,dtype=float32_t]) : Small pupil
        scale : (float) : 2*pi/lambda_target with lambda_target expressed in microns
    """
    if Btt is None:
        return Psfrecs(type,nactus, niter, IFvalue, IFrowind, IFcolind, TT, spupil, scale)
    else:
        return Psfrecs(type,nactus, niter, IFvalue, IFrowind, IFcolind, TT, spupil, scale, nmodes, Btt, covmodes)

cdef class Psfrecs:

    def __cinit__(self, bytes type, int nactus, int niter,
                  np.ndarray[ndim=1, dtype=np.float32_t] IFvalue,
                  np.ndarray[ndim=1, dtype=np.int32_t] IFrowind,
                  np.ndarray[ndim=1, dtype=np.int32_t] IFcolind,
                  np.ndarray[ndim=2, dtype=np.float32_t] TT,
                  np.ndarray[ndim=2, dtype=np.float32_t] spupil,
                  float scale, int nmodes=2,
                  np.ndarray[ndim=2, dtype=np.float32_t] Btt=None,
                  np.ndarray[ndim=2, dtype=np.float32_t] covmodes=None):

        cdef carma_context * context = &carma_context.instance()
        cdef int device
        device = context.get_activeDevice()
        self.device = device
        cdef int size = spupil.shape[0]
        cdef int Npts = np.where(spupil)[0].shape[0]
        cdef int IFnz = IFvalue.shape[0]
        cdef np.ndarray[dtype = np.float32_t] pup_F = spupil.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] TT_F = TT.flatten("F")
        if Btt is None:
            Btt = np.zeros((2,2),dtype=np.float32)
            covmodes = np.zeros((2,2),dtype=np.float32)
        cdef np.ndarray[dtype = np.float32_t] Btt_F = Btt.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] cov_F = covmodes.flatten("F")

        self.psfrecs = new sutra_psfrecs(context, device, type, nactus, nmodes,
                                      niter,< float *> IFvalue.data,
                                      < int *> IFrowind.data, < int *> IFcolind.data,
                                      IFvalue.shape[0],<float *> TT_F.data, < float *> pup_F.data,
                                      size, Npts, scale, <float*> Btt_F.data, <float*> cov_F.data)

    def __dealloc__(self):
        del self.psfrecs

    cpdef psf_rec_roket(self,np.ndarray[ndim=2, dtype=np.float32_t] err):
        """ Retrieve the psf from err buffer of command error from a roket output

        :parameter err: (np.ndarray[ndim=2, dtype=np.float32_t]) : buffer of command error
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.psfrecs.device, 1)
        cdef np.ndarray[dtype = np.float32_t] err_F = err.flatten("F")

        self.psfrecs.psf_rec_roket(<float *> err_F.data)

    cpdef psf_rec_Vii(self):
        """ Compute Telescope OTF and residual OTF using Vii algorithm
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.psfrecs.device, 1)

        self.psfrecs.psf_rec_Vii()

    def get_psf(self):
        """Return the reconsructed psf from a sutra_psfrecs object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : psf
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.psfrecs.d_psf.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.psfrecs.d_psf.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return np.fft.fftshift(data)

    def get_otftel(self):
        """Return the telescope OTF computed by a sutra_psfrecs object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : otftel
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        if(self.psfrecs.d_otftel):
            dims = self.psfrecs.d_otftel.getDims()
            data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            self.psfrecs.d_otftel.device2host(< float *> data_F.data)
            data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

            return data
        else:
            raise ValueError("otftel was not initialized")

    def get_mask(self):
        """Return the telescope OTF computed by a sutra_psfrecs object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : otftel
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        if(self.psfrecs.d_mask):
            dims = self.psfrecs.d_otftel.getDims()
            data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            self.psfrecs.d_mask.device2host(< float *> data_F.data)
            data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

            return data
        else:
            raise ValueError("otftel was not initialized")

    def get_otfVii(self):
        """Return the parallel OTF computed by a sutra_psfrecs object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : otfvii
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        if(self.psfrecs.d_otfVii):
            dims = self.psfrecs.d_otfVii.getDims()
            data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            self.psfrecs.d_otfVii.device2host(< float *> data_F.data)
            data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

            return data
        else:
            raise ValueError("otfVii was not initialized")

    def get_err(self):
        """Return the command error buffer from a sutra_psfrecs object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : err
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.psfrecs.d_err.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.psfrecs.d_err.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_phase(self):
        """Return the phase in the pupil computed by a sutra_psfrecs object.

        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : phase
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.psfrecs.d_phase.getDims()
        data = np.zeros(dims[1], dtype=np.float32)
        self.psfrecs.d_phase.device2host(< float *> data.data)

        return data

    def get_wherephase(self):
        """Return the wherephase vector of a sutra_psfrecs object.

        :return:
            data : (np.ndarray[ndim=1,dtype=np.int32_t]) : wherephase
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)
        cdef np.ndarray[ndim = 1, dtype = np.int32_t] data
        cdef const long * dims = NULL

        dims = self.psfrecs.d_wherephase.getDims()
        data = np.zeros(dims[1], dtype=np.int32)
        self.psfrecs.d_wherephase.device2host(< int *> data.data)

        return data

    def get_IF(self):
        """Get the influence functions matrix from a sutra_psfrecs object
        Return a scipy.sparse object which shape is (nactus,Npts in the pupil)

        :return:
            IF : (scipy.sparse) : influence functions matrix
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        sparse = naga_sparse_obj_Float()

        sparse.copy(self.psfrecs.d_IF)

        return sparse.get_sparse()

    def get_amplipup(self):
        """Return the complex amplitude in the pupil plane.
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef const long * dims
        dims = self.psfrecs.d_amplipup.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.complex64)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.complex64)
        self.psfrecs.d_amplipup.device2host(< cuFloatComplex *> data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def get_pupfft(self):
        """Return the pupfft object of a sutra_psfrecs
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef const long * dims
        dims = self.psfrecs.d_pupfft.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.complex64)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.complex64)
        self.psfrecs.d_pupfft.device2host(< cuFloatComplex *> data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def get_newmodek(self):
        """Return the newmodek object of a sutra_psfrecs
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef const long * dims
        dims = self.psfrecs.d_newmodek.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.complex64)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.complex64)
        self.psfrecs.d_newmodek.device2host(< cuFloatComplex *> data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def get_covmodes(self):
        """Return the covmodes object of a sutra_psfrecs
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef const long * dims
        dims = self.psfrecs.d_covmodes.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.float32)
        self.psfrecs.d_covmodes.device2host(< float *> data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def get_TT(self):
        """Return the Tip-Tilt influence functions object of a sutra_psfrecs
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef const long * dims
        dims = self.psfrecs.d_TT.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.float32)
        self.psfrecs.d_TT.device2host(< float *> data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def get_term1(self):
        """Return the term1 object of a sutra_psfrecs
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef const long * dims
        dims = self.psfrecs.d_term1.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.float32)
        self.psfrecs.d_term1.device2host(< float *> data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def get_term2(self):
        """Return the term2 object of a sutra_psfrecs
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef const long * dims
        dims = self.psfrecs.d_term2.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.float32)
        self.psfrecs.d_term2.device2host(< float *> data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def get_eigenvals(self):
        """Return the eigen values of the covariance matrix
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims

        dims = self.psfrecs.h_eigenvals.getDims()
        data = np.zeros((dims[1]), dtype=np.float32)
        self.psfrecs.h_eigenvals.fill_into(< float *> data.data)

        return data

    def get_Dphi(self):
        """Return the Dphi
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef const long * dims
        dims = self.psfrecs.d_Dphi.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.complex64)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.complex64)
        self.psfrecs.d_Dphi.device2host(< cuFloatComplex *> data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def set_covmodes(self, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """Set the covariance matrix
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        cdef np.ndarray[dtype = np.float32_t] data_F = data.flatten("F")
        self.psfrecs.d_covmodes.host2device(< float *> data_F.data)

    def set_eigenvals(self, np.ndarray[ndim=1, dtype=np.float32_t] eigenvals):
        """Set the eigen values
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.psfrecs.device, 1)

        self.psfrecs.h_eigenvals.fill_from(< float *> eigenvals.data)
