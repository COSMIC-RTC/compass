import numpy as np
cimport numpy as np
# np.import_array()


def roket_init(Rtc rtc, Sensors sensors, Target target, Dms dms, Telescope tel,
               Atmos atm, int loopcontroller, int geocontroller, int nactus, int nmodes,
               int nfilt, int niter,
               np.ndarray[ndim=2, dtype=np.float32_t] Btt,
               np.ndarray[ndim=2, dtype=np.float32_t] P,
               np.ndarray[ndim=2, dtype=np.float32_t] gRD,
               np.ndarray[ndim=2, dtype=np.float32_t] RD):
    """
    Initialize a Roket object

    :parameters:
        rtc: (Rtc): Rtc object
        sensors: (Sensors): Sensors object
        target: (Target) : Target object
        dms: (Dms) : Dms object
        tel: (Telescope) : Telescope object
        atm: (Atmos) : Atmos object
        loopcontroller: (int) : index of the loop controller
        geocontroller: (int) : index of the geo controller
        nactus: (int) : total number of actuators
        nmodes: (int) : total number of modes
        nfilt: (int) : number of modes to filter
        niter: (int): number of iterations
        Btt: (np.ndarray[ndim=2, dtype=np.float32_t]) : Volts to Btt modes matrix
        P: (np.ndarray[ndim=2, dtype=np.float32_t]) : Btt modes to Volts matrix
        gRD: (np.ndarray[ndim=2, dtype=np.float32_t]) : (1-g*cmat*imat) matrix
        RD: (np.ndarray[ndim=2, dtype=np.float32_t]) cmat * imat matrix
    """

    return Roket(rtc, sensors, target, dms, tel, atm, loopcontroller, geocontroller,
                 nactus, nmodes, nfilt, niter, Btt, P, gRD, RD)


cdef class Roket:

    def __cinit__(self, Rtc rtc, Sensors sensors, Target target, Dms dms, Telescope tel,
                  Atmos atm, int loopcontroller, int geocontroller, int nactus, int nmodes,
                  int nfilt, int niter,
                  np.ndarray[ndim=2, dtype=np.float32_t] Btt,
                  np.ndarray[ndim=2, dtype=np.float32_t] P,
                  np.ndarray[ndim=2, dtype=np.float32_t] gRD,
                  np.ndarray[ndim=2, dtype=np.float32_t] RD):
        """
        Initialize a sutra_roket object on the GPU

        :parameters:
          rtc: (Rtc): Rtc object
          sensors: (Sensors): Sensors object
          target: (Target) : Target object
          dms: (Dms) : Dms object
          tel: (Telescope) : Telescope object
          atm: (Atmos) : Atmos object
          loopcontroller: (int) : index of the loop controller
          geocontroller: (int) : index of the geo controller
          nactus: (int) : total number of actuators
          nmodes: (int) : total number of modes
          nfilt: (int) : number of modes to filter
          niter: (int): number of iterations
          Btt: (np.ndarray[ndim=2, dtype=np.float32_t]) : Volts to Btt modes matrix
          P: (np.ndarray[ndim=2, dtype=np.float32_t]) : Btt modes to Volts matrix
          gRD: (np.ndarray[ndim=2, dtype=np.float32_t]) : (1-g*cmat*imat) matrix
          RD: (np.ndarray[ndim=2, dtype=np.float32_t]) cmat * imat matrix
        """
        cdef carma_context * context = &carma_context.instance()
        cdef int device
        device = 1
        cdef np.ndarray[dtype = np.float32_t] Btt_F = Btt.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] P_F = P.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] gRD_F = gRD.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] RD_F = RD.flatten("F")

        self.roket = new sutra_roket(context, device, rtc.rtc, sensors.sensors,
                                     target.target, dms.dms, tel.telescope, atm.s_a,
                                     loopcontroller, geocontroller, nactus, nmodes,
                                     nfilt, niter, < float * > Btt_F.data,
                                     < float * > P_F.data, < float * > gRD_F.data,
                                     < float * > RD_F.data)

    def __dealloc__(self):
        del self.roket

    def getContributor(self, bytes type):
        """Return the error buffer define by "type" from a sutra_roket object.
        :parameters:
            type: (str) : error buffer needed ("noise", "nonlinear", "tomo",
                            "filtered", "aliasing", "bandwidth" or "fitting")
        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : error buffer
                    (just a float value if type == b"fitting")
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef carma_obj[float] * buffer
        cdef const long * dims = NULL

        if(type == b"fitting"):
            return self.roket.fitting

        if(type == b"noise"):
            buffer = self.roket.d_noise
        elif(type == b"nonlinear"):
            buffer = self.roket.d_nonlinear
        elif(type == b"tomo"):
            buffer = self.roket.d_tomo
        elif(type == b"filtered"):
            buffer = self.roket.d_filtered
        elif(type == b"aliasing"):
            buffer = self.roket.d_alias
        elif(type == b"bandwidth"):
            buffer = self.roket.d_bandwidth
        else:
            raise "type unknown"

        dims = buffer.getDims()
        data_F = np.zeros((dims[1], dims[2]), dtype=np.float32)
        buffer.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data_F

    def get_Btt(self):
        """Return the Btt matrix from a sutra_roket object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : Btt
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.roket.d_Btt.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.roket.d_Btt.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_covv(self):
        """Return the command covariance matrix from a sutra_roket object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : covv
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.roket.d_covv.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.roket.d_covv.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_covm(self):
        """Return the slopes covariance matrix from a sutra_roket object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : Btt
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.roket.d_covm.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.roket.d_covm.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_psfortho(self):
        """Return the PSF orthogonal from a sutra_roket object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : psf
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.roket.d_psfortho.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.roket.d_psfortho.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_filtmodes(self):
        """Return the filt modes from a sutra_roket object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : filtmodes
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.roket.d_filtmodes.getDims()
        data = np.zeros(dims[1], dtype=np.float32)
        self.roket.d_filtmodes.device2host( < float * > data.data)

        return data

    def get_modes(self):
        """Return the modes from a sutra_roket object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) :modes
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.roket.d_modes.getDims()
        data = np.zeros(dims[1], dtype=np.float32)
        self.roket.d_modes.device2host( < float * > data.data)

        return data

    def get_RD(self):
        """Return the RD matrix from a sutra_roket object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : RD
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.roket.d_RD.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.roket.d_RD.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_gRD(self):
        """Return the gRD matrix from a sutra_roket object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : gRD
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.roket.d_gRD.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.roket.d_gRD.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def get_P(self):
        """Return the P matrix from a sutra_roket object.

        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : P
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        dims = self.roket.d_P.getDims()
        data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.roket.d_P.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    def computeBreakdown(self):
        """ Launch the computation of the error breakdown for the current iteration
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.roket.device, 1)
        self.roket.compute_breakdown()
