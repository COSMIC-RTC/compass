import numpy as np
cimport numpy as np
np.import_array()

def roket_init(Rtc rtc, Sensors sensors, Target target, Dms dms, Telescope tel,
                  Atmos atm, int loopcontroller, int geocontroller, int nactus, int nmodes,
                  int nfilt, int niter,
                  np.ndarray[ndim=2, dtype=np.float32_t] Btt,
                  np.ndarray[ndim=2, dtype=np.float32_t] P,
                  np.ndarray[ndim=2, dtype=np.float32_t] gRD,
                  np.ndarray[ndim=2, dtype=np.float32_t] RD):

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
        cdef carma_context * context = &carma_context.instance()
        cdef int device
        device = 1
        cdef np.ndarray[dtype = np.float32_t] Btt_F = Btt.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] P_F = P.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] gRD_F = gRD.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] RD_F = RD.flatten("F")

        self.roket = new sutra_roket(context, device, rtc.rtc, sensors.sensors,
                                      target.target,dms.dms,tel.telescope, atm.s_a,
                                      loopcontroller,geocontroller,nactus,nmodes,
                                      nfilt,niter,< float *> Btt_F.data,
                                      < float *> P_F.data, < float *> gRD_F.data,
                                      < float *> RD_F.data)


    def __dealloc__(self):
        del self.roket
        
    cpdef getContributor(self, bytes type):
        """Return the error buffer define by "type" from a sutra_roket object.
        :parameters:
            type: (str) : error buffer needed ("noise", "nonlinear", "tomo",
                            "filtered", "aliasing", "bandwidth" or "fitting")
        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : error buffer
                    (just a float value if type == "fitting")
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.roket.device, 1)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef carma_obj[float] * buffer
        cdef const long * dims = NULL

        if(type == "fitting"):
            return self.roket.fitting

        if(type == "noise"):
            buffer = self.roket.d_noise
        elif(type == "nonlinear"):
            buffer = self.roket.d_nonlinear
        elif(type == "tomo"):
            buffer = self.roket.d_tomo
        elif(type == "filtered"):
            buffer = self.roket.d_filtered
        elif(type == "aliasing"):
            buffer = self.roket.d_alias
        elif(type == "bandwidth"):
            buffer = self.roket.d_bandwidth

        dims = buffer.getDims()
        data_F = np.zeros((dims[1], dims[2]), dtype=np.float32)
        buffer.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data_F

    cpdef get_Btt(self):
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
        self.roket.d_Btt.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    cpdef get_RD(self):
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
        self.roket.d_RD.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    cpdef get_gRD(self):
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
        self.roket.d_gRD.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data
    cpdef get_P(self):
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
        self.roket.d_P.device2host(< float *> data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

        return data

    cpdef computeBreakdown(self):
        """ Launch the computation of the error breakdown for the current iteration
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.roket.device, 1)
        self.roket.compute_breakdown()
