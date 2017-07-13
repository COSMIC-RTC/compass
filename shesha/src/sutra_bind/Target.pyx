import numpy as np
cimport numpy as np
# np.import_array()

include "../par.pxi"

#################################################
# P-Class target
#################################################
cdef class Target:

    def __cinit__(self, naga_context ctxt, Telescope telescope, int ntargets,
                  np.ndarray[ndim=1, dtype=np.float32_t] xpos,
                  np.ndarray[ndim=1, dtype=np.float32_t] ypos,
                  np.ndarray[ndim=1, dtype=np.float32_t] Lambda,
                  np.ndarray[ndim=1, dtype=np.float32_t] mag,
                  float zerop,
                  np.ndarray[ndim=1, dtype=np.int64_t] size,
                  int Npts,
                  int device=-1
                  ):

        if(device < 0):
            device = ctxt.get_activeDevice()
        else:
            ctxt.set_activedevice(device)

        self.ntargets = ntargets
        self.xpos = xpos
        self.ypos = ypos
        self.Lambda = Lambda
        self.mag = mag
        self.device = device
        self.context = ctxt

        self.target = new sutra_target(ctxt.c, telescope.telescope, ntargets,
                                       < float * > xpos.data, < float * > ypos.data,
                                       < float * > Lambda.data, < float * > mag.data,
                                       zerop, < long * > size.data,
                                       Npts, device)

    def __dealloc__(self):
        del self.target

    def add_layer(self, int n, bytes l_type, float alt, float xoff, float yoff):
        """
            Add a phase screen dm or atmos as layers of turbulence

        :parameters:
            n: (int) : index of the target
            l_type: (str) : "atmos" or "dm"
            alt: (float) : altitude
            xoff: (float) : x-offset
            yoff: (float) : y-offset
        """
        self.context.set_activeDevice(self.device)
        self.target.d_targets[n].add_layer( < char * > l_type, alt, xoff, yoff)

    def init_strehlmeter(self, int nTarget):
        """
            Initialise target's strehl

        :param nTarget: (int) : index of the target
        """
        cdef sutra_source * s_s_ptr
        self.target.d_targets[nTarget].init_strehlmeter()

    def atmos_trace(self, int nTarget, Atmos atm, Telescope tel):
        """
            Raytracing of the target through the atmosphere

        :parameters:
            nTarget: (int)   : index of the target
            atm: (atmos)     : atmos to get through
            tel: (Telescope) : telescope
        """

        self.context.set_activeDevice(self.device)
        self.target.d_targets[nTarget].raytrace(atm.s_a)
        cdef carma_obj[float] * d_screen = self.target.d_targets[nTarget].d_phase.d_screen
        cdef carma_obj[float] * d_tel = tel.telescope.d_phase_ab_M1

        d_screen.axpy(1.0, d_tel, 1, 1)

    def get_image(self, int nTarget, bytes type_im, long puponly=0, bool comp_le=False, bool fluxNorm=True):
        """
            Return the image from the target (or long exposure image according to the requested type)

        :parameters:
            nTarget: (int) : index of the target
            type_im: (str) : type of the image to get ("se" or "le")
            puponly: (int) : if 1, image computed from phase on the pupil only
            comp_le: (bool) : if False (default), the computed image is not taken into account in the LE image
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[nTarget]

        cdef const long * dims = src.d_image.getDims()
        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.float32)
        cdef float flux
        cdef carma_obj[float] * tmp_img
        if(fluxNorm):
            src.comp_image(puponly, comp_le)
            flux = src.zp * 10 ** (-0.4 * src.mag)
            tmp_img = new carma_obj[float](self.context.c, dims)
            if(type_im == b"se"):
                roll_mult[float](
                    tmp_img.getData(),
                    src.d_image.getData(), src.d_image.getDims(1), src.d_image.getDims(2),
                    flux,
                    self.context.c.get_device(src.device))
            elif(type_im == b"le"):
                roll_mult[float](
                    tmp_img.getData(),
                    src.d_leimage.getData(), src.d_leimage.getDims(1), src.d_leimage.getDims(2),
                    flux,
                    self.context.c.get_device(src.device))

            tmp_img.device2host( < float * > data_F.data)
            del tmp_img
        else:
            if(type_im == b"se"):
                src.d_image.device2host( < float*> data_F.data)
            if(type_im == b"le"):
                src.d_leimage.device2host( < float*> data_F.data)
                data_F /= src.strehl_counter

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    cpdef set_ncpa_phase(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data_F = data.flatten("F")
        self.target.d_targets[n].set_ncpa_phase( < float*>data_F.data, data.size)

    cpdef get_ncpa_phase(self, int n):
        cdef carma_obj[float] * ph
        cdef const long * cdims
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data

        ph = self.target.d_targets[n].d_phase.d_screen
        cdims = ph.getDims()
        data = np.empty((cdims[2], cdims[1]), dtype=np.float32)
        self.target.d_targets[n].get_ncpa_phase( < float*>data.data, data.size)
        return np.reshape(data.flatten("F"), (cdims[1], cdims[2]))

    def reset_phase(self, int nTarget):
        """
            Reset the phase's screen of the target

        :param nTarget: (int) : index of the target
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[nTarget]
        src.d_phase.d_screen.reset()

    def get_phase(self, int nTarget):
        """
            Return the phase's screen of the target

        :param nTarget: (int) : index of the target
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[nTarget]

        cdef const long * dims
        dims = src.d_phase.d_screen.getDims()
        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.float32)
        src.d_phase.d_screen.device2host( < float * > data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def set_phase(self, int nTarget, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            Set the phase's screen of the target

        :param nTarget: (int) : index of the target
        :param data: (np.ndarray[ndim=2,dtype=np.float32_t]) : phase screen
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[nTarget]

        cdef np.ndarray[dtype = np.float32_t] data_F = data.flatten("F")

        src.d_phase.d_screen.host2device( < float * > data_F.data)

    def set_pupil(self, int nTarget, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            Set the pupil used for PSF computation

        :param nTarget: (int) : index of the target
        :param data: (np.ndarray[ndim=2,dtype=np.float32_t]) : pupil
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[nTarget]
        cdef np.ndarray[dtype = np.float32_t] data_F
        cdef np.ndarray[dtype = np.int32_t] wherephase
        cdef long Npts
        cdef const long * dims
        cdef long dims1[2]
        dims = src.d_pupil.getDims()
        if((dims[2], dims[1]) == np.shape(data)):
            data_F = data.flatten("F")
            src.d_pupil = new carma_obj[float](src.current_context, dims)
            src.d_pupil.host2device(< float * > data_F.data)
            wherephase = np.where(data_F)[0].astype(np.int32)
            Npts = wherephase.size
            dims1[0] = 1
            dims1[1] = Npts
            del src.d_phasepts
            del src.d_wherephase
            src.d_phasepts = new carma_obj[float](src.current_context, dims1)
            src.d_wherephase = new carma_obj[int](src.current_context, dims1)
            src.d_wherephase.host2device(< int * > wherephase.data)
        else:
            raise IndexError("Pupil dimension mismatch")

    def get_phasetele(self, int nTarget):
        """
            Return the telemetry phase of the target

        :param nTarget: (int) : index of the target
        :return data: (np.ndarray(ndim=2,np.float32)) : phase screen
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[nTarget]

        cdef const long * dims
        dims = src.phase_telemetry.getDims()
        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.float32)

        src.phase_telemetry.fill_into( < float * > data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def get_amplipup(self, int nTarget):
        """
            Return the complex amplitude in the pupil plane of the target.

        :param nTarget: (int) : index of the target
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[nTarget]

        cdef const long * dims
        dims = src.d_amplipup.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.complex64)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.complex64)
        src.d_amplipup.device2host( < cuFloatComplex * > data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def reset_strehl(self, int nTarget):
        """
            Reset the target's strehl

        :param nTarget: (int) : index of the target
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[nTarget]
        src.init_strehlmeter()

    def get_strehl(self, int nTarget, bool comp_strehl=True):
        """
            Compute and return the target's strehl

        :param nTarget: (int) : index of the target
        :return strehl: (np.array(4,dtype=np.float32)) :
            [Strehl SE, Strehl LE, instantaneous phase variance over the pupil,
                average  phase variance over the pupil]
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[nTarget]

        cdef np.ndarray strehl = np.empty(4, dtype=np.float32)
        if(comp_strehl):
            src.comp_image(0, True)
            src.comp_strehl()

        strehl[0] = src.strehl_se
        strehl[1] = src.strehl_le
        strehl[2] = src.phase_var
        if(src.phase_var_count > 0):
            strehl[3] = src.phase_var_avg / float(src.phase_var_count)
        else:
            strehl[3] = 0.

        return strehl

    def dmtrace(self, int ntar, Dms dms, int reset=0, int do_phase_var=1):
        """
            Raytracing of the target through the dms

        :parameters:
            ntar: (int)   : index of the target
            dms: (Dms)     : dms to go through
            reset: (int) : if >0, reset the screen before raytracing
            do_phase_var: (int) : if 0, doesn't take the screen into account in the phase average (unused)
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.target.d_targets[ntar].device, 1)

        self.target.d_targets[ntar].raytrace(dms.dms, reset, do_phase_var)

    def ncpatrace(self, int ntar, int reset=0):
        """
            Raytracing of the target through NCPA

        :parameters:
            ntar: (int)   : index of the target
            reset: (int) : if >0, reset the screen before raytracing
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.target.d_targets[ntar].device, 1)

        self.target.d_targets[ntar].raytrace(reset)

    def __str__(self):
        info = "Target objects:\n"
        info += "Contains " + str(self.target.ntargets) + " Target(s):\n"
        cdef int i
        cdef sutra_source * target
        info += "Source # | position(\") |  Mag  | Lambda (mic.)\n"
        for i in range(self.target.ntargets):
            target = self.target.d_targets.at(i)
            info += "%8d" % (i + 1) + " | " + "%4d" % target.tposx + " , " + "%-4d" % target.tposy + \
                " | " + "%5.2f" % target.mag + " | " + "%5.3f" % target.Lambda + "\n"
        info += "--------------------------------------------------------"
        return info


IF USE_BRAMA == 1:
    # child constructor must have the same prototype (same number of
    # non-optional arguments)
    cdef class Target_brama(Target):
        def __cinit__(self, naga_context ctxt, Telescope telescope, int ntargets,
                      np.ndarray[ndim=1, dtype=np.float32_t] xpos,
                      np.ndarray[ndim=1, dtype=np.float32_t] ypos,
                      np.ndarray[ndim=1, dtype=np.float32_t] Lambda,
                      np.ndarray[ndim=1, dtype=np.float32_t] mag,
                      float zerop,
                      np.ndarray[ndim=1, dtype=np.int64_t] size,
                      int Npts,
                      int device=-1
                      ):
            del self.target

            cdef carma_context * context = &carma_context.instance()
            self.target = new sutra_target_brama(context, "target_brama", telescope.telescope, -1, ntargets,
                                                 < float * > xpos.data, < float * > ypos.data,
                                                 < float * > Lambda.data, < float * > mag.data,
                                                 zerop, < long * > size.data,
                                                 Npts, context.get_activeDevice())

        def __dealloc__(self):
            pass  # del self.target

        cpdef publish(self):
            cdef sutra_target_brama * target = < sutra_target_brama * > (self.target)
            target.publish()
