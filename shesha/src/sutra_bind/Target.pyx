include "../../par.pxi"
import numpy as np
cimport numpy as np

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

        self.device = device
        self.context = ctxt

        self.target = new sutra_target(ctxt.c, telescope.telescope, ntargets,
                                       < float * > xpos.data, < float * > ypos.data,
                                       < float * > Lambda.data, < float * > mag.data,
                                       zerop, < long * > size.data,
                                       Npts, device)

    def __dealloc__(self):
        del self.target

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

    def add_layer(self, int n, bytes l_type, float alt, float xoff, float yoff):
        """Add a phase screen dm or atmos as layers of turbulence

        :parameters:
            n: (int) : index of the target

            l_type: (str) : "atmos" or "dm"

            alt: (float) : altitude

            xoff: (float) : x-offset

            yoff: (float) : y-offset
        """
        self.context.set_activeDevice(self.device)
        self.target.d_targets[n].add_layer(< char * > l_type, alt, xoff, yoff)

    def init_strehlmeter(self, int n):
        """Initialise target's strehl

        :param n: (int) : index of the target
        """
        self.target.d_targets[n].init_strehlmeter()

    def raytrace(self, int n, bytes type_trace, Telescope tel=None, Atmos atmos=None, Dms dms=None, bool rst=False, bool ncpa=False, int do_phase_var=1):
        """
            Does the raytracing for the target phase screen in sutra_target

        :parameters:
            n: (int) :
            type_trace: (str) : "all" : raytracing across atmos and dms seen
                                "dm"  : raytracing across dms seen only
                                "atmos" : raytracing across atmos only
                                "ncpa" : raytracing across ncpa
            tel: (Telescope) :(optional) Telescope object
            atmos: (Atmos) :(optional) Atmos object
            dms: (Dms) : (optional) Dms object
            rst: (int) : (optional) reset before raytracing if rst = 1
            ncpa: (int) : (optional) NCPA raytracing if ncpa = 1
            do_phase_var: (int) : if 0, doesn't take the screen into account in the phase average (unused)
        """
        self.context.set_activeDevice(self.device)
        cdef carma_obj[float] * d_screen = self.target.d_targets[n].d_phase.d_screen
        cdef sutra_atmos * atm = atmos.s_a

        if(type_trace == b"all"):
            rst = 0

        if(type_trace == b"all" or type_trace == b"atmos"):
            self.target.d_targets[n].raytrace(atm)

        if(type_trace == b"all" or type_trace == b"dm"):
            self.target.d_targets[n].raytrace(dms.dms, rst, do_phase_var)

        if(ncpa):
            self.target.d_targets[n].raytrace(rst)

        if tel is not None:
            d_screen.axpy(1.0, tel.telescope.d_phase_ab_M1_m, 1, 1)

    def reset_strehl(self, int n):
        """Reset the target's strehl

        :param n: (int) : index of the target
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]
        src.init_strehlmeter()

    def reset_phase(self, int n):
        """Reset the phase's screen of the target

        :param n: (int) : index of the target
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]
        src.d_phase.d_screen.reset()

    def comp_image(self, int n, bool puponly=False, bool comp_le=True, bool comp_strehl=True):
        """
        Compute the target image (ie. PSF) short and long exposure
        It also compute the Strehl ratio of the obtained images

        :parameters:
            n: (int): index of the target
            puponly: (bool): (optional) compute the PSF of the pupil (Airy) (default: False)
            comp_le: (bool): (optional) compute the long exposure image (default: True)
            comp_strehl: (bool): (optional) compute also the Strehl ratio (default: True)
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]

        src.comp_image(puponly, comp_le)
        if(comp_strehl):
            src.comp_strehl()

    def comp_strehl(self, int n):
        """
        Compute the Strehl Ratio of the PSF short exposure and long exposure

        :param n: (int): index of the target
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]

        src.comp_strehl()
    """
           _______. _______ .___________.                _______  _______ .___________.
          /       ||   ____||           |     ___       /  _____||   ____||           |
         |   (----`|  |__   `---|  |----`    ( _ )     |  |  __  |  |__   `---|  |----`
          \   \    |   __|      |  |         / _ \/\   |  | |_ | |   __|      |  |
      .----)   |   |  |____     |  |        | (_>  <   |  |__| | |  |____     |  |
      |_______/    |_______|    |__|         \___/\/    \______| |_______|    |__|

    """

    def get_image(self, int n, bytes type_im):
        """Return the image from the target (or long exposure image according to the requested type)

        :parameters:
            n: (int) : index of the target
            type_im: (str) : type of the image to get ("se" or "le")
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]

        cdef const long * dims = src.d_image.getDims()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F = np.empty((dims[2], dims[1]), dtype=np.float32)

        if(type_im == b"se"):
            src.d_image.device2host(< float*> data_F.data)
        if(type_im == b"le"):
            src.d_leimage.device2host(< float*> data_F.data)
            data_F /= src.strehl_counter

        return np.fft.fftshift(data_F.T.copy())

    def set_ncpa_phase(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            TODO : docstring

            :parameters:
                n : (int) : target number
                data: (np.ndarray[ndim=2, dtype=np.float32_t]) : ncpa phase
        """
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data_F = data.T.copy()
        self.target.d_targets[n].set_ncpa_phase(< float*>data_F.data, data.size)

    def get_ncpa_phase(self, int n):
        """
            TODO : docstring

            :parameters:
                n : (int) : target number
        """

        cdef carma_obj[float] * ph
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data

        ph = self.target.d_targets[n].d_phase.d_screen
        cdims = ph.getDims()
        data = np.empty((cdims[2], cdims[1]), dtype=np.float32)
        self.target.d_targets[n].get_ncpa_phase(< float*>data.data, data.size)
        return data.T.copy()

    def get_phase(self, int n):
        """Return the phase's screen of the target

        :param n: (int) : index of the target
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]

        cdef const long * dims
        dims = src.d_phase.d_screen.getDims()
        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        src.d_phase.d_screen.device2host(< float * > data_F.data)

        return data_F.T.copy()

    def set_phase(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """Set the phase's screen of the target

        :param n: (int) : index of the target
        :param data: (np.ndarray[ndim=2,dtype=np.float32_t]) : phase screen
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]

        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F = data.T.copy()

        src.d_phase.d_screen.host2device(< float * > data_F.data)

    def set_pupil(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """Set the pupil used for PSF computation

        :param n: (int) : index of the target
        :param data: (np.ndarray[ndim=2,dtype=np.float32_t]) : pupil
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef np.ndarray[dtype= np.int32_t] wherephase
        cdef long Npts
        cdef const long * dims
        cdef long dims1[2]
        dims = src.d_pupil.getDims()
        if((dims[2], dims[1]) == np.shape(data)):
            data_F = data.T.copy()
            src.d_pupil = new carma_obj[float](src.current_context, dims)
            src.d_pupil.host2device( < float * > data_F.data)
            wherephase = np.where(data_F)[0].astype(np.int32)
            Npts = wherephase.size
            dims1[0] = 1
            dims1[1] = Npts
            del src.d_phasepts
            del src.d_wherephase
            src.d_phasepts = new carma_obj[float](src.current_context, dims1)
            src.d_wherephase = new carma_obj[int](src.current_context, dims1)
            src.d_wherephase.host2device( < int * > wherephase.data)
        else:
            raise IndexError("Pupil dimension mismatch")

    def get_phasetele(self, int n):
        """Return the telemetry phase of the target

        :param n: (int) : index of the target
        :return data: (np.ndarray(ndim=2,np.float32)) : phase screen
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]

        cdef const long * dims
        dims = src.phase_telemetry.getDims()
        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)

        src.phase_telemetry.fill_into(< float * > data_F.data)

        return data_F.T.copy()

    def get_amplipup(self, int n):
        """Return the complex amplitude in the pupil plane of the target.

        :param n: (int) : index of the target
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]

        cdef const long * dims
        dims = src.d_amplipup.getDims()

        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.complex64)
        src.d_amplipup.device2host(< cuFloatComplex * > data_F.data)

        return data_F.T.copy()

    def get_strehl(self, int n):
        """Compute and return the target's strehl

        :param n: (int) : index of the target
        :return strehl: (np.array(4,dtype=np.float32)) : [Strehl SE, Strehl LE,
                                                        instantaneous phase variance over the pupil,
                                                        average  phase variance over the pupil]
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.target.d_targets[n]

        cdef np.ndarray strehl = np.empty(4, dtype=np.float32)
        strehl[0] = src.strehl_se
        strehl[1] = src.strehl_le
        strehl[2] = src.phase_var
        if(src.phase_var_count > 0):
            strehl[3] = src.phase_var_avg / float(src.phase_var_count)
        else:
            strehl[3] = 0.

        return strehl

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
        IF USE_BRAMA == 0:
            raise EnvironmentError(
                "Error: Trying to run a BRAMA Target over NO-BRAMA compile.")
        ELSE:
            Target.__init__(
                self,
                ctxt,
                telescope,
                ntargets,
                xpos,
                ypos,
                Lambda,
                mag,
                zerop,
                size,
                Npts,
                device)
            del self.target

            self.target = new sutra_target_brama(ctxt.c, "target_brama", telescope.telescope, 1, ntargets,
                                                 < float * > xpos.data, < float * > ypos.data,
                                                 < float * > Lambda.data, < float * > mag.data,
                                                 zerop, < long * > size.data,
                                                 Npts, ctxt.c.get_activeDevice())

    IF USE_BRAMA == 1:
        def __dealloc__(self):
            pass  # del self.target

        cpdef publish(self):
            cdef sutra_target_brama * target = < sutra_target_brama * > (self.target)
            target.publish()
