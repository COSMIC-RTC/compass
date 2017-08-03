import numpy as np
cimport numpy as np
# np.import_array()

cdef class Telescope:

    def __cinit__(self, naga_context context, long pup_size, long num_eleme_pup, np.ndarray[ndim=2, dtype=np.float32_t] pupil,
                  np.ndarray[ndim=2, dtype=np.float32_t] phase_ab_M1,
                  long pup_size_m,
                  np.ndarray[ndim=2, dtype=np.float32_t] pupil_m,
                  np.ndarray[ndim=2, dtype=np.float32_t] phase_ab_m1_m):
        self.context = context
        cdef np.ndarray[dtype = np.float32_t] pupil_F = pupil.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] phase_ab_M1_F = phase_ab_M1.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] pupil_m_F = pupil_m.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] phase_ab_m1_m_F = phase_ab_m1_m.flatten("F")

        self.telescope = new sutra_telescope(self.context.c, pup_size, num_eleme_pup,
                         < float *> pupil_F.data, < float *> phase_ab_M1_F.data, pup_size_m,
                         < float *> pupil_m_F.data, < float *> phase_ab_m1_m_F.data)


    def __dealloc__(self):
        del self.telescope

    def set_pupil(self, np.ndarray[ndim=2, dtype=np.float32_t] pup):
        """ Set the small pupil

        :parameters:
            pup: (np.ndarray[ndim=2,dtype=np.float32_t]) :  Telescope pupil (small size) to be sent to for target image computation
        """

        cdef np.ndarray[dtype = np.float32_t] pup_F = pup.flatten("F")
        self.telescope.d_pupil.host2device(< float *> pup_F.data)

    def set_pupil_m(self, np.ndarray[ndim=2, dtype=np.float32_t] pup):
        """ Set the medium pupil

        :parameters:
            pup: (np.ndarray[ndim=2,dtype=np.float32_t]) :  Telescope pupil (medium size) to be sent to for wfs image computation
        """

        cdef np.ndarray[dtype = np.float32_t] pup_F = pup.flatten("F")
        self.telescope.d_pupil_m.host2device(< float *> pup_F.data)

    def set_phase_ab_M1(self, np.ndarray[ndim=2, dtype=np.float32_t] phase_ab):
        """ Set the M1 phase aberration in the small pupil

        :parameters:
            phase_ab: (np.ndarray[ndim=2,dtype=np.float32_t]) : M1 phase aberration in the small pupil
        """

        cdef np.ndarray[dtype = np.float32_t] phase_ab_F = phase_ab.flatten("F")
        self.telescope.d_phase_ab_M1.host2device(< float *> phase_ab_F.data)

    def set_phase_ab_M1_m(self, np.ndarray[ndim=2, dtype=np.float32_t] phase_ab):
        """ Set the M1 phase aberration in the medium pupil

        :parameters:
            phase_ab: (np.ndarray[ndim=2,dtype=np.float32_t]) : M1 phase aberration in the medium pupil
        """

        cdef np.ndarray[dtype = np.float32_t] phase_ab_F = phase_ab.flatten("F")
        self.telescope.d_phase_ab_M1_m.host2device(< float *> phase_ab_F.data)

    def get_pupil(self):
        """Return the small pupil of the sutra_telescope object

        :return:
           pup : (np.ndarray[ndim=2,dtype=np.float32_t]) :pupil

        """
        self.context.set_activeDeviceForCpy(self.telescope.device, 1)

        cdef const long * dims = NULL
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] pup_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] pup

        dims = self.telescope.d_pupil.getDims()
        pup_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.telescope.d_pupil.device2host(< float *> pup_F.data)


        pup = np.reshape(pup_F.flatten("F"), (dims[1], dims[2]))
        return pup

    def get_mpupil(self):
        """Return the medium pupil of the sutra_telescope object

        :return:
           mpup : (np.ndarray[ndim=2,dtype=np.float32_t]) : pupil (medium size)

        """
        self.context.set_activeDeviceForCpy(self.telescope.device, 1)

        cdef const long * dims = NULL
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] pup_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] pup

        dims = self.telescope.d_pupil_m.getDims()
        pup_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.telescope.d_pupil_m.device2host(< float *> pup_F.data)


        pup = np.reshape(pup_F.flatten("F"), (dims[1], dims[2]))
        return pup

    def get_phase_ab(self):
        """Return the M1 phase aberration of the sutra_telescope object

        :return:
           phase_ab : (np.ndarray[ndim=2,dtype=np.float32_t]) : M1 phase aberration

        """
        self.context.set_activeDeviceForCpy(self.telescope.device, 1)

        cdef const long * dims = NULL
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] phase_ab_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] phase_ab

        dims = self.telescope.d_phase_ab_M1.getDims()
        phase_ab_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.telescope.d_phase_ab_M1.device2host(< float *> phase_ab_F.data)


        phase_ab = np.reshape(phase_ab_F.flatten("F"), (dims[1], dims[2]))
        return phase_ab

    def get_phase_ab_m(self):
        """Return the M1 phase aberration of the sutra_telescope object (medium size)

        :return:
           phase_ab : (np.ndarray[ndim=2,dtype=np.float32_t]) : M1 phase aberration (medium size)

        """
        self.context.set_activeDeviceForCpy(self.telescope.device, 1)

        cdef const long * dims = NULL
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] phase_ab_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] phase_ab

        dims = self.telescope.d_phase_ab_M1_m.getDims()
        phase_ab_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
        self.telescope.d_phase_ab_M1_m.device2host(< float *> phase_ab_F.data)


        phase_ab = np.reshape(phase_ab_F.flatten("F"), (dims[1], dims[2]))
        return phase_ab
