

cdef class Telescope:

    def __cinit__(self, long pup_size, long num_eleme_pup, np.ndarray[ndim=2,dtype=np.float32_t] pupil,
                  np.ndarray[ndim=2,dtype=np.float32_t] phase_ab_M1,
                  long pup_size_m,
                  np.ndarray[ndim=2,dtype=np.float32_t] pupil_m,
                  np.ndarray[ndim=2,dtype=np.float32_t] phase_ab_m1_m):
        cdef carma_context *context=carma_context.instance()

        cdef np.ndarray[ndim=1,dtype=np.float32_t] pupil_F= pupil.flatten("F")
        cdef np.ndarray[ndim=1,dtype=np.float32_t] phase_ab_M1_F= phase_ab_M1.flatten("F")
        cdef np.ndarray[ndim=1,dtype=np.float32_t] pupil_m_F= pupil_m.flatten("F")
        cdef np.ndarray[ndim=1,dtype=np.float32_t] phase_ab_m1_m_F= phase_ab_m1_m.flatten("F")

        self.telescope= new sutra_telescope(context, pup_size,  num_eleme_pup,
                         <float*>pupil_F.data, <float*>phase_ab_M1_F.data, pup_size_m,
                         <float*>pupil_m_F.data, <float*>phase_ab_m1_m_F.data)


    cpdef set_pupil(self,np.ndarray[ndim=2,dtype=np.float32_t] pup):
        """ Set the small size pupil

        :parameters:
            pup: (np.ndarray[ndim=2,dtype=np.float32_t]) :  Telescope pupil (small size) to be sent to for target image computation
        """

        cdef np.ndarray[ndim=1,dtype=np.float32_t] pup_F=pup.flatten("F")
        self.telescope.d_pupil.host2device(<float*>pup_F.data)

#    cpdef set_pupil_m(self,np.ndarray[ndim=2,dtype=np.float32_t] pup_m)
#    cpdef set_phase_ab_M1(self,np.ndarray[ndim=2,dtype=np.float32_t] phase)
#    cpdef set_phase_ab_M1_m(self,np.ndarray[ndim=2,dtype=np.float32_t] phase_m)
