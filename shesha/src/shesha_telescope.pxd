from naga cimport *

import numpy as np
cimport numpy as np

from shesha_param import *
from shesha_param cimport *

cdef class Telescope:

    cdef sutra_telescope *telescope
    
    cpdef set_pupil(self,np.ndarray[ndim=2,dtype=np.float32_t] pup)
#    cpdef set_pupil_m(np.ndarray[ndim=2,dtype=np.float32_t] pup_m)
#    cpdef set_phase_ab_M1(np.ndarray[ndim=2,dtype=np.float32_t] phase)
#    cpdef set_phase_ab_M1_m(np.ndarray[ndim=2,dtype=np.float32_t] phase_m)
