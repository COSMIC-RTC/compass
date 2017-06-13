cimport numpy as np
from naga_context cimport *
from naga_obj cimport *
from naga_sparse_obj cimport *

include "sutra.pxd"

cdef class Gamora:

    cdef sutra_gamora *gamora
    cdef int device

    cpdef psf_rec_roket(self,np.ndarray[ndim=2, dtype=np.float32_t] err)
    cpdef psf_rec_Vii(self)
