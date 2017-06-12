cimport numpy as np
from naga_context cimport *
from naga_obj cimport *

include "sutra.pxd"

cdef class Groot:

    cdef sutra_groot *groot
    cdef int device

    cpdef compute_Cerr(self)
