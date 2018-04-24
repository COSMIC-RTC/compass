cimport numpy as np
from naga_context cimport *
from naga_obj cimport *

include "sutra.pxd"

cdef class GrootS:

    cdef sutra_groot[float] *groot
    cdef int device

cdef class GrootD:

    cdef sutra_groot[double] *groot
    cdef int device
