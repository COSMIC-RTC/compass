cimport numpy as np
from naga.context cimport context as naga_context
from naga.obj cimport *

include "sutra.pxd"

cdef class GrootS:

    cdef sutra_groot[float] *groot
    cdef int device

cdef class GrootD:

    cdef sutra_groot[double] *groot
    cdef int device
