cimport numpy as np
from naga.context cimport context as naga_context
from naga.obj cimport *

include "sutra.pxd"

cdef class Groot:

    cdef sutra_groot *groot
    cdef int device

cdef class GrootAlias:

    cdef sutra_groot *groot
    cdef int device
