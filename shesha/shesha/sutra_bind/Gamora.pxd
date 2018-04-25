cimport numpy as np
from naga.context cimport context as naga_context
from naga.obj cimport *
from naga.sparse_obj cimport *

include "sutra.pxd"

cdef class Gamora:

    cdef sutra_gamora *gamora
    cdef int device
