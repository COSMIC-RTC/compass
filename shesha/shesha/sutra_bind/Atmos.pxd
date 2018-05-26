include "sutra.pxd"
from naga.context cimport context as naga_context

cimport numpy as np

#################################################
# P-Class Atmos
#################################################
cdef class Atmos:
    cdef sutra_atmos * s_a
    """ sutra_atmos object"""
    cdef naga_context context
    """ naga_context """
