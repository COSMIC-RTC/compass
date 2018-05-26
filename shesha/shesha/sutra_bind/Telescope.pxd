include "sutra.pxd"
from naga.context cimport context as naga_context

cimport numpy as np

#################################################
# P-Class Telescope
#################################################
cdef class Telescope:

    cdef sutra_telescope * telescope
    """sutra_telescope object"""
    cdef naga_context context
    """ naga_context """
