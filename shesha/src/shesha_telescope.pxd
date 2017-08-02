from naga_context cimport *

cimport numpy as np

from shesha_param import *
from shesha_param cimport *

cdef class Telescope:

    cdef sutra_telescope * telescope
    """sutra_telescope object"""
    cdef naga_context context
