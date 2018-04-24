include "sutra.pxd"

from naga_context cimport *
cimport numpy as np

from libc.stdlib cimport malloc, free
from libcpp.cast cimport dynamic_cast
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from libc.stdint cimport uintptr_t

from cpython.string cimport PyString_AsString

from libc.math cimport sin

from naga_sparse_obj cimport *

from Sensors import *
from Sensors cimport *


#################################################
# P-Class Dms
#################################################
cdef class Dms:
    cdef naga_context context
    cdef sutra_dms *dms
    """sutra_dms object"""
    cdef int device
    """ GPU device number"""
