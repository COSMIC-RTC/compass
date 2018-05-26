include "sutra.pxd"

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

from naga.sparse_obj cimport sparse_obj_Double, sparse_obj_Float
from naga.context cimport context as naga_context

from shesha.sutra_bind.Sensors cimport *


#################################################
# P-Class Dms
#################################################
cdef class Dms:
    cdef naga_context context
    cdef sutra_dms *dms
    """sutra_dms object"""
    cdef int device
    """ GPU device number"""
