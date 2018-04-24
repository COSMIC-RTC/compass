include "sutra.pxd"
from naga_context cimport *

cimport numpy as np

from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair

from cpython.string cimport PyString_AsString

#################################################
# P-Class Atmos
#################################################
cdef class Atmos:
    cdef sutra_atmos * s_a
    """ sutra_atmos object"""
    cdef naga_context context
    """ context """
