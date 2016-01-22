# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 11:19:55 2016

@author: sevin
"""

from naga cimport *

import numpy as np
cimport numpy as np

include "../par.pxi"

from libc.stdlib cimport malloc, free
from libcpp.cast cimport dynamic_cast
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from libc.stdint cimport uintptr_t

from cpython.string cimport PyString_AsString

from libc.math cimport sin

cdef np.float32_t dtor = (np.pi/180.)

from shesha_telescope import *
from shesha_telescope cimport *
from shesha_param import *
from shesha_param cimport *
from shesha_sensors cimport *
from shesha_atmos cimport *
from shesha_dms cimport *

#################################################
# P-Class Target_brama
#################################################
cdef class Target_brama(Target):
    cpdef publish(self)
