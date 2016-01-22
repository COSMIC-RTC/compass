# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 11:32:51 2016

@author: sevin
"""

include "../par.pxi"
include "sutra.pxd"

from naga cimport *

import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc, free
from libcpp.cast cimport dynamic_cast
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from libc.stdint cimport uintptr_t

from cpython.string cimport PyString_AsString

#################################################
# C-Class sutra_rtc_brama
#################################################
cdef extern from "sutra_rtc_brama.h":
    cdef cppclass sutra_rtc_brama(sutra_rtc):
        
        sutra_rtc_brama(carma_context *context, sutra_sensors *wfs, sutra_target *target, char* name)

        void publish();

#################################################
# C-Class sutra_target_brama
#################################################
cdef extern from "sutra_target_brama.h":
    cdef cppclass sutra_target_brama(sutra_target):
        
        sutra_target_brama(carma_context *context, char* name, sutra_telescope *d_tel, int subsample_, int ntargets, float *xpos, 
                           float *ypos, float *zlambda, float *mag, float zerop, long *sizes,
                           int Npts, int device)

        void publish();
