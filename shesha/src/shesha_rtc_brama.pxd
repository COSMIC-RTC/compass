# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:27:07 2016

@author: sevin
"""

include "sutra_brama.pxd"

from shesha_param import *
from shesha_param cimport *
from shesha_rtc import *
from shesha_rtc cimport *

cdef extern from *:
    sutra_rtc_brama* dynamic_cast_rtc_brama_ptr "dynamic_cast<sutra_rtc_brama*>" (sutra_rtc*) except NULL

#################################################
# P-Class Rtc_brama
#################################################
cdef class Rtc_brama(Rtc):
    cpdef publish(self)
