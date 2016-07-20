cimport numpy as np

from shesha_param import *
from shesha_param cimport *
from shesha_atmos import *
from shesha_atmos cimport *
from shesha_telescope import *
from shesha_telescope cimport *
from shesha_rtc import *
from shesha_rtc cimport *
from shesha_target import *
from shesha_target cimport *
from shesha_sensors import *
from shesha_sensors cimport *

cdef class Roket:

    cdef sutra_roket *roket
    cdef int device

    cpdef getContributor(self, bytes type)
    cpdef get_Btt(self)
    cpdef get_P(self)
    cpdef get_RD(self)
    cpdef get_gRD(self)
    cpdef computeBreakdown(self)
