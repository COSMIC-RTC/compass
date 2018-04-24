include "sutra.pxd"

from Telescope import *
from Telescope cimport *
from Atmos import *
from Atmos cimport *
from Dms import *
from Dms cimport *



#################################################
# P-Class Target
#################################################
cdef class Target:
    cdef sutra_target * target
    cdef int device
    cdef naga_context context


IF USE_BRAHMA == 1:
    cdef extern from * :
        sutra_target_brahma * dynamic_cast_target_brahma_ptr "dynamic_cast<sutra_target_brahma*>" (sutra_target *) except NULL

    #################################################
    # P-Class Target_brahma
    #################################################
    cdef class Target_brahma(Target):
        cpdef publish(self)
