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


IF USE_BRAMA == 1:
    cdef extern from * :
        sutra_target_brama * dynamic_cast_target_brama_ptr "dynamic_cast<sutra_target_brama*>" (sutra_target *) except NULL

    #################################################
    # P-Class Target_brama
    #################################################
    cdef class Target_brama(Target):
        cpdef publish(self)
