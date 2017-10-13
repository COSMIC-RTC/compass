from naga_context cimport *

cimport numpy as np

include "sutra.pxd"

from libc.stdlib cimport malloc, free
from libcpp.cast cimport dynamic_cast
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from libc.stdint cimport uintptr_t

from cpython.string cimport PyString_AsString

from libc.math cimport sin

from naga_obj cimport naga_obj_Float2D 

from Telescope import *
from Telescope cimport *
#from shesha_param import *
#from shesha_param cimport *
from Atmos import *
from Atmos cimport *
from Dms cimport *


#################################################
# Dynamic casts
#################################################
cdef extern from * :
    sutra_wfs_geom * dynamic_cast_wfs_geom_ptr "dynamic_cast<sutra_wfs_geom*>" (sutra_wfs *) except NULL
    sutra_wfs_sh * dynamic_cast_wfs_sh_ptr "dynamic_cast<sutra_wfs_sh*>" (sutra_wfs *) except NULL
    sutra_wfs_pyr * dynamic_cast_wfs_pyr_ptr "dynamic_cast<sutra_wfs_pyr*>" (sutra_wfs *) except NULL
    sutra_wfs_pyr_pyr4 * dynamic_cast_wfs_pyr_pyr4_ptr "dynamic_cast<sutra_wfs_pyr_pyr4*>" (sutra_wfs *) except NULL
    sutra_wfs_pyr_pyrhr * dynamic_cast_wfs_pyr_pyrhr_ptr "dynamic_cast<sutra_wfs_pyr_pyrhr*>" (sutra_wfs *) except NULL
    sutra_wfs_pyr_roof * dynamic_cast_wfs_pyr_roof_ptr "dynamic_cast<sutra_wfs_pyr_roof*>" (sutra_wfs *) except NULL


ctypedef fused float_or_complex:
    np.ndarray[ndim=2, dtype=np.float32_t]
    np.ndarray[ndim=2, dtype=np.complex64_t]

#################################################
# P-Class Sensors
#################################################
cdef class Sensors:

    cdef sutra_sensors * sensors
    cdef naga_context context
    cdef _get_bincube(self, int n)
    cdef _get_pyrimg(self, int n)
    cdef _set_pyrimg(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data)
    cdef _set_submask(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data)
    cdef _get_submask(self, int n)
    cdef _get_pyrimghr(self, int n)
    cdef _get_slopesDims(self, int n)
    cdef _get_slopes(self, int n)
    cdef _get_hrmap(self, int n)
    cpdef copy_pyrimg(self, int n, naga_obj_Float2D data,  naga_obj_Int1D validx, naga_obj_Int1D validy)
    cdef _set_validpix(self, int n, np.ndarray[ndim=1, dtype=np.int32_t] datax, np.ndarray[ndim=1, dtype=np.int32_t] datay) 
