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

from shesha_telescope import *
from shesha_telescope cimport *
from shesha_param import *
from shesha_param cimport *
from shesha_atmos import *
from shesha_atmos cimport *
from shesha_dms cimport *

cdef np.float32_t dtor = (np.pi / 180.)


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


#################################################
# P-Class Sensors
#################################################
cdef class Sensors:

    cdef sutra_sensors * sensors
    cpdef sensors_initgs(self, np.ndarray[ndim=1, dtype=np.float32_t] xpos,
                             np.ndarray[ndim=1, dtype=np.float32_t] ypos,
                             np.ndarray[ndim=1, dtype=np.float32_t] Lambda,
                             np.ndarray[ndim=1, dtype=np.float32_t] mag,
                             float zerop,
                             np.ndarray[ndim=1, dtype=np.int64_t  ] size,
                             np.ndarray[ndim=1, dtype=np.float32_t] noise,
                             np.ndarray[ndim=1, dtype=np.int64_t  ] seed,
                             np.ndarray[ndim=1, dtype=np.float32_t] G,
                             np.ndarray[ndim=1, dtype=np.float32_t] thetaML,
                             np.ndarray[ndim=1, dtype=np.float32_t] dx,
                             np.ndarray[ndim=1, dtype=np.float32_t] dy)



    cpdef sensors_initarr(self, int n, Param_wfs wfs)
    cpdef sensors_addlayer(self, int i, bytes type, float alt, float xoff, float yoff)
    cpdef comp_modulation(self, int n, int cpt)
    cdef _get_bincube(self, int n)
    cdef _get_pyrimg(self, int n)
    cdef _set_pyrimg(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data)
    cdef _set_submask(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data)
    cdef _get_submask(self, int n)
    cdef _get_pyrimghr(self, int n)
#    cdef _set_pyrimghr(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data)
    cpdef get_binimg(self, int n, Telescope tel=?, Atmos atmos=?, Dms dms=?)
    cpdef _get_slopesDims(self, int n)
    cdef _get_slopes(self, int n)
    cpdef slopes_geom(self, int nsensors, int t)
    cpdef sensors_trace(self, int n, str type_trace, Telescope tel=?, Atmos atmos=?, Dms dms=?, int rst=?)
    cpdef get_bincubeNotNoisy(self, int n)
    cpdef get_binimg_notnoisy(self, int n)
    cpdef set_bincube(self, int n, np.ndarray[ndim=3, dtype=np.float32_t] data)
    IF USE_MPI:
        cpdef gather_bincube(self, int n)
        cpdef gather_bincube_cuda_aware(self, int n)
        cpdef Bcast_dscreen(self)
        cpdef Bcast_dscreen_cuda_aware(self)
    cdef _get_rank(self, int n)

    # for profiling purpose
    '''
    cdef gather_bincube_prof(self,int n)
    cdef wait1_prof(self)
    cdef wait2_prof(self)
    cdef d2h_prof(self,float *ptr,n)
    cdef h2d_prof(self,float *ptr,n)
    cdef gather_prof(self,float *ptr, int size, int *count, int *disp)
    '''

    cdef  _get_hrmap(self, int n)
    # cdef getDims(self)

cpdef noise_cov(int nw, Param_wfs p_wfs, Param_atmos p_atmos, Param_tel p_tel)

cpdef prep_lgs_prof(Param_wfs p_wfs, int nsensors, Param_tel p_tel,
                    np.ndarray[dtype=np.float32_t] prof,
                    np.ndarray[dtype=np.float32_t] h,
                    float beam, Sensors sensors,
                    bytes center=?, int imat=?)
