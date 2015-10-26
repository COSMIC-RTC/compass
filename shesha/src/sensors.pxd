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

IF USE_MPI == 1:
    print "using mpi"
    #from mpi4py import MPI
    from mpi4py cimport MPI
    #cimport mpi4py.MPI as MPI
    # C-level cdef, typed, Python objects
    from mpi4py cimport mpi_c as mpi
    #from mpi4py cimport libmpi as mpi
    #cimport mpi4py.libmpi as mpi


from libc.math cimport sin

from param import *
from param cimport *
from atmos import *
from atmos cimport *
from dms cimport *

cdef np.float32_t dtor = (np.pi/180.)


#################################################
# Dynamic casts
#################################################
cdef extern from *:
    sutra_wfs_geom* dynamic_cast_wfs_geom_ptr "dynamic_cast<sutra_wfs_geom*>" (sutra_wfs*) except NULL
    sutra_wfs_sh* dynamic_cast_wfs_sh_ptr "dynamic_cast<sutra_wfs_sh*>" (sutra_wfs*) except NULL
    sutra_wfs_pyr_pyr4* dynamic_cast_wfs_pyr_pyr4_ptr "dynamic_cast<sutra_wfs_pyr_pyr4*>" (sutra_wfs*) except NULL
    sutra_wfs_pyr_roof* dynamic_cast_wfs_pyr_roof_ptr "dynamic_cast<sutra_wfs_pyr_roof*>" (sutra_wfs*) except NULL








#################################################
# P-Class Sensors
#################################################
cdef class Sensors:
    
    cdef sutra_sensors *sensors
    cpdef sensors_initgs(self,np.ndarray[ndim=1,dtype=np.float32_t] xpos,
                             np.ndarray[ndim=1,dtype=np.float32_t] ypos,
                             np.ndarray[ndim=1,dtype=np.float32_t] Lambda,
                             np.ndarray[ndim=1,dtype=np.float32_t] mag,
                             float zerop,
                             np.ndarray[ndim=1,dtype=np.int64_t  ] size,
                             np.ndarray[ndim=1,dtype=np.float32_t] noise,
                             np.ndarray[ndim=1,dtype=np.int64_t  ] seed)


    cpdef sensors_initarr(self,int n, Param_wfs wfs, Param_geom geom)
    cpdef sensors_addlayer(self,int i, bytes type, float alt, float xoff, float yoff)
    cdef _get_bincube(self, int n)
    cdef _get_pyrimg(self,int n)
    cdef _get_binimg(self, int n)
    cdef _get_slopesDims(self,int n)
    cdef _get_slopes(self, int n)
    cpdef slopes_geom(self,int nsensors, int t)
    cpdef sensors_trace(self,int n, str type_trace, Atmos atmos=?, Dms dms=?, int rst=?)
    IF USE_MPI==1:
        cpdef gather_bincube(self,int n)
        cpdef gather_bincube_cuda_aware(self,int n)
        cpdef Bcast_dscreen(self)
        cpdef Bcast_dscreen_cuda_aware(self)
    cdef _get_rank(self,int n)

    #for profiling purpose
    '''
    cdef gather_bincube_prof(self,int n)
    cdef wait1_prof(self)
    cdef wait2_prof(self)
    cdef d2h_prof(self,float *ptr,n)
    cdef h2d_prof(self,float *ptr,n)
    cdef gather_prof(self,float *ptr, int size, int *count, int *disp)
    '''

    cdef  _get_hrmap(self, int n)
    #cdef getDims(self)

cpdef noise_cov(int nw, Param_wfs p_wfs, Param_atmos p_atmos, Param_tel p_tel)

cpdef prep_lgs_prof(Param_wfs p_wfs,int nsensors, Param_tel p_tel, 
                    np.ndarray[dtype=np.float32_t] prof,
                    np.ndarray[dtype=np.float32_t] h,
                    float beam, Sensors sensors,
                    bytes center=?, int imat=?)
