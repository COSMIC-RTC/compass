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

from param cimport *
from param import *
from sensors cimport *


#################################################
# P-Class Dms
#################################################
cdef class Dms:
    cdef sutra_dms *dms
    """sutra_dms object"""
    cdef int device
    """ GPU device number"""
    cpdef add_dm(self, bytes type_dm, float alt, long dim, long ninflu, 
                long influsize, long ninflupos, long npts, float puhs4imat,
                int device=?)
    cpdef remove_dm(self,bytes type_dm, float alt)

    cpdef load_pzt(self, float alt,
                    np.ndarray[ndim=3,dtype=np.float32_t] influ,
                    np.ndarray[ndim=1,dtype=np.int32_t] influpos,
                    np.ndarray[ndim=1,dtype=np.int32_t] npoints,
                    np.ndarray[ndim=1,dtype=np.int32_t] istart,
                    np.ndarray[ndim=1,dtype=np.int32_t] xoff,
                    np.ndarray[ndim=1,dtype=np.int32_t] yoff,
                    np.ndarray[ndim=2,dtype=np.float32_t] kern)

    #TODO dims of arrays
    cpdef load_kl(self,float alt, np.ndarray[ndim=1,dtype=np.float32_t] rabas,
                    np.ndarray[ndim=1,dtype=np.float32_t] azbas,
                    np.ndarray[ndim=1,dtype=np.int32_t] ord,
                    np.ndarray[ndim=1,dtype=np.float32_t] cr,
                    np.ndarray[ndim=1,dtype=np.float32_t] cp)


    cpdef load_tt(self,float alt, np.ndarray[ndim=3,dtype=np.float32_t] influ)

    cpdef set_comm(self,bytes type_dm,float alt,
                    np.ndarray[ndim=1,dtype=np.float32_t] comm)
    cpdef shape_dm(self,bytes type_dm,float alt)

    cpdef computeKLbasis(self, bytes type_dm, float alt, 
        np.ndarray[ndim=1,dtype=np.float32_t] xpos, np.ndarray[ndim=1,dtype=np.float32_t] ypos,
        np.ndarray[ndim=1,dtype=np.int32_t] indx_pup, long dim, float norm, float ampli)
    cpdef get_KLbasis(self,bytes type_dm, float alt)

    cpdef getComm(self,bytes type_dm,float alt)
    cpdef getInflu(self,bytes type_dm,float alt)
    cpdef comp_oneactu(self,bytes type_dm, float alt, int nactu, float ampli)


cpdef comp_dmgeom(Param_dm dm, Param_geom geom)
cpdef compute_klbasis(Dms g_dm,Param_dm p_dm, Param_geom p_geom,Param_atmos p_atmos,Param_tel p_tel)
