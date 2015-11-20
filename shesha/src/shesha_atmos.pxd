from naga cimport *

import numpy as np
cimport numpy as np

include "../par.pxi"

from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair

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


from shesha_param import *
from shesha_param cimport *





#################################################
# P-Class Atmos
#################################################
cdef class Atmos:
    cdef sutra_atmos *s_a
    """ sutra_atmos object"""
    cdef naga_context context
    """ context """
    cdef realinit(self,naga_context ctxt,int nscreens,
                np.ndarray[ndim=1,dtype=np.float32_t] r0,
                np.ndarray[ndim=1,dtype=np.int64_t] size,
                np.ndarray[ndim=1,dtype=np.float32_t] altitude,
                np.ndarray[ndim=1,dtype=np.float32_t] windspeed,
                np.ndarray[ndim=1,dtype=np.float32_t] winddir,
                np.ndarray[ndim=1,dtype=np.float32_t] deltax,
                np.ndarray[ndim=1,dtype=np.float32_t] deltay,
                int device )


