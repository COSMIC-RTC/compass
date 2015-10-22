from naga cimport *

import numpy as np
cimport numpy as np

include "../par.pxi"
#cpdef long RASC = 180.*3600./np.pi

from libc.stdlib cimport malloc, free
from libcpp.cast cimport dynamic_cast
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from libc.stdint cimport uintptr_t

from cpython.string cimport PyString_AsString

from param import *
from param cimport *
from sensors cimport *
from dms cimport *
from target cimport *


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

#################################################
# Dynamic casts
#################################################
cdef extern from *:
    sutra_centroider_tcog* dynamic_cast_centroider_tcog_ptr "dynamic_cast<sutra_centroider_tcog*>" (sutra_centroider*) except NULL
    sutra_centroider_corr* dynamic_cast_centroider_corr_ptr "dynamic_cast<sutra_centroider_corr*>" (sutra_centroider*) except NULL
    sutra_centroider_wcog* dynamic_cast_centroider_wcog_ptr "dynamic_cast<sutra_centroider_wcog*>" (sutra_centroider*) except NULL
    sutra_controller_generic* dynamic_cast_controller_generic_ptr "dynamic_cast<sutra_controller_generic*>" (sutra_controller*) except NULL
    sutra_controller_geo* dynamic_cast_controller_geo_ptr "dynamic_cast<sutra_controller_geo*>" (sutra_controller*) except NULL
    sutra_controller_ls* dynamic_cast_controller_ls_ptr "dynamic_cast<sutra_controller_ls*>" (sutra_controller*) except NULL
    sutra_controller_mv* dynamic_cast_controller_mv_ptr "dynamic_cast<sutra_controller_mv*>" (sutra_controller*) except NULL
    sutra_controller_cured* dynamic_cast_controller_cured_ptr "dynamic_cast<sutra_controller_cured*>" (sutra_controller*) except NULL
    sutra_controller_kalman* dynamic_cast_controller_kl_ptr "dynamic_cast<sutra_controller_kalman*>" (sutra_controller*) except NULL





#################################################
# P-Class Rtc
#################################################
cdef class Rtc:
    cdef sutra_rtc *rtc
    cdef int device
    cdef int use_brama
    #cdef sensors_initbcube(self,int ncentro)
    cpdef getcentroids(self,int ncontrol, Sensors g_wfs=?, int nwfs=?)
    cpdef docentroids(self,int ncontrol=?)
    cpdef docentroids_geom(self,int ncontrol=?)
    cdef init_proj(self,int i,Dms dms,np.ndarray[ndim=1,dtype=np.int32_t] indx_dm,
            np.ndarray[ndim=1,dtype=np.float32_t] unitpervolt, 
            np.ndarray[ndim=1,dtype=np.int32_t] indx_pup)
    cdef init_modalOpti(self,int ncontro,int nmodes,int nrec, np.ndarray[ndim=2,dtype=np.float32_t] M2V,
            float gmin, float gmax, int ngain, float Fs)
    cdef loadOpenLoop(self,int ncontro, np.ndarray[ndim=2, dtype=np.float32_t] ol_slopes)
    cdef modalControlOptimization(self,int ncontro)
    cdef set_gain(self,int ncontro, float gain)
    cdef set_mgain(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] mgain)
    cpdef get_mgain(self,int ncontro)
    cpdef set_imat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data)
    cpdef get_imat(self, int ncontro)
    cpdef set_cmat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data)
    cpdef get_cmat(self,int ncontro)
    cpdef get_cphim(self, int ncontro)
    cpdef get_cmm(self,int ncontro)
    cpdef set_cmm(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data)
    cdef set_decayFactor(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] decay)
    cdef set_matE(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] matE)


    cdef doimat_geom(self, int ncontro, Dms g_dms,int geom)
    cdef doimat(self, int ncontro, Dms g_dms)
    cdef sensors_compslopes(self, int ncentro, int nmax=?, float thresh=?)
    cdef add_controller(self, int nactu, float delay, bytes type_control, Dms dms,
                 char **type_dmseen, np.ndarray[ndim=1,dtype=np.float32_t] alt,
                 int ndm, long Nphi=?)


    cdef imat_svd(self,int ncontro)
    cdef setU(self,int ncontro,np.ndarray[ndim=2,dtype=np.float32_t] U)
    cdef setEigenvals(self, int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] eigenvals)
    cdef getU(self, int ncontro)
    cpdef getEigenvals(self,int ncontro)
    cpdef getCmmEigenvals(self,int ncontro)
    cpdef getCenbuff(self, int ncontro)
    cdef getErr(self,int ncontro)
    cpdef getCom(self,int ncontro)
    cpdef getolmeas(self,int ncontro)
    cpdef getVoltage(self,int ncontro)
    cpdef getCentroids(self,int ncontro)
    cdef buildcmat(self,int ncontro,int nfilt, int filt_tt=?)
    cpdef buildcmatmv(self,int ncontro,float cond)
    cdef loadnoisemat(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] N)
    cpdef docontrol(self,int ncontro)
    cpdef docontrol_geo(self, int ncontro, Dms dms, Target target, int ntarget)
    cpdef applycontrol(self,int ncontro,Dms dms)


