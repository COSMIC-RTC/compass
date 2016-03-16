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

from shesha_param import *
from shesha_param cimport *
from shesha_sensors cimport *
from shesha_dms cimport *
from shesha_target cimport *

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

    #cdef sensors_initbcube(self,int ncentro)
    cpdef getcentroids(self,int ncontrol, Sensors g_wfs=?, int nwfs=?)
    cpdef docentroids(self,int ncontrol=?)
    cpdef docentroids_geom(self,int ncontrol=?)
    cpdef init_proj(self,int i,Dms dms,np.ndarray[ndim=1,dtype=np.int32_t] indx_dm,
            np.ndarray[ndim=1,dtype=np.float32_t] unitpervolt,
            np.ndarray[ndim=1,dtype=np.int32_t] indx_pup,
            np.ndarray[ndim=1,dtype=np.int32_t] indx_mpup)
    cpdef init_modalOpti(self,int ncontro,int nmodes,int nrec, np.ndarray[ndim=2,dtype=np.float32_t] M2V,
            float gmin, float gmax, int ngain, float Fs)
    cpdef loadOpenLoop(self,int ncontro, np.ndarray[ndim=2, dtype=np.float32_t] ol_slopes)
    cpdef modalControlOptimization(self,int ncontro)
    cpdef set_gain(self,int ncontro, float gain)
    cpdef set_mgain(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] mgain)
    cpdef setCom(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] comvec)
    cpdef setCentroids(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] centro)
    cpdef get_mgain(self,int ncontro)
    cpdef set_imat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data)
    cpdef get_imat(self, int ncontro)
    cpdef set_cmat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data)
    cpdef get_cmat(self,int ncontro)
    cpdef get_cphim(self, int ncontro)
    cpdef get_cmm(self,int ncontro)
    cpdef set_cmm(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data)
    cpdef set_decayFactor(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] decay)
    cpdef set_matE(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] matE)
    cpdef set_openloop(self,int ncontro, int openloop)


    cpdef doimat_geom(self, int ncontro, Dms g_dms,int geom)
    cpdef doimat(self, int ncontro, Dms g_dms)
    cpdef sensors_compslopes(self, int ncentro, int nmax=?, float thresh=?)
    cdef add_controller(self, int nactu, float delay, bytes type_control, Dms dms,
                 char **type_dmseen, np.ndarray[ndim=1,dtype=np.float32_t] alt,
                 int ndm, long Nphi=?, bool wfs_direction=?)


    cpdef imat_svd(self,int ncontro)
    cpdef setU(self,int ncontro,np.ndarray[ndim=2,dtype=np.float32_t] U)
    cpdef setEigenvals(self, int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] eigenvals)
    cpdef getU(self, int ncontro)
    cpdef getEigenvals(self,int ncontro)
    cpdef getCmmEigenvals(self,int ncontro)
    cpdef getCenbuff(self, int ncontro)
    cpdef getErr(self,int ncontro)
    cpdef getCom(self,int ncontro)
    cpdef getolmeas(self,int ncontro)
    cpdef getVoltage(self,int ncontro)
    cpdef getCentroids(self,int ncontro)
    cpdef buildcmat(self,int ncontro,int nfilt, int filt_tt=?)
    cpdef buildcmatmv(self,int ncontro,float cond)
    cpdef loadnoisemat(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] N)
    cpdef docontrol(self,int ncontro)
    cpdef docontrol_geo(self, int ncontro, Dms dms, Target target, int ntarget)
    cpdef docontrol_geo_onwfs(self, int ncontro, Dms dms, Sensors wfs, int nwfs)
    cpdef applycontrol(self,int ncontro,Dms dms)
    cpdef get_nfiltered(self,int ncontro,Param_rtc p_rtc)

IF USE_BRAMA == 1:
        cdef extern from *:
            sutra_rtc_brama* dynamic_cast_rtc_brama_ptr "dynamic_cast<sutra_rtc_brama*>" (sutra_rtc*) except NULL

        #################################################
        # P-Class Rtc_brama
        #################################################
        cdef class Rtc_brama(Rtc):
            cpdef publish(self)

