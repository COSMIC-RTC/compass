include "sutra.pxd"
from naga_context cimport *
from naga_obj cimport *
from naga_sparse_obj cimport *
from naga_magma cimport *

cimport numpy as np

# cpdef long RASC = 180.*3600./np.pi

from libc.stdlib cimport malloc, free
from libcpp.cast cimport dynamic_cast
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from libc.stdint cimport uintptr_t

from cpython.string cimport PyString_AsString

from Sensors cimport *
from Dms cimport *
from Target cimport *


#################################################
# Dynamic casts
#################################################
cdef extern from * :
    sutra_centroider_tcog * dynamic_cast_centroider_tcog_ptr "dynamic_cast<sutra_centroider_tcog*>" (sutra_centroider *) except NULL
    sutra_centroider_bpcog * dynamic_cast_centroider_bpcog_ptr "dynamic_cast<sutra_centroider_bpcog*>" (sutra_centroider *) except NULL
    sutra_centroider_corr * dynamic_cast_centroider_corr_ptr "dynamic_cast<sutra_centroider_corr*>" (sutra_centroider *) except NULL
    sutra_centroider_wcog * dynamic_cast_centroider_wcog_ptr "dynamic_cast<sutra_centroider_wcog*>" (sutra_centroider *) except NULL
    sutra_centroider_pyr * dynamic_cast_centroider_pyr_ptr "dynamic_cast<sutra_centroider_pyr*>" (sutra_centroider *) except NULL
    sutra_controller_generic * dynamic_cast_controller_generic_ptr "dynamic_cast<sutra_controller_generic*>" (sutra_controller *) except NULL
    sutra_controller_geo * dynamic_cast_controller_geo_ptr "dynamic_cast<sutra_controller_geo*>" (sutra_controller *) except NULL
    sutra_controller_ls * dynamic_cast_controller_ls_ptr "dynamic_cast<sutra_controller_ls*>" (sutra_controller *) except NULL
    sutra_controller_mv * dynamic_cast_controller_mv_ptr "dynamic_cast<sutra_controller_mv*>" (sutra_controller *) except NULL
    sutra_controller_cured * dynamic_cast_controller_cured_ptr "dynamic_cast<sutra_controller_cured*>" (sutra_controller *) except NULL
    sutra_controller_kalman * dynamic_cast_controller_kl_ptr "dynamic_cast<sutra_controller_kalman*>" (sutra_controller *) except NULL





#################################################
# P-Class Rtc
#################################################
cdef class Rtc:
    cdef naga_context context
    cdef sutra_rtc * rtc
    cdef int device

IF USE_BRAHMA == 1:
        cdef extern from * :
            sutra_rtc_brahma * dynamic_cast_rtc_brahma_ptr "dynamic_cast<sutra_rtc_brahma*>" (sutra_rtc *) except NULL

        #################################################
        # P-Class Rtc_brahma
        #################################################
        cdef class Rtc_brahma(Rtc):
            cpdef publish(self)
