cimport numpy as np

from Rtc cimport *
from Telescope cimport *
from Atmos cimport *
from Sensors cimport *
from Target cimport *
from Dms cimport *

cdef class Roket:

    cdef sutra_roket *roket
    cdef int device
