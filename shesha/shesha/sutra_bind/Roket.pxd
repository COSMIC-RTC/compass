cimport numpy as np

from shesha.sutra_bind.Rtc cimport *
from shesha.sutra_bind.Telescope cimport *
from shesha.sutra_bind.Atmos cimport *
from shesha.sutra_bind.Sensors cimport *
from shesha.sutra_bind.Target cimport *
from shesha.sutra_bind.Dms cimport *

cdef class Roket:

    cdef sutra_roket *roket
    cdef int device
