cimport numpy as np

from shesha_param import *
from shesha_param cimport *

from shesha_sensors cimport *

cdef class Acquisition:

    cdef sutra_acquisim * acquisition
    """sutra_acquisim object"""
