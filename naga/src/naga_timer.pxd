cdef extern from "carma_utils.h":
    void carmaSafeDeviceSynchronize()

#################################################
# C-Class carma_timer
#################################################
cdef extern from "carma_timer.h":
    cdef cppclass carma_timer:
        pass
        void start()
        void reset()
        void stop()
        double elapsed()

#################################################
# P-Class naga_timer
#################################################
cdef class naga_timer:
    cdef carma_timer *timer



