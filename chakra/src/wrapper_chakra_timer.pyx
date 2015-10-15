#################################################
# P-Class chakra_timer
#################################################

cdef class chakra_timer:

    def __cinit__(self):
        self.timer=new carma_timer()

    def start(self):
        self.timer.start()

    def reset(self):
        self.timer.reset()

    def stop(self):
        self.timer.stop()
        return self.timer.elapsed()


def threadSync():
    carmaSafeDeviceSynchronize()
