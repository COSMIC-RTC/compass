cdef class Param_loop:
    def set_niter(self,long n):
        """Set the number of iteration
        
        :param n: (long) : number of iteration
        """
        self.niter=n


    def set_ittime(self,float t):
        """Set iteration time

        :param t: (float) :iteration time
        """
        self.ittime=t
