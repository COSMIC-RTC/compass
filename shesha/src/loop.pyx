cdef class Param_loop:
    def set_niter(self,long n):
        """Set the number of iteration
        
        :parameters:
            n: (long) : number of iteration
        """
        self.niter=n


    def set_ittime(self,float t):
        """Set iteration time

        :parameters:
            t: (float) :iteration time
        """
        self.ittime=t
