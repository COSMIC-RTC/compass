cdef class Param_centroider:
    def set_type(self, bytes t):
        """Set the centroider type
        :param t: (str) : centroider type
        """
        self.type_centro=t

    def set_type_fct(self, bytes f):
        """Set the type of ref function

        :param f: (str) : type of ref function
        """
        self.type_fct=f

    def set_nwfs(self, long n):
        """Set the index of wfs

        :param n: (int) : index of wfs
        """
        self.nwfs=n

    def set_nmax(self, long n):
        """TODO doc
        
        :param n: (int) :
        """ 
        self.nmax=n

    def set_thresh(self, float t):
        """TODO doc

        :param t: (float) :
        """
        self.thresh=t

    def set_width(self, float w):
        """Set the width of the Gaussian

        :param w: (float) : width of the gaussian
        """
        self.width=w

    def set_sizex(self, long s):
        """TODO doc

        :param s: (long) :
        """
        self.sizex=s

    def set_sizey(self, long s):
        """TODO doc

        :param s: (long) :
        """
        self.sizey=s

    def set_weights(self,np.ndarray[ndim=3 ,dtype=np.float32_t] w):
        """TODO doc

        :param w: (np.ndarray[ndim=3 ,dtype=np.float32_t]) : 
        """
        self.weights=w
