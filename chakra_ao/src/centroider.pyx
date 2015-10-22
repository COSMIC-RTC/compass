cdef class Param_centroider:
    def set_type(self, bytes t):
        """Set the centroider type
        :paramers:
            t: (str) : centroider type
        """
        self.type_centro=t

    def set_type_fct(self, bytes f):
        """Set the type of ref function

        :parameters:
            f: (str) : type of ref function
        """
        self.type_fct=f

    def set_nwfs(self, long n):
        """Set the index of wfs

        :parameters:
            n: (int) : index of wfs
        """
        self.nwfs=n

    def set_nmax(self, long n):
        """Set the number of brightest pixels to use for bpcog
        
        :parameters:
            n: (int) : number of brightest pixels
        """ 
        self.nmax=n

    def set_thresh(self, float t):
        """Set the threshold for tcog

        :parameters:
            t: (float) : threshold
        """
        self.thresh=t

    def set_width(self, float w):
        """Set the width of the Gaussian

        :parameters:
            w: (float) : width of the gaussian
        """
        self.width=w

    def set_sizex(self, long s):
        """Set sizex parameters for corr centroider (interp_mat size)

        :parameters:
            s: (long) : x size
        """
        self.sizex=s

    def set_sizey(self, long s):
        """Set sizey parameters for corr centroider (interp_mat size)

        :parameters:
            s: (long) : y size
        """
        self.sizey=s

    def set_weights(self,np.ndarray[ndim=3 ,dtype=np.float32_t] w):
        """Set the weights to use with wcog or corr

        :parameters:
            w: (np.ndarray[ndim=3 ,dtype=np.float32_t]) : weights
        """
        self.weights=w
