cdef class Param_controller:
    def set_type(self,bytes b):
        self.type_control=b

    def set_nwfs(self,l):
        """Set the indices of wfs seen by the controller

        :parameters:
            l: (list of int) : indices of wfs
        """
        self.nwfs=np.array(l,dtype=np.int32)

    def set_nvalid(self,l):
        """Set the number of valid subaps

        :parameters:
            l: (list of int) : number of valid subaps per wfs
        """
        self.nvalid=np.array(l,dtype=np.int32)

    def set_ndm(self,l):
        """Set the indices of dms controlled

        :parameters:
            l: (list of int) : indices of dms
        """
        self.ndm=np.array(l,dtype=np.int32)

    def set_nactu(self,l):
        """Set the number of controlled actuators

        :parameters:
            l: (list of int) : number of controled actuator per dm
        """
        self.nactu=np.array(l,dtype=np.int32)

    def set_maxcond(self,float m):
        """Set the conditioning for imat inversion

        :parameters:
            : (float) : max condition number
        """
        self.maxcond=m

    def set_delay(self,float d):
        """Set the loop delay expressed in frames

        :parameters:
            d: (float) :delay [frames]
        """
        self.delay=d

    def set_gain(self,float g):
        """Set the loop gain 

        :parameters:
            g: (float) : loop gain
        """
        self.gain=g

    def set_nkl(self, long n):
        """Set the number of KL modes used for computation of covmat in case of minimum variance controller

        :parameters:
            n: (long) : number of KL modes
        """
        self.nkl=n

    def set_cured_ndivs(self,long c):
        """Set the subdivision levels in cured

        :parameters:
            c: (long) : subdivision levels in cured
        """
        self.cured_ndivs=c

    def set_modopti(self, int m):
        """Set the flag for modal optimization

        :parameters:
            m: (int) : flag for modal optimization
        """
        self.modopti=m

    def set_nrec(self, int n):
        """Set the number of sample of open loop slopes for modal optimization computation

        :parameters:
            n: (int) : number of sample
        """
        self.nrec=n

    def set_nmodes(self, int n):
        """Set the number of modes for M2V matrix (modal optimization)

        :parameters:
            n: (int) : number of modes
        """
        self.nmodes=n

    def set_gmin(self, float g):
        """Set the minimum gain for modal optimization

        :parameters:
            g: (float) : minimum gain for modal optimization
        """
        self.gmin=g

    def set_gmax(self, float g):
        """Set the maximum gain for modal optimization

        :parameters:
            g: (flaot) : maximum gain for modal optimization
        """
        self.gmax=g

    def set_ngain(self, int n):
        """Set the number of tested gains

        :parameters:
            n: (int) : number of tested gains
        """
        self.ngain=n

    def set_imat(self,np.ndarray[ndim=2,dtype=np.float32_t] imat):
        """Set the full interaction matrix

        :parameters:
            imat: (np.ndarray[ndim=2,dtype=np.float32_t]) : full interaction matrix
        """
        self.imat=imat

    def set_cmat(self,np.ndarray[ndim=2,dtype=np.float32_t] cmat):
        """Set the full control matrix

        :parameters:
            cmat: (np.ndarray[ndim=2,dtype=np.float32_t]) : full control matrix
        """
        self.cmat=cmat

