cdef class Param_controller:
    def set_type(self,bytes b):
        self.type_control=b
    def set_nwfs(self,l):
        self.nwfs=np.array(l,dtype=np.int32)
    def set_nvalid(self,l):
        self.nvalid=np.array(l,dtype=np.int32)
    def set_ndm(self,l):
        self.ndm=np.array(l,dtype=np.int32)
    def set_nactu(self,l):
        self.nactu=np.array(l,dtype=np.int32)
    def set_maxcond(self,float m):
        self.maxcond=m
    def set_delay(self,float d):
        self.delay=d
    def set_gain(self,float g):
        self.gain=g
    def set_nkl(self, long n):
        self.nkl=n
    def set_cured_ndivs(self,long c):
        self.cured_ndivs=c
    def set_modopti(self, int m):
        self.modopti=m
    def set_nrec(self, int n):
        self.nrec=n
    def set_nmodes(self, int n):
        self.nmodes=n
    def set_gmin(self, float g):
        self.gmin=g
    def set_gmax(self, float g):
        self.gmax=g
    def set_ngain(self, int n):
        self.ngain=n
    def set_imat(self,np.ndarray[ndim=2,dtype=np.float32_t] imat):
        self.imat=imat
    def set_cmat(self,np.ndarray[ndim=2,dtype=np.float32_t] cmat):
        self.cmat=cmat

