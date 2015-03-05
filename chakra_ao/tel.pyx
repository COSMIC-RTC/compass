
#################################################
# P-Class (parametres) param_tel
#################################################
cdef class param_tel:
    def __cinit__(self):
        self.diam=0
        self.type_ap="Generic"
        self.t_spiders=0.01
        self.spiders_type="six"
        self.nbrmissing=0
        self.referr=0

    def set_diam(self,float d):
        """set attribute diam to d

        d -- float : telescope diameter (in meters)."""
        self.diam=d

    def set_cobs(self,float c):
        """set attribute cobs to c

        c -- float : central obstruction ratio."""
        self.cobs=c

    def set_type_ap(self,str t):
        """set attribute type_ap to t

        t -- str : EELT aperture type."""
        self.type_ap=t

    def set_t_spiders(self,float spider):
        """set attribute t_spiders to spider

        spider -- float : secondary supports ratio."""
        self.t_spiders=spider

    def set_spiders_type(self, str spider):
        """set attribute spiders_type to spider

        spider -- str : secondary supports type."""
        self.spiders_type= spider

    def set_pupangle(self,float p):
        """set attribute pupangle to p

        p -- float : rotation angle of pupil."""
        self.pupangle=p

    def set_nbrmissing(self,long nb):
        """set attribute nbrmissing to nb

        nb -- long : number of missing segments for EELT pupil (max is 20)."""
        self.nbrmissing=nb

    def set_referr(self, float ref):
        """set attribute referr to ref

        ref -- float : std of reflectivity errors for EELT segments (fraction)."""
        self.referr=ref
