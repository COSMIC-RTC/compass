
#################################################
# P-Class (parametres) param_tel
#################################################
cdef class Param_tel:
    def __cinit__(self):
        self.diam=0
        self.type_ap="Generic"
        self.t_spiders=-1
        self.spiders_type=None
        self.nbrmissing=0
        self.referr=0

    def set_diam(self,float d):
        """set the telescope diameter

        :param d: (float) : telescope diameter (in meters)
        """
        self.diam=d

    def set_cobs(self,float c):
        """set the central obstruction ratio

        :param c: (float) : central obstruction ratio
        """
        self.cobs=c

    def set_type_ap(self,str t):
        """set the EELT aperture type

        :param t: (str) : EELT aperture type
        """
        self.type_ap=t

    def set_t_spiders(self,float spider):
        """set the secondary supports ratio

        :param spider: (float) : secondary supports ratio
        """
        self.t_spiders=spider

    def set_spiders_type(self, str spider):
        """set the secondary supports type

        :param spider: (str) : secondary supports type
        """
        self.spiders_type= spider

    def set_pupangle(self,float p):
        """set the rotation angle of pupil

        :param p: (float) : rotation angle of pupil
        """
        self.pupangle=p

    def set_nbrmissing(self,long nb):
        """set the number of missing segments for EELT pupil

        :param nb: (long) : number of missing segments for EELT pupil (max is 20)
        """
        self.nbrmissing=nb

    def set_referr(self, float ref):
        """set the std of reflectivity errors for EELT segments

        :param ref: (float) : std of reflectivity errors for EELT segments (fraction)
        """
        self.referr=ref
