import make_pupil as mkP

#################################################
# P-Class (parametres) param_geom
#################################################
cdef class param_geom:
    #def __cinit__(self):


    def geom_init(self, param_tel tel, long pupdiam):
        self.pupdiam=pupdiam
        #first poxer of 2 greater than pupdiam
        self.ssize=long(2**np.ceil(np.log2(pupdiam)+1))
        # using images centered on 1/2 pixels
        self.cent=self.ssize/2+0.5
        # valid pupil geometry
        self._p1=long(np.ceil(self.cent-pupdiam/2.))
        self._p2=long(np.floor(self.cent+pupdiam/2.))
        self.pupdiam=self._p2-self._p1+1
        self._n = self.pupdiam+4
        self._n1= self._p1-2
        self._n2= self._p2+2
       

        #useful pupil

        self._spupil=mkP.make_pupil(self.pupdiam,self.pupdiam,tel)


        # large pupil (used for image formation)
        self._ipupil=mkP.pad_array(self._spupil,self.ssize)
        
        # useful pupil + 4 pixels
        self._mpupil=mkP.pad_array(self._spupil,self._n)
        """TODO
        from yorick, func geom_init, file yoga_ao/yorick/yoga_ao.i, l 102
        if (y_target.apod==1):
            apodizer = float(make_apodizer(y_geom.pupdiam,y_geom.pupdiam,"/root/compass/trunk/yoga_ao/data/apodizer/SP_HARMONI_I4_C6_N1024.fits",180/12.))
            y_geom._apodizer = &apodizer

        """    

    def set_ssize( self,long s):
        """Set the attribute ssize to s
         s -- long : linear size of full image (in pixels)."""
        self.ssize= s

    def set_zenithangle( self, float z):
        """Set the attribute zenithangle to z
         z -- float : observations zenith angle (in deg)."""
        self.zenithangle=z
 
    def set_pupdiam( self,long p ):
        """Set the attribute pupdiam to p
         p -- long : linear size of total pupil (in pixels)."""
        self.pupdiam=p
 
    def set_cent( self, float c ):
        """Set the attribute cent to c
         c -- float : central point of the simulation."""
        self.cent=c


    def get_ipupil(self):
        """return the full pupil support"""
        return self._ipupil

    def get_mpupil(self):
        """return the padded pupil"""
        return self._mpupil

    def get_spupil(self):
        """return the pupil"""
        return self._spupil


    def get_n(self):
        return self._n

    def get_n1(self):
        return self._n1

    def get_n2(self):
        return self._n2

    def get_p1(self):
        return self._p1

    def get_p2(self):
        return self._p2
