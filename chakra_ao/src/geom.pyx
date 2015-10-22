
#################################################
# P-Class (parametres) Param_geom
#################################################
cdef class Param_geom:
    #def __cinit__(self):


    def geom_init(self, Param_tel tel, long pupdiam,apod):
        """Initialize the geometry of the system

        :parameters:
            tel: (Param_tel) : telescope settings

            pupdiam: (long) : linear size of total pupil 

            apod: (int) : apodizer
        """
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
      
        #TODO check filename 
        cdef bytes filename = <bytes> chakra_ao_savepath+\
                              <bytes>"apodizer/SP_HARMONI_I4_C6_N1024.npy"

        cdef float cent=self.pupdiam/2.+0.5

        #useful pupil
        self._spupil=mkP.make_pupil(self.pupdiam,self.pupdiam,tel,cent,cent)
        
        # large pupil (used for image formation)
        self._ipupil=mkP.pad_array(self._spupil,self.ssize).astype(np.float32)
        
        # useful pupil + 4 pixels
        self._mpupil=mkP.pad_array(self._spupil,self._n).astype(np.float32)
        if(apod==1):
            self.apodizer=make_apodizer(self.pupdiam,self.pupdiam,filename,180./12.)


    def set_ssize( self,long s):
        """Set linear size of full image

         :parameters:
             s: (long) : linear size of full image (in pixels)."""
        self.ssize= s

    def set_zenithangle( self, float z):
        """Set observations zenith angle 

         :parameters:
             z: (float) : observations zenith angle (in deg)."""
        self.zenithangle=z
 
    def set_pupdiam( self,long p ):
        """Set the linear size of total pupil

        :parameters:
            p: (long) : linear size of total pupil (in pixels)."""
        self.pupdiam=p
 
    def set_cent( self, float c ):
        """Set the central point of the simulation

         :parameters:
             c: (float) : central point of the simulation."""
        self.cent=c


    def get_ipupil(self):
        """return the full pupil support"""
        return self._ipupil

    def get_mpupil(self):
        """return the padded pupil"""
        return self._mpupil

    def get_spupil(self):
        """return the small pupil"""
        return self._spupil


    def get_n(self):
        """Return the linear size of the medium pupil"""
        return self._n

    def get_n1(self):
        """Return the min(x,y) for valid points for the total pupil"""
        return self._n1

    def get_n2(self):
        """Return the max(x,y) for valid points for the total pupil"""
        return self._n2

    def get_p1(self):
        """Return the min(x,y) for valid points for the medium pupil"""
        return self._p1

    def get_p2(self):
        """Return the max(x,y) for valid points for the medium pupil"""
        return self._p2

