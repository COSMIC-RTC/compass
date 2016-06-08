import os
import sys
import inspect

import numpy as np
cimport numpy as np
np.import_array()

shesha_dir = os.environ.get('SHESHA_ROOT')
if(shesha_dir is None):
    raise EnvironmentError("Environment variable 'SHESHA_ROOT' must be define")
shesha_savepath = < bytes > (shesha_dir + "/data/")
print "shesha_savepath:", shesha_savepath

sys.path.append(shesha_dir + '/src')
import make_pupil as mkP

from scipy.ndimage import interpolation as interp

RASC = 180.*3600. / np.pi

#################################################
# P-Class (parametres) Param_loop
#################################################
cdef class Param_loop:
    def set_niter(self, long n):
        """Set the number of iteration

        :parameters:
            n: (long) : number of iteration
        """
        self.niter = n


    def set_ittime(self, float t):
        """Set iteration time

        :parameters:
            t: (float) :iteration time
        """
        self.ittime = t



#################################################
# P-Class (parametres) param_tel
#################################################
cdef class Param_tel:
    def __cinit__(self):
        self.diam = 0
        self.type_ap = "Generic"
        self.t_spiders = -1
        self.spiders_type = None
        self.nbrmissing = 0
        self.referr = 0
        self.std_piston = 0
        self.std_tt = 0

    def set_diam(self, float d):
        """set the telescope diameter

        :param d: (float) : telescope diameter (in meters)
        """
        self.diam = d

    def set_cobs(self, float c):
        """set the central obstruction ratio

        :param c: (float) : central obstruction ratio
        """
        self.cobs = c

    def set_type_ap(self, str t):
        """set the EELT aperture type

        :param t: (str) : EELT aperture type
        """
        self.type_ap = t

    def set_t_spiders(self, float spider):
        """set the secondary supports ratio

        :param spider: (float) : secondary supports ratio
        """
        self.t_spiders = spider

    def set_spiders_type(self, str spider):
        """set the secondary supports type

        :param spider: (str) : secondary supports type
        """
        self.spiders_type = spider

    def set_pupangle(self, float p):
        """set the rotation angle of pupil

        :param p: (float) : rotation angle of pupil
        """
        self.pupangle = p

    def set_nbrmissing(self, long nb):
        """set the number of missing segments for EELT pupil

        :param nb: (long) : number of missing segments for EELT pupil (max is 20)
        """
        self.nbrmissing = nb

    def set_referr(self, float ref):
        """set the std of reflectivity errors for EELT segments

        :param ref: (float) : std of reflectivity errors for EELT segments (fraction)
        """
        self.referr = ref

    def set_std_piston(self, float piston):
        """set the std of piston errors for EELT segments

        :param piston: (float) : std of piston errors for EELT segments
        """
        self.std_piston = piston

    def set_std_tt(self, float tt):
        """set the std of tip-tilt errors for EELT segments

        :param tt: (float) : std of tip-tilt errors for EELT segments
        """
        self.std_tt = tt


#################################################
# P-Class (parametres) Param_geom
#################################################
cdef class Param_geom:

    def geom_init(self, Param_tel tel, long pupdiam, apod):
        """Initialize the system geometry

        :parameters:
            tel: (Param_tel) : telescope settings

            pupdiam: (long) : linear size of total pupil

            apod: (int) : apodizer
        """
        self.pupdiam = pupdiam
        # first poxer of 2 greater than pupdiam
        self.ssize = long(2 ** np.ceil(np.log2(pupdiam) + 1))
        # using images centered on 1/2 pixels
        self.cent = self.ssize / 2 + 0.5
        # valid pupil geometry
        self._p1 = long(np.ceil(self.cent - pupdiam / 2.))
        self._p2 = long(np.floor(self.cent + pupdiam / 2.))
        self.pupdiam = self._p2 - self._p1 + 1
        self._n = self.pupdiam + 4
        self._n1 = self._p1 - 2
        self._n2 = self._p2 + 2

        # TODO check filename
        cdef bytes filename = < bytes > shesha_savepath + \
                              < bytes > "apodizer/SP_HARMONI_I4_C6_N1024.npy"

        cdef float cent = self.pupdiam / 2. + 0.5

        # useful pupil
        self._spupil = mkP.make_pupil(self.pupdiam, self.pupdiam, tel, cent, cent).astype(np.float32)

        self._phase_ab_M1 = mkP.make_phase_ab(self.pupdiam, self.pupdiam, tel, self._spupil).astype(np.float32)

        # large pupil (used for image formation)
        self._ipupil = mkP.pad_array(self._spupil, self.ssize).astype(np.float32)

        # useful pupil + 4 pixels
        self._mpupil = mkP.pad_array(self._spupil, self._n).astype(np.float32)
        self._phase_ab_M1_m = mkP.pad_array(self._phase_ab_M1, self._n).astype(np.float32)

        if(apod == 1):
            self._apodizer = make_apodizer(self.pupdiam, self.pupdiam, filename, 180. / 12.).astype(np.float32)
        else:
            self._apodizer = np.ones((self._spupil.shape[0], self._spupil.shape[1])).astype(np.float32)



    def set_ssize(self, long s):
        """Set linear size of full image

         :param s: (long) : linear size of full image (in pixels)."""
        self.ssize = s

    def set_zenithangle(self, float z):
        """Set observations zenith angle

         :param z: (float) : observations zenith angle (in deg)."""
        self.zenithangle = z

    def set_pupdiam(self, long p):
        """Set the linear size of total pupil

        :param p: (long) : linear size of total pupil (in pixels)."""
        self.pupdiam = p

    def set_cent(self, float c):
        """Set the central point of the simulation

         :param c: (float) : central point of the simulation."""
        self.cent = c


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



#################################################
# P-Class (parametres) Param_wfs
#################################################
cdef class Param_wfs:

    def __cinit__(self, bool error_budget=False):
        self.error_budget = error_budget
        self.nphotons4imat = 1.e5

    def set_type(self, str t):
        """Set the type of wfs

        :param t: (str) : type of wfs ("sh" or "pyr")
        """
        self.type_wfs = t

    def set_nxsub(self, long n):
        """Set the linear number of subaps

        :param n: (long) : linear number of subaps
        """
        self.nxsub = n

    def set_npix(self, long  n):
        """Set the number of pixels per subap

        :param n: (long) : number of pixels per subap
        """
        self.npix = n

    def set_pixsize(self, float p):
        """Set the pixel size

        :param p: (float) : pixel size (in arcsec) for a subap
        """
        self.pixsize = p

    def set_Lambda(self, float L):
        """Set the observation wavelength

        :param L: (float) : observation wavelength (in ��m) for a subap
        """
        self.Lambda = L

    def set_optthroughput(self, float o):
        """Set the wfs global throughput

        :param o: (float) : wfs global throughput
        """
        self.optthroughput = o

    def set_fracsub(self, float f):
        """Set the minimal illumination fraction for valid subaps

        :param f: (float) : minimal illumination fraction for valid subaps
        """
        self.fracsub = f

    def set_openloop(self, long o):
        """Set the loop state (open or closed)

        :param o: (long) : 1 if in "open-loop" mode (i.e. does not see dm)
        """
        self.openloop = o

    def set_fssize(self, float f):
        """Set the size of field stop

        :param f: (float) : size of field stop in arcsec
        """
        self.fssize = f

    def set_fstop(self,str f):
        """Set the size of field stop

        :param f: (str) : size of field stop in arcsec
        """
        self.fstop = f

    def set_atmos_seen(self, int i):
        """Tells if the wfs sees the atmosphere layers

        :param i: (int) :1 if the WFS sees the atmosphere layers
        """
        self.atmos_seen = i

    def set_xpos(self, float x):
        """Set the guide star x position on sky

        :param x: (float) : guide star x position on sky (in arcsec)
        """
        self.xpos = x

    def set_ypos(self, float y):
        """Set the guide star y position on sky

        :param y: (float) : guide star y position on sky (in arcsec)
        """
        self.ypos = y

    def set_gsalt(self, float g):
        """Set the altitude of guide star

        :param g: (float) : altitude of guide star (in m) 0 if ngs
        """
        self.gsalt = g

    def set_gsmag(self, float g):
        """Set the magnitude of guide star

        :param g: (float) : magnitude of guide star
        """
        self.gsmag = g

    def set_zerop(self, float z):
        """Set the detector zero point

        :param z: (float) : detector zero point
        """
        self.zerop = z

    def set_noise(self, float n):
        """Set the desired noise

        :param n: (float) : desired noise : < 0 = no noise / 0 = photon only / > 0 photon + ron
        """
        self.noise = n

    def set_nphotons4imat(self, float nphot):
        """Set the desired numner of photons used for doing imat

        :param nphot: (float) : desired number of photons
        """
        self.nphotons4imat = nphot

    def set_kernel(self, float k):
        """Set the attribute kernel

        :param k: (float) :
        """
        self.kernel = k

    def set_laserpower(self, float l):
        """Set the laser power

        :param l: (float) : laser power in W
        """
        self.laserpower = l

    def set_lltx(self, float l):
        """Set the x position of llt

        :param l: (float) : x position (in meters) of llt
        """
        self.lltx = l

    def set_llty(self, float l):
        """Set the y position of llt

        :param l: (float) : y position (in meters) of llt
        """
        self.llty = l

    def set_proftype(self, str p):
        """Set the type of sodium profile

        :param p: (str) : type of sodium profile "gauss", "exp", etc ...
        """
        self.proftype = p

    def set_beamsize(self, float b):
        """Set the laser beam fwhm on-sky

        :param b: (float) : laser beam fwhm on-sky (in arcsec)
        """
        self.beamsize = b

    def set_pyr_ampl(self, float p):
        """Set the pyramid wfs modulation amplitude radius

        :param p: (float) : pyramid wfs modulation amplitude radius (in arsec)
        """
        self.pyr_ampl = p

    def set_pyr_npts(self, long p):
        """Set the total number of point along modulation circle

        :param p: (long) : total number of point along modulation circle
        """
        self.pyr_npts = p

    def set_pyr_loc(self, str p):
        """Set the location of modulation

        :param p: (str) : location of modulation, before/after the field stop.
                          valid value are "before" or "after" (default "after")
        """
        self.pyr_loc = p

    def set_pyrtype(self, str p):
        """Set the type of pyramid,

        :param p: (str) : type of pyramid, either 0 for "Pyramid" or 1 for "RoofPrism"
        """
        self.pyrtype = p

    def set_dms_seen(self, np.ndarray[ndim=1, dtype=np.int32_t] dms_seen):
        """Set the index of dms seen by the WFS

        :param dms_seen: (np.ndarray[ndim=1,dtype=np.int32_t) : index of dms seen by the WFS
        """
        self.dms_seen = dms_seen

    def set_lgsreturnperwatt(self, float lpw):
        """Set the return per watt factor

        :param lpw: (float) : return per watt factor (high season : 10 ph/cm2/s/W)
        """
        self.lgsreturnperwatt = lpw


    def set_altna(self, np.ndarray[ndim=1, dtype=np.float32_t] a):
        """Set the corresponding altitude

        :param a: (np.ndarray[ndim=1,dtype=np.float32]) : corresponding altitude
        """
        self._altna = a


    def set_profna(self, np.ndarray[ndim=1, dtype=np.float32_t] p):
        """Set the sodium profile

        :param p: (np.ndarray[ndim=1,dtype=np.float32]) : sodium profile
        """
        self._profna = p

    def set_errorBudget(self, bool error_budget):
        """ Set the error budget flag : if True, enable error budget analysis
        for this simulation

        :param error_budget: (bool) : error budget flag
        """
        self.error_budget = error_budget

    def get_validsub(self):
        return self._validsubsx, self._validsubsy



'''
    cdef make_lgs_prof1d(self, Param_tel p_tel,
            np.ndarray[dtype=np.float32_t] prof, np.ndarray[dtype=np.float32_t] h,
            float beam, bytes center=<bytes>""):
        """same as prep_lgs_prof but cpu only. original routine from rico

        :parameters:
            p_tel: (Param_tel) : telescope settings

            prof: (np.ndarray[dtype=np.float32]) : Na profile intensity, in arbitrary units

            h: (np.ndarray[dtype=np.float32]) : altitude, in meters. h MUST be an array with EQUALLY spaced elements.

            beam: (float) : size in arcsec of the laser beam

            center: (string) : either "image" or "fourier" depending on where the centre should be.
        """

        self._prof1d=prof
        self._profcum=np.zeros(prof.shape[1]+1,dtype=np.float32)
        self._profcum[1:]=prof.cumsum()


        cdef float subapdiam=p_tel.diam/self.nxsub #diam of subap
        cdef np.ndarray[dtype=np.float32_t] xsubs, ysubs,dOffAxis
        cdef np.ndarray[dtype=np.float32_t] g
        if(self.nxsub>1):
            xsubs=np.linspace( (subapdiam-p_tel.diam)/2,(p_tel.diam-subapdiam)/2,
                                self.nxsub).astype(np.float32)
        else:
            xsubs=np.zeros(1,dtype=np.float32)
        ysubs=xsubs.copy().astype(np.float32)


        #cdef int nP=prof.shape[0] #UNUSED
        cdef float hG=np.sum(h*prof)/np.sum(prof)
        cdef np.ndarray[dtype=np.float32_t] x=np.arange(self._Ntot).astype(np.float32)-self._Ntot/2
        # x expressed in pixels. (0,0) is in the fourier-center.
        x=x*self._qpixsize #x expressed in arcseconds
        #cdef float dx=x[1]-x[0] #UNUSED
        #cdef float dh=h[1]-h[0] #UNUSED

        if(self.nxsub>1):
            dOffAxis=np.sqrt((xsubs[self._validsubsy]-self.lltx)**2+
                             (ysubs[self._validsubsx]-self.llty)**2)
        else:
            dOffAxis=np.sqrt((xsubs-self.lltx)**2+(ysubs-self.llty)**2)

        cdef np.ndarray[ndim=2,dtype=np.float32_t] profi,
        cdef np.ndarray[ndim=1,dtype=np.float32_t] zhc, avg_zhc,avg_x
        profi=np.zeros((self._Ntot,self._nvalid),dtype=np.float32)

        #cdef np.ndarray[ndim=1,dtype=np.int32_t] subsdone, dif2do #UNUSED
        cdef np.ndarray[ndim=1,dtype=np.int64_t] inds
        subsdone=np.ones(self._nvalid,dtype=np.int32)
        #dif2do =np.zeros(self._nvalid,dtype=np.int32) #UNUSED
        cdef float tmp,dzhc

        cdef int i

        zhc=(h.astype(np.float32)-hG)*(206265.*tmp/hG**2)
        dzhc = zhc[1]-zhc[0]

        tmp=dOffAxis[np.where(subsdone)][0]
        zhc=(h-hG)*(206265.*tmp/hG**2)
        avg_zhc=np.zeros(zhc.size+1,dtype=np.float32)
        avg_zhc[0]=zhc[0]
        avg_zhc[avg_zhc.size-1]=zhc[zhc.size-1]
        avg_zhc[1:-1]=0.5*(zhc[1:]+zhc[:-1])
        avg_x=np.zeros(x.size+1,dtype=np.float32)
        avg_x[0]=x[0]
        avg_x[avg_x.size-1]=x[x.size-1]
        avg_x[1:-1]=0.5*(x[1:]+x[:-1])


        while(np.any(subsdone)):
            tmp=dOffAxis[np.where(subsdone)][0]
            inds=np.where(dOffAxis==tmp)[0]
            # height, translated in arcsec due to perspective effect
            zhc=(h-hG)*(206265.*tmp/hG**2)
            dzhc = zhc[1]-zhc[0]

            if(self._qpixsize>dzhc):
                avg_zhc=np.zeros(zhc.size+1,dtype=np.float32)
                avg_zhc[0]=zhc[0]
                avg_zhc[avg_zhc.size-1]=zhc[zhc.size-1]
                avg_zhc[1:-1]=0.5*(zhc[1:]+zhc[:-1])
                avg_x=np.zeros(x.size+1,dtype=np.float32)
                avg_x[0]=x[0]
                avg_x[avg_x.size-1]=x[x.size-1]
                avg_x[1:-1]=0.5*(x[1:]+x[:-1])

                for i in range(inds.size):
                    profi[:,inds[i]]=np.diff(np.interp(avg_x,avg_zhc,self._profcum)).astype(np.float32)

            else:
                for i in range(inds.size):
                    profi[:,inds[i]]=np.interp(x,zhc,prof)
            subsdone[inds]=0



        cdef w=beam/2.35482005
        if(w==0):
            #TODO what is n
            n=1
            g=np.zeros(n,dtype=np.float32)
            if(center=="image"):
                g[n/2-1]=0.5
                g[n/2]=0.5
            else:
                g[n/2]=1

        else:
            if(center=="image"):
                if( (self.npix*self._nrebin)%2 != self._Nfft%2):
                    g=np.exp(-(x+self._qpixsize)**2/(2*w**2.))
                else:
                    g=np.exp(-(x+self._qpixsize/2)**2/(2*w**2.))

            else:
                g=np.exp(-x**2/(2*w**2.))


        self._ftbeam=np.fft.fft(g).astype(np.complex64)
        self._beam=g.astype(np.float32)
        #convolved profile in 1D.

        cdef np.ndarray[ndim=2,dtype=np.float32_t] p1d
        cdef np.ndarray[ndim=2,dtype=np.float32_t] g_extended = np.tile(g,(self._nvalid,1)).T


        p1d=np.fft.ifft( np.fft.fft(profi,axis=0)*
               np.fft.fft(g_extended,axis=0) ,
              axis=0).real.astype(np.float32)
        p1d=p1d*p1d.shape[0]
        p1d=np.roll(p1d,int(self._Ntot/2.+0.5),axis=0)
        p1d=np.abs(p1d)


        cdef np.ndarray[ndim=3,dtype=np.float32_t] im
        im=np.zeros((p1d.shape[1],p1d.shape[0],p1d.shape[0]),dtype=np.float32)


        cdef int l,c
        for i in range(p1d.shape[1]):
            for l in range(p1d.shape[0]):
                for c in range(p1d.shape[0]):
                    im[i,l,c]=g[l]*p1d[c,i]

        if(ysubs.size>1):
            azimuth=np.arctan2(ysubs[self._validsubsy]-self.llty,
                                xsubs[self._validsubsx]-self.lltx)
        else:
            azimuth=np.arctan2(ysubs-self.llty,
                                xsubs-self.lltx)

        self._azimuth=azimuth

        cdef float xcent,ycent

        if(center=="image"):
            xcent=self._Ntot/2.+0.5
            ycent=xcent
        else:
            xcent=self._Ntot/2.+1
            ycent=xcent

        cdef np.ndarray[ndim=1,dtype=np.float32_t] max_im
        if(ysubs.size>0):
        #TODO rotate
            im=rotate3d(im,azimuth*180/np.pi,xcent,ycent)
            max_im=np.max(im,axis=(1,2))
            im=(im.T/max_im).T
        else:
            im=rotate(im,azimuth*180/np.pi,xcent,ycent)
            im=im/np.max(im)

        self._lgskern=im.T
'''


#################################################
# P-Class (parametres) Param_atmos
#################################################
cdef class Param_atmos:
    def __cinit__(self):
        self.nscreens = 0

    def set_nscreens(self, long n):
        """Set the number of turbulent layers

        :param n: (long) number of screens.
        """
        self.nscreens = n

    def set_r0(self, float r):
        """Set the global r0

        :param r: (float) : global r0
        """
        self.r0 = r

    def set_pupixsize(self, float xsize):
        """Set the pupil pixel size

        :param xsize: (float) : pupil pixel size
        """
        self.pupixsize = xsize

    def set_L0(self, l):
        """Set the L0 per layers

        :param l: (lit of float) : L0 for each layers
        """
        self.L0 = np.array(l, dtype=np.float32)

    def set_dim_screens(self, l):
        """Set the size of the phase screens

        :param l: (lit of float) : phase screens sizes
        """
        self.dim_screens = np.array(l, dtype=np.float32)

    def set_alt(self, l):
        """Set the altitudes of each layer

        :param l: (lit of float) : altitudes
        """
        self.alt = np.array(l, dtype=np.float32)

    def set_winddir(self, l):
        """Set the wind direction for each layer

        :param l: (lit of float) : wind directions
        """
        self.winddir = np.array(l, dtype=np.float32)

    def set_windspeed(self, l):
        """Set the the wind speed for each layer

        :param l: (lit of float) : wind speeds
        """
        self.windspeed = np.array(l, dtype=np.float32)

    def set_frac(self, l):
        """Set the fraction of r0 for each layers

        :param l: (lit of float) : fraction of r0
        """
        self.frac = np.array(l, dtype=np.float32)

    def set_deltax(self, l):
        """Set the translation speed on axis x for each layer

        :param l: (lit of float) : translation speed
        """
        self.deltax = np.array(l, dtype=np.float32)

    def set_deltay(self, l):
        """Set the translation speed on axis y for each layer

        :param l: (lit of float) : translation speed
        """
        self.deltay = np.array(l, dtype=np.float32)

    def set_seeds(self, l):
        """Set the seed for each layer

        :param l: (lit of int) : seed
        """
        self.seeds = np.array(l, dtype=np.int64)



#################################################
# P-Class (parametres) Param_dm
#################################################
cdef class Param_dm:

    def set_type(self, bytes t):
        """set the dm type

        :param t: (str) : type of dm
        """
        self.type_dm = t
        
    def set_pattern(self,bytes t):
        """set the pattern type

        :param t: (str) : type of pattern
        """
        self.type_pattern=t
        
    def set_file_influ_hdf5(self,bytes f):
        """set the name of hdf5 influence file 

        :param filename: (str) : Hdf5 file influence name
        """
        self.file_influ_hdf5=f
        
    def set_center_name(self,bytes f):
        """set the name of hdf5 influence file 

        :param filename: (str) : Hdf5 file influence name
        """
        self.center_name=f

    def set_cube_name(self,bytes cubename):
        """set the name of influence cube in hdf5 

        :param cubename: (str) : name of influence cube
        """
        self.cube_name=cubename
        
    def set_x_name(self,bytes xname):
        """set the name of x coord of influence fonction in file 

        :param t: (str) : name of x coord of influence
        """
        self.x_name=xname

    def set_y_name(self,bytes yname):
        """set the name of y coord of influence fonction in file 

        :param yname: (str) : name of y coord of influence
        """
        self.y_name=yname
    
    def set_nact(self, long n):

        """set the number of actuator

        :param n: (long) : number of actuators in the diameter
        """
        self.nact = n

    def set_margin(self,float n):
        """set the margin for outside actuator select

        :param n: (float) : pupille diametre ratio for actuator select 
        """
        self.margin=n  
        
    def set_margin_out(self,float n):
        """set the margin for outside actuator select

        :param n: (float) : pupille diametre ratio for actuator select 
        """
        self.margin_out=n       

    def set_margin_in(self,float n):
        """set the margin for inside actuator select (central obstruction)

        :param n: (float) : pupille diametre ratio for actuator select 
        """
        self.margin_in=n    

    def set_alt(self, float a):
        """set the conjugaison altitude

        :param a: (float) : conjugaison altitude (im m)
        """
        self.alt = a

    def set_thresh(self, float t):
        """set the threshold on response for selection

        :param t: (float) : threshold on response for selection (<1)
        """
        self.thresh = t

    def set_coupling(self, float c):
        """set the actuators coupling

        :param c: (float) : actuators coupling (<0.3)
        """
        self.coupling = c

    def set_unitpervolt(self, float u):
        """set the Influence function sensitivity

        :param u: (float) : Influence function sensitivity in unit/volt
        """
        self.unitpervolt = u

    def set_push4imat(self, p):
        """set the nominal voltage for imat

        :param p: (float) : nominal voltage for imat
        """
        self.push4imat = p

    def set_ntotact(self, long n):
        """set the total number of actuators

        :param n: (long) : total number of actuators
        """
        self._ntotact = n

    def set_xpos(self, np.ndarray[ndim=1, dtype=np.float32_t] xpos):
        """Set the x positions of influ functions (lower left corner)

        :param xpos: (np.ndarray[ndim=1,dtype=np.float32_t]) : x positions of influ functions
        """
        self._xpos = xpos

    def set_ypos(self, np.ndarray[ndim=1, dtype=np.float32_t] ypos):
        """Set the y positions of influ functions (lower left corner)

        :param ypos: (np.ndarray[ndim=1,dtype=np.float32_t]) : y positions of influ functions
        """
        self._ypos = ypos

    def set_i1(self, np.ndarray[ndim=1, dtype=np.int32_t] i1):
        """TODO doc

        :param i1: (np.ndarray[ndim=1,dtype=np.int32_t]) :
        """
        self._i1 = i1

    def set_j1(self, np.ndarray[ndim=1, dtype=np.int32_t] j1):
        """TODO doc

        :param j1: (np.ndarray[ndim=1,dtype=np.int32_t]) :
        """
        self._j1 = j1

    def set_influ(self, np.ndarray[ndim=3, dtype=np.float32_t] influ):
        """Set the influence function

        :param influ: (np.ndarray[ndim=3,dtype=np.float32_t]) : influence function
        """
        self._influ = influ



#################################################
# P-Class (parametres) Param_target
#################################################
cdef class Param_target:
    def __cinit__(self):
        self.ntargets = 0
        self.zerop = 1.


    def set_nTargets(self, int n):
        """Set the number of targets

        :param n: (int) : number of targets
        """
        self.ntargets = n
    def set_apod(self, int a):
        """Tells if the apodizer is used

        The apodizer is used if a is not 0
        :param a: (int) boolean for apodizer
        """
        self.apod = a

    def set_Lambda(self, l):
        """Set the observation wavelength

        :param l: (list of float) : observation wavelength for each target
        """
        self.Lambda = np.array(l, dtype=np.float32)

    def set_xpos(self, l):
        """Set the x positions on sky (in arcsec)

        :param l: (list of float) : x positions on sky for each target
        """
        self.xpos = np.array(l, dtype=np.float32)

    def set_ypos(self, l):
        """Set the y positions on sky (in arcsec)

        :param l: (list of float) : y positions on sky for each target
        """
        self.ypos = np.array(l, dtype=np.float32)

    def set_mag(self, l):
        """set the magnitude

        :param l: (list of float) : magnitude for each target
        """
        self.mag = np.array(l, dtype=np.float32)

    def set_dms_seen(self, l):
        """set the dms seen by the target

        :param l: (list of int) : index for each dm
        """
        self.dms_seen = np.array(l, dtype=np.int32)

    def set_zerop(self, float z):
        """Set the detector zero point

        :param z: (float) : detector zero point
        """
        self.zerop = z



#################################################
# P-Class (parametres) Param_rtc
#################################################
cdef class Param_rtc:

    def set_nwfs(self, n):
        """Set the number of wfs

        :param n: (int) number of wfs
        """
        self.nwfs = n
    def set_centroiders(self, l):
        """Set the centroiders

        :param l: (list of Param_centroider) : centroiders settings
        """
        self.centroiders = l

    def set_controllers(self, l):
        """Set the controller

        :param l: (list of Param_controller) : controllers settings
        """
        self.controllers = l



#################################################
# P-Class (parametres) Param_centroider
#################################################
cdef class Param_centroider:

    def set_type(self, bytes t):
        """Set the centroider type
        :param t: (str) : centroider type
        """
        self.type_centro = t

    def set_type_fct(self, bytes f):
        """Set the type of ref function

        :param f: (str) : type of ref function
        """
        self.type_fct = f

    def set_nwfs(self, long n):
        """Set the index of wfs

        :param n: (int) : index of wfs
        """
        self.nwfs = n

    def set_nmax(self, long n):
        """Set the number of brightest pixels to use for bpcog

        :parameters:
            n: (int) : number of brightest pixels
        """
        self.nmax = n

    def set_thresh(self, float t):
        """Set the threshold for tcog

        :parameters:
            t: (float) : threshold
        """
        self.thresh = t

    def set_width(self, float w):
        """Set the width of the Gaussian

        :param w: (float) : width of the gaussian
        """
        self.width = w

    def set_sizex(self, long s):
        """Set sizex parameters for corr centroider (interp_mat size)

        :parameters:
            s: (long) : x size
        """
        self.sizex = s

    def set_sizey(self, long s):
        """Set sizey parameters for corr centroider (interp_mat size)

        :parameters:
            s: (long) : y size
        """
        self.sizey = s

    def set_weights(self, np.ndarray[ndim=3 , dtype=np.float32_t] w):
        """Set the weights to use with wcog or corr

        :parameters:
            w: (np.ndarray[ndim=3 ,dtype=np.float32_t]) : weights
        """
        self.weights = w



#################################################
# P-Class (parametres) Param_controller
#################################################
cdef class Param_controller:

    def set_type(self, bytes b):
        self.type_control = b

    def set_nwfs(self, l):
        """Set the indices of wfs

        :param l: (list of int) : indices of wfs
        """
        self.nwfs = np.array(l, dtype=np.int32)

    def set_nvalid(self, l):
        """Set the number of valid subaps

        :param l: (list of int) : number of valid subaps per wfs
        """
        self.nvalid = np.array(l, dtype=np.int32)

    def set_ndm(self, l):
        """Set the indices of dms

        :param l: (list of int) : indices of dms
        """
        self.ndm = np.array(l, dtype=np.int32)

    def set_nactu(self, l):
        """Set the number of controled actuator

        :param l: (list of int) : number of controled actuator per dm
        """
        self.nactu = np.array(l, dtype=np.int32)

    def set_maxcond(self, float m):
        """Set the max condition number

        :param : (float) : max condition number
        """
        self.maxcond = m

    def set_TTcond(self, float m):
        """Set the tiptilt condition number for cmat filtering with mv controller

        :param : (float) : tiptilt condition number
        """
        self.TTcond = m

    def set_delay(self, float d):
        """Set the loop delay expressed in frames

        :parameters:
            d: (float) :delay [frames]
        """
        self.delay = d

    def set_gain(self, float g):
        """Set the loop gain

        :parameters:
            g: (float) : loop gain
        """
        self.gain = g

    def set_nkl(self, long n):
        """Set the number of KL modes used for computation of covmat in case of minimum variance controller

        :param n: (long) : number of KL modes
        """
        self.nkl = n

    def set_cured_ndivs(self, long c):
        """Set the subdivision levels in cured

        :param c: (long) : subdivision levels in cured
        """
        self.cured_ndivs = c

    def set_modopti(self, int m):
        """Set the flag for modal optimization

        :param m: (int) : flag for modal optimization
        """
        self.modopti = m

    def set_nrec(self, int n):
        """Set the number of sample of open loop slopes for modal optimization computation

        :param n: (int) : number of sample
        """
        self.nrec = n

    def set_nmodes(self, int n):
        """Set the number of modes for M2V matrix (modal optimization)

        :param n: (int) : number of modes
        """
        self.nmodes = n

    def set_gmin(self, float g):
        """Set the minimum gain for modal optimization

        :param g: (float) : minimum gain for modal optimization
        """
        self.gmin = g

    def set_gmax(self, float g):
        """Set the maximum gain for modal optimization

        :param g: (flaot) : maximum gain for modal optimization
        """
        self.gmax = g

    def set_ngain(self, int n):
        """Set the number of tested gains

        :param n: (int) : number of tested gains
        """
        self.ngain = n

    def set_imat(self, np.ndarray[ndim=2, dtype=np.float32_t] imat):
        """Set the full interaction matrix

        :param imat: (np.ndarray[ndim=2,dtype=np.float32_t]) : full interaction matrix
        """
        self.imat = imat

    def set_cmat(self, np.ndarray[ndim=2, dtype=np.float32_t] cmat):
        """Set the full control matrix

        :param cmat: (np.ndarray[ndim=2,dtype=np.float32_t]) : full control matrix
        """
        self.cmat = cmat






cpdef make_apodizer(int dim, int pupd, bytes filename, float angle):
    """TODO doc

    :parameters:
        (int) : im:

        (int) : pupd:

        (str) : filename:

        (float) : angle:
    """

    print "Opening apodizer"
    print "reading file:", filename
    cdef np.ndarray pup = np.load(filename)
    cdef int A = pup.shape[0]

    if(A > dim):
        raise ValueError("Apodizer dimensions must be smaller.")

    if (A != pupd):
        # use misc.imresize (with bilinear)
        print "TODO pup=bilinear(pup,pupd,pupd)"

    if (angle != 0):
        # use ndimage.interpolation.rotate
        print "TODO pup=rotate2(pup,angle)"
        pup = interp.rotate(pup, angle, reshape=False, order=2)

    reg = np.where(mkP.dist(pupd) > pupd / 2.)
    pup[reg] = 0.

    cdef np.ndarray[ndim = 2, dtype = np.float32_t] pupf = np.zeros((dim, dim), dtype=np.float32)

    if (dim != pupd):
        if ((dim - pupd) % 2 != 0):
            pupf[(dim - pupd + 1) / 2:(dim + pupd + 1) / 2, (dim - pupd + 1) / 2:(dim + pupd + 1) / 2] = pup

        else:
            pupf[(dim - pupd) / 2:(dim + pupd) / 2, (dim - pupd) / 2:(dim + pupd) / 2] = pup

    else:
        pupf = pup

    pupf = np.abs(pupf).astype(np.float32)

    return pupf


cpdef rotate3d(np.ndarray[ndim=3, dtype=np.float32_t] im,
              np.ndarray[ndim=1, dtype=np.float32_t] ang,
              float cx=-1, float cy=-1, float zoom=1.0):
    """Rotates an image of an angle "ang" (in DEGREES).

    The center of rotation is cx,cy.
    A zoom factor can be applied.

    (cx,cy) can be omitted :one will assume one rotates around the
    center of the image.
    If zoom is not specified, the default value of 1.0 is taken.

    modif dg : allow to rotate a cube of images with one angle per image

    :parameters:
        im: (np.ndarray[ndim=3,dtype=np.float32_t]) : array to rotate

        ang: (np.ndarray[ndim=1,dtype=np.float32_t]) : rotation angle  (in degrees)

        cx: (float) : (optional) rotation center on x axis (default: image center)

        cy: (float) : (optional) rotation center on x axis (default: image center)

        zoom: (float) : (opional) zoom factor (default =1.0)
"""

# TODO test it
    if(zoom == 0):
        zoom = 1.0
    if(ang.size == 1):
        if(zoom == 1.0 and ang[0] == 0.):
            return im

    ang *= np.pi / 180

    cdef int nx = im.shape[1]
    cdef int ny = im.shape[2]
    cdef np.ndarray[ndim = 2, dtype = np.float32_t] matrot

    cdef int i

    if(cx < 0):
        cx = nx / 2 + 1
    if(cy < 0):
        cy = ny / 2 + 1

    cdef np.ndarray x = np.tile(np.arange(nx) - cx + 1, (ny, 1)).T / zoom
    cdef np.ndarray y = np.tile(np.arange(ny) - cy + 1, (nx, 1)) / zoom

    cdef np.ndarray[ndim = 3, dtype = np.int64_t] rx = np.zeros((nx, ny, ang.size)).astype(np.int64)
    cdef np.ndarray[ndim = 3, dtype = np.int64_t] ry = np.zeros((nx, ny, ang.size)).astype(np.int64)
    cdef np.ndarray[ndim = 3, dtype = np.float32_t] wx = np.zeros((nx, ny, ang.size)).astype(np.float32)
    cdef np.ndarray[ndim = 3, dtype = np.float32_t] wy = np.zeros((nx, ny, ang.size)).astype(np.float32)

    cdef np.ndarray[ndim = 3, dtype = np.int64_t] ind = np.zeros((nx, ny, ang.size)).astype(np.int64)

    cdef np.ndarray[ndim = 3, dtype = np.float32_t] imr = np.zeros((im.shape[0], im.shape[1], im.shape[2])).\
                                                    astype(np.float32)

    for i in range(ang.size):
        matrot = np.array([[np.cos(ang[i]), -np.sin(ang[i])],
                         [np.sin(ang[i]), np.cos(ang[i])]], dtype=np.float32)
        wx[:, :, i] = x * matrot[0, 0] + y * matrot[1, 0] + cx
        wy[:, :, i] = x * matrot[0, 1] + y * matrot[1, 1] + cy

    nn = np.where(wx < 1)
    wx[nn] = 1.
    nn = np.where(wy < 1)
    wy[nn] = 1.

    nn = np.where(wx > (nx - 1))
    wx[nn] = nx - 1
    nn = np.where(wy > (ny - 1))
    wy[nn] = ny - 1

    rx = wx.astype(np.int64)  # partie entiere
    ry = wy.astype(np.int64)
    wx -= rx  # partie fractionnaire
    wy -= ry

    ind = rx + (ry - 1) * nx
    if(ang.size > 1):
        for i in range(ang.size):
            ind[:, :, i] += i * nx * ny

    imr.flat = \
            (im.flatten()[ind.flatten()] * 
                    (1 - wx.flatten()) + \
                im.flatten()[ind.flatten() + 1] * wx.flatten())\
             *(1 - wy.flatten()) + \
             (im.flatten()[ind.flatten() + nx] * (1 - wx.flatten()) + im.flatten()[ind.flatten() + nx + 1] * wx.flatten()) * wy.flatten()

    return imr


cpdef rotate(np.ndarray[ndim=3, dtype=np.float32_t] im,
            float ang, float cx=-1, float cy=-1, float zoom=1.0):
    """Rotates an image of an angle "ang" (in DEGREES).

    The center of rotation is cx,cy.
    A zoom factor can be applied.

    (cx,cy) can be omitted :one will assume one rotates around the
    center of the image.
    If zoom is not specified, the default value of 1.0 is taken.

    :parameters:
        im: (np.ndarray[ndim=3,dtype=np.float32_t]) : array to rotate

        ang: (float) : rotation angle (in degrees)

        cx: (float) : (optional) rotation center on x axis (default: image center)

        cy: (float) : (optional) rotation center on x axis (default: image center)

        zoom: (float) : (opional) zoom factor (default =1.0)

    """
# TODO test it
    if(zoom == 0):
        zoom = 1.0
    if(zoom == 1.0 and ang == 0.):
        return im

    ang *= np.pi / 180
    cdef int nx = im.shape[1]
    cdef int ny = im.shape[2]
    cdef np.ndarray[ndim = 2, dtype = np.float32_t] matrot

    cdef int i

    if(cx < 0):
        cx = nx / 2 + 1
    if(cy < 0):
        cy = ny / 2 + 1

    cdef np.ndarray x = np.tile(np.arange(nx) - cx + 1, (ny, 1)).T / zoom
    cdef np.ndarray y = np.tile(np.arange(ny) - cy + 1, (nx, 1)) / zoom

    cdef np.ndarray[ndim = 3, dtype = np.int64_t] rx = np.zeros((nx, ny, ang.size)).astype(np.int64)
    cdef np.ndarray[ndim = 3, dtype = np.int64_t] ry = np.zeros((nx, ny, ang.size)).astype(np.int64)
    cdef np.ndarray[ndim = 3, dtype = np.float32_t] wx = np.zeros((nx, ny, ang.size)).astype(np.float32)
    cdef np.ndarray[ndim = 3, dtype = np.float32_t] wy = np.zeros((nx, ny, ang.size)).astype(np.float32)

    cdef np.ndarray[ndim = 3, dtype = np.int32_t] ind = np.zeros((nx, ny, ang.size))

    cdef np.ndarray[ndim = 3, dtype = np.float32_t] imr = np.zeros((im.shape[0], im.shape[1], im.shape[2])).\
                                                    astype(np.float32)

    matrot = np.array([[np.cos(ang), -np.sin(ang)],
                     [np.sin(ang), np.cos(ang)]])

    wx[:, :] = x * matrot[0, 0] + y * matrot[1, 0] + cx
    wy[:, :] = x * matrot[0, 1] + y * matrot[1, 1] + cy

    nn = np.where(wx < 1)
    wx[nn] = 1.
    nn = np.where(wy < 1)
    wy[nn] = 1.

    nn = np.where(wx > (nx - 1))
    wx[nn] = nx - 1
    nn = np.where(wy > (ny - 1))
    wy[nn] = ny - 1

    rx = wx.astype(np.int64)  # partie entiere
    ry = wy.astype(np.int64)
    wx -= rx  # partie fractionnaire
    wy -= ry

    ind = rx + (ry - 1) * nx

    imr = (im[ind] * (1 - wx) + im[ind + 1] * wx) * (1 - wy) + (im[ind + nx] * (1 - wx) + im[ind + nx + 1] * wx) * wy

    return imr





cpdef  indices(int dim1, int dim2=-1):
    """DOCUMENT indices(dim)

    Return a dimxdimx2 array. First plane is the X indices of the pixels in the dimxdim array. Second plane contains the Y indices.
    Inspired by the Python scipy routine of the same name.
    New (June 12 2002), dim can either be:

    * a single number N (e.g. 128) in which case the returned array are square (NxN)
    * a Yorick array size, e.g. [#dimension,N1,N2], in which case the returned array are N1xN2
    * a vector [N1,N2], same result as previous case

    F.Rigaut 2002/04/03

    :parameters:
        dim1: (int) : first dimension

        dim2: (int) : (optional) second dimension
    """


    if (dim2 < 0):
        y = np.tile((np.arange(dim1, dtype=np.float32) + 1), (dim1, 1))
        x = np.copy(y.T)
        return y, x
    else :
        x = np.tile((np.arange(dim1, np.float32) + 1), (dim2, 1))
        y = np.tile((np.arange(dim2, np.float32) + 1), (dim1, 1)).T
        return y, x



cpdef makegaussian(int size, float fwhm, int xc=-1, int yc=-1, int norm=0):
    """makegaussian(size,fwhm,xc,yc)
    Returns a centered gaussian of specified size and fwhm.
    norm returns normalized 2d gaussian

    :parameters:
        size: (int) :

        fwhm: (float) :

        xc: (int) : (optional) center position on x axis

        yc: (int) : (optional) center position on y axis

        norm: (int) : (optional) normalization
    """
    cdef np.ndarray tmp
    tmp = np.exp(-(mkP.dist(size, xc, yc) / (fwhm / 1.66)) ** 2.)
    if (norm > 0):
        tmp = tmp / (fwhm ** 2.*1.140075)
    return tmp

def get_classAttributes(Param_class):
    """ get_classAttributes(Param_class)
    Return all the attribute names of the given Param_class

    :param Param_class: shesha parameters class
    :return: list of strings (attributes names)

    :example:
        import shesha as ao
        get_classAttributes(ao.Param_wfs)
    """
    d = inspect.getmembers(Param_class)
    return [ i[0] for i in d if inspect.isgetsetdescriptor(i[1])]
