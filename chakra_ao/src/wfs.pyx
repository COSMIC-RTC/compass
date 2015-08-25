
cdef class Param_wfs:

    def set_type(self, str t):
        """Set the attribute  to t
        t -- str"""
        self.type_wfs=t

    def set_nxsub(self, long n):
        """Set the attribute nxsub to n
        n -- long"""
        self.nxsub=n

    def set_npix(self,long  n):
        """Set the attribute npix to n
        n -- long"""
        self.npix=n

    def set_pixsize(self,float p):
        """Set the attribute pixsize to p
        p -- float"""
        self.pixsize=p

    def set_Lambda(self,float L):
        """Set the attribute Lambda to L
        L -- float"""
        self.Lambda=L

    def set_optthroughput(self,float o):
        """Set the attribute optthroughput to o
        o -- float"""
        self.optthroughput=o

    def set_fracsub(self,float f):
        """Set the attribute fracsub to f
        f -- float"""
        self.fracsub=f

    def set_openloop(self,long o):
        """Set the attribute openloop to o
        o -- long"""
        self.openloop=o

    def set_fssize(self,float f):
        """Set the attribute fssize to f
        f -- float"""
        self.fssize=f

    def set_fstop(self,str f):
        """Set the attribute fstop to f
        f -- str"""
        self.fstop=f

    def set_atmos_seen(self, int i):
        """Set attribute atmos_seen to i
        i -- int 
        """
        self.atmos_seen=i

    def set_xpos(self,float x):
        """Set the attribute xpos to x
        x -- float"""
        self.xpos=x

    def set_ypos(self,float y):
        """Set the attribute ypos to y
        y -- float"""
        self.ypos=y

    def set_gsalt(self,float g):
        """Set the attribute gsalt to g
        g -- float"""
        self.gsalt=g

    def set_gsmag(self,float g):
        """Set the attribute gsmag to g
        g -- float"""
        self.gsmag=g

    def set_zerop(self,float z):
        """Set the attribute zerop to 
        z -- float"""
        self.zerop=z

    def set_noise(self,float n):
        """Set the attribute noise to n
        n -- float"""
        self.noise=n

    def set_kernel(self,float k):
        """Set the attribute kernel to k
        k -- float"""
        self.kernel=k

    def set_laserpower(self,float l):
        """Set the attribute laserpower to l
        l -- float"""
        self.laserpower=l

    def set_lltx(self,float l):
        """Set the attribute lltx to l
        l -- float"""
        self.lltx=l

    def set_llty(self,float l):
        """Set the attribute llty to l
        l -- float"""
        self.llty=l

    def set_proftype(self,str p):
        """Set the attribute proftype to p
        p -- str"""
        self.proftype=p

    def set_beamsize(self,float b):
        """Set the attribute beamsize to b
        b -- float"""
        self.beamsize=b

    def set_pyr_ampl(self,float p):
        """Set the attribute pyr_ampl to p
        p -- float"""
        self.pyr_ampl=p

    def set_pyr_npts(self,long p):
        """Set the attribute pyr_npts to p
        p -- long"""
        self.pyr_npts=p

    def set_pyr_loc(self,str p):
        """Set the attribute pyr_loc to p
        p -- str"""
        self.pyr_loc=p

    def set_pyrtype(self,str p):
        """Set the attribute pyrtype to p
        p -- str"""
        self.pyrtype=p

    def set_dms_seen(self, np.ndarray[ndim=1,dtype=np.int32_t] dms_seen):
        """Set the attribute dms_seen to dms_seen
        dms_seen -- numpy array
        """
        self.dms_seen=dms_seen

    def set_lgsreturnperwatt(self, float lpw):
        """Set attribute lgsreturnperwatt to lpw
        lpw -- float """
        self.lgsreturnperwatt=lpw


    def set_altna(self,np.ndarray[ndim=1,dtype=np.float32_t] a):
        """Set attribute _altna to a
        a -- np.ndarray[ndim=1,dtype=np.float32]"""
        self._altna=a


    def set_profna(self,np.ndarray[ndim=1,dtype=np.float32_t] p):
        """Set attribute _profna to p
        p -- np.ndarray[ndim=1,dtype=np.float32]"""
        self._profna=p


    cpdef prep_lgs_prof(self,int nsensors, Param_tel p_tel, 
                        np.ndarray[dtype=np.float32_t] prof,
                        np.ndarray[dtype=np.float32_t] h,
                        float beam, Sensors sensors,
                        bytes center=<bytes>"", int imat=0):
        """The function returns an image array(double,n,n) of a laser beacon elongated by perpective
 effect. It is obtaind by convolution of a gaussian of width "lgsWidth" arcseconds, with the
 line of the sodium profile "prof". The altitude of the profile is the array "h".
 prof     -- np.ndarray[dtype=np.float32] : Na profile intensity, in arbitrary units
 h        -- np.ndarray[dtype=np.float32] : altitude, in meters. h MUST be an array with EQUALLY spaced elements.
 beam     -- float : size in arcsec of the laser beam
 center   -- string: either "image" or "fourier" depending on where the centre should be.
 
 Computation of LGS spot from the sodium profile:
 Everything is done here in 1D, because the Na profile is the result of the convolution of a function
 P(x,y) = profile(x) . dirac(y)
 by a gaussian function, for which variables x and y can be split :
 exp(-(x^2+y^2)/2.s^2)  =  exp(-x^2/2.s^2) * exp(-y^2/2.s^2)
 The convolution is (symbol $ denotes integral)
 C(X,Y) = $$ exp(-x^2/2.s^2) * exp(-y^2/2.s^2) * profile(x-X) * dirac(y-Y)  dx  dy
 First one performs the integration along y
 C(X,Y) = exp(-Y^2/2.s^2)  $ exp(-x^2/2.s^2) * profile(x-X)  dx
 which shows that the profile can be computed by
 - convolving the 1-D profile
 - multiplying it in the 2nd dimension by a gaussian function
 
 If one has to undersample the inital profile, then some structures may be "lost". In this case,
 it's better to try to "save" those structures by re-sampling the integral of the profile, and
 then derivating it afterwards.
 Now, if the initial profile is a coarse one, and that one has to oversample it, then a
 simple re-sampling of the profile is adequate."""

        self._prof1d=prof
        self._profcum=np.zeros(prof.size+1,dtype=np.float32)
        self._profcum[1:]=prof.cumsum()
        cdef float subapdiam=p_tel.diam/self.nxsub  #diam of subap
        cdef np.ndarray[dtype=np.float32_t] xsubs, ysubs,dOffAxis,g
        if(self.nxsub>1):
            xsubs=np.linspace( (subapdiam-p_tel.diam)/2,(p_tel.diam-subapdiam)/2,
                                self.nxsub).astype(np.float32)
        else:
            xsubs=np.zeros(1,dtype=np.float32)
        ysubs=xsubs.copy().astype(np.float32)


        cdef int nP=prof.shape[0]  #number of points of the profile
        cdef float hG=np.sum(h*prof)/np.sum(prof) #center of gravity of the profile
        cdef np.ndarray[dtype=np.float32_t] x=np.arange(self._Ntot).astype(np.float32)\
                                                -self._Ntot/2
        #x expressed in pixels. (0,0) is in the fourier-center

        x=x*self._qpixsize #x expressed in arcseconds
        cdef float dx=x[1]-x[0]
        cdef float dh=h[1]-h[0]

        if(self.nxsub>1):
            dOffAxis=np.sqrt((xsubs[self._validsubsx]-self.lltx)**2+
                             (ysubs[self._validsubsy]-self.llty)**2)
        else:
            dOffAxis=np.sqrt((xsubs-self.lltx)**2+(ysubs-self.llty)**2)

        if(imat>0):
            dOffAxis*=0.

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
                g=np.exp(-(x+self._qpixsize/2)**2/(2*w**2.))
            else:
                g=np.exp(-x**2/(2*w**2.))

       
        self._ftbeam=np.fft.fft(g,axis=0).astype(np.complex64)
        self._beam=g
        #convolved profile in 1D.

        if(xsubs.size>1):
            azimuth=np.arctan2(ysubs[self._validsubsy]-self.llty,
                                xsubs[self._validsubsx]-self.lltx)
        else:
            azimuth=np.arctan2(ysubs-self.llty,
                                xsubs-self.lltx)

        self._azimuth=azimuth

        cdef carma_context *context=carma_context.instance()
        context.set_activeDevice(sensors.sensors.device,1)

        cdef sutra_lgs *lgs=sensors.sensors.d_wfs[nsensors].d_gs.d_lgs

        lgs.lgs_init( self._prof1d.size,hG,h[0],dh,self._qpixsize,
            <float*>dOffAxis.data, <float*>self._prof1d.data,
            <float*>self._profcum.data, <float*>self._beam.data,
            <cuFloatComplex*>self._ftbeam.data, <float*>self._azimuth.data)

        lgs.lgs_update(context.get_device(sensors.sensors.device))
        lgs.lgs_makespot(context.get_device(sensors.sensors.device),0)



    cdef make_lgs_prof1d(self,int nsensors, Param_tel p_tel,
            np.ndarray[dtype=np.float32_t] prof, np.ndarray[dtype=np.float32_t] h,
            float beam, bytes center=<bytes>""):
        """same as prep_lgs_prof but cpu only. original routine from rico"""

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


        cdef int nP=prof.shape[0]
        cdef float hG=np.sum(h*prof)/np.sum(prof)
        cdef np.ndarray[dtype=np.float32_t] x=np.arange(self._Ntot).astype(np.float32)-self._Ntot/2
        # x expressed in pixels. (0,0) is in the fourier-center.
        x=x*self._qpixsize #x expressed in arcseconds
        cdef float dx=x[1]-x[0]
        cdef float dh=h[1]-h[0]

        if(self.nxsub>1):
            dOffAxis=np.sqrt((xsubs[self._validsubsy]-self.lltx)**2+
                             (ysubs[self._validsubsx]-self.llty)**2)
        else:
            dOffAxis=np.sqrt((xsubs-self.lltx)**2+(ysubs-self.llty)**2)

        cdef np.ndarray[ndim=2,dtype=np.float32_t] profi,
        cdef np.ndarray[ndim=1,dtype=np.float32_t] zhc, avg_zhc,avg_x
        profi=np.zeros((self._Ntot,self._nvalid),dtype=np.float32)

        cdef np.ndarray[ndim=1,dtype=np.int32_t] subsdone, dif2do
        cdef np.ndarray[ndim=1,dtype=np.int64_t] inds
        subsdone=np.ones(self._nvalid,dtype=np.int32)
        dif2do =np.zeros(self._nvalid,dtype=np.int32)
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
            print "MAX_IM",max_im
            im=(im.T/max_im).T
        else:
            im=rotate(im,azimuth*180/np.pi,xcent,ycent)
            im=im/np.max(im)

        self._lgskern=im.T


cdef type_present( liste,int pyr, int roof, int sh, int geo):
    """for each type of wfs,
    return 1 if the wfs type is present (0 else)
    """

    cdef int i,l

    cdef bytes wfs_type

    l=len(liste)
    for i in range(l):
        wfs_type=liste[i]
        if(wfs_type=="pyr"):
            pyr=1
        elif (wfs_type=="roof"):
            roof=1
        elif(wfs_type=="sh"):
            sh=1
        elif(wfs_type=="geo"):
            geo=1


cdef wheremax(liste):
    """return the index of the maximum value of the list

    liste -- list of value
    """
    cdef int i,j
    cdef int l=len(liste)
    cdef float m=0
    cdef float val
    j=-1
    for i in range(l):
        val=liste[i]
        if(m<val):
            j=i
            m=val
    return j


def wfs_init( wfs, Param_atmos p_atmos, Param_tel p_tel, Param_geom p_geom,
            Param_target p_target, Param_loop p_loop,int comm_size, int rank, dm=None):
    """
    wfs         : list of Param_wfs
    p_atmos     : Param_atmos
    p_tel       : Param_tel
    p_geom      : Param_geom
    p_target    : Param_target
    p_loop      : Param_loop
    comm_size   : int communicator size
    rank        : int process rank
    dm          : list of Param_dm
    """
    cdef int nsensors=len(wfs)
    cdef int i
    cdef Param_wfs o
    cdef int any_pyr=0
    cdef int any_roof=0
    cdef int any_sh=0
    cdef int any_geo=0

    
    liste=[o.type_wfs for o in wfs]
    type_present(liste,any_pyr,any_roof, any_sh,any_geo)

    #dm = None
    if(wfs[0].dms_seen is None):
        for i in range(nsensors):
            if(not wfs[i].openloop and dm is not None):
                wfs[i].set_dms_seen(np.arange(len(dm),dtype=np.int32)+1)


    cdef int indmax
    # first get the wfs with max # of subaps
    # we'll derive the geometry from the requirements in terms of sampling
    if(any_sh):
        indmax=wheremax([o.nxsub for o in wfs if o.type_wfs=="sh"])
    else:
        indmax=wheremax([o.nxsub for o in wfs])

    #init geometry
    cdef float pixsize0=wfs[0].pixsize
    if(any_geo==0):
        init_wfs_geom(wfs[indmax], wfs[0], indmax, p_atmos, p_tel, p_geom, p_target,
                      p_loop, init=1,comm_size=comm_size, verbose=rank)
    else:
        init_wfs_geom(wfs[indmax], wfs[0], indmax, p_atmos, p_tel, p_geom, p_target,
                      p_loop, init=0,comm_size=comm_size, verbose=rank)

    ##do the same for other wfs
    for i in range(nsensors):
        if(i!=indmax):
            init_wfs_geom(wfs[i], wfs[i], i, p_atmos, p_tel, p_geom, p_target,
                        p_loop, init=0,comm_size=comm_size, verbose=rank)
    # create sensor object on gpu
    # and init sensor gs object on gpu

    #arrays needed to call Sensors constructor
    t_wfs  = [o.type_wfs  for o in wfs]
    #cdef np.ndarray t_wfs  = np.array([o.type_wfs  for o in wfs],dtype=np.str)
    cdef np.ndarray nxsub  = np.array([o.nxsub     for o in wfs],dtype=np.int64)
    cdef np.ndarray nvalid = np.array([o._nvalid   for o in wfs],dtype=np.int64)
    cdef np.ndarray nphase = np.array([o._pdiam    for o in wfs],dtype=np.int64)
    cdef np.ndarray pdiam  = np.array([o._subapd   for o in wfs],dtype=np.float32)
    cdef np.ndarray npix   = np.array([o.npix      for o in wfs],dtype=np.int64)
    cdef np.ndarray nrebin = np.array([o._nrebin   for o in wfs],dtype=np.int64)
    cdef np.ndarray nfft   = np.array([o._Nfft     for o in wfs],dtype=np.int64)
    cdef np.ndarray ntota  = np.array([o._Ntot     for o in wfs],dtype=np.int64)
    cdef np.ndarray nphot  = np.array([o._nphotons for o in wfs],dtype=np.float32)
    cdef np.ndarray lgs    = np.array([o.gsalt>0   for o in wfs],dtype=np.int32)

    #arrays needed to call sensors_initgs 
    cdef np.ndarray xpos   = np.array([o.xpos   for o in wfs], dtype=np.float32)
    cdef np.ndarray ypos   = np.array([o.ypos   for o in wfs], dtype=np.float32)
    cdef np.ndarray Lambda = np.array([o.Lambda for o in wfs], dtype=np.float32)
    cdef np.ndarray mag
    cdef float zerop= wfs[0].zerop
    cdef np.ndarray size   = np.zeros(nsensors,dtype=np.int64)+p_geom._n
    cdef np.ndarray noise
    cdef np.ndarray seed = np.array([],dtype=np.int64)
    cdef np.ndarray npup = (np.zeros((nsensors))+p_geom._n).astype(np.int64)

    cdef np.ndarray tmp


    if(wfs[0].type_wfs=="sh"):
        g_wfs= Sensors(nsensors,t_wfs,npup,nxsub,nvalid,nphase,pdiam,npix,nrebin,
                nfft,ntota,nphot,lgs,comm_size=comm_size, rank=rank)

        mag=np.array([o.gsmag    for o in wfs], dtype=np.float32)
        noise=np.array([o.noise    for o in wfs], dtype=np.float32)
        g_wfs.sensors_initgs(xpos,ypos,Lambda,mag,zerop,size,noise,seed)

    elif(wfs[0].type_wfs=="pyr" or wfs[0].type_wfs=="roof"):
        npup=np.array([wfs[0].pyr_npts])
        g_wfs= Sensors(nsensors,t_wfs,npup,nxsub,nvalid,nphase,pdiam,npix,nrebin,
                nfft,ntota,nphot,lgs,comm_size=comm_size, rank=rank)

        mag=np.array([o.gsmag    for o in wfs], dtype=np.float32)
        noise=np.array([o.noise    for o in wfs], dtype=np.float32)
        g_wfs.sensors_initgs(xpos,ypos,Lambda,mag,zerop,size,noise,seed)
    
        
    elif(wfs[0].type_wfs=="geo"):
        npup=np.array([wfs[0].p_geom._n])
        g_wfs= Sensors(nsensors, wfs[0].type_wfs,npup,nxsub,nvalid,nphase,pdiam,
                       comm_size=comm_size, rank=rank)

        mag=np.zeros(nsensors-1,dtype=np.float32)
        noise=np.zeros(nsensors-1,dtype=np.float32)-1
        g_wfs.sensors_initgs(xpos,ypos,Lambda,mag,zerop,size,noise,seed)

    # fill sensor object with data
    cdef int *ph=NULL
    for i in range(nsensors):
        g_wfs.sensors_initarr(i,wfs[i],p_geom)

    #lgs case
    for i in range(nsensors):
        if(wfs[i].gsalt>0):
            #lgs mode requested
            if(wfs[i].proftype is None and wfs[i].proftype==""):
                wfs[i].set_proftype("Gauss1")

            if(wfs[i].proftype=="Gauss1"):
                profilename = "allProfileNa_withAltitude_1Gaussian.npy"
            elif(wfs[i].proftype=="Gauss2"):
                profilename = "allProfileNa_withAltitude_2Gaussian.npy"
            elif(wfs[i].proftype=="Gauss3"):
                profilename = "allProfileNa_withAltitude_3Gaussian.npy"
            elif(wfs[i].proftype=="Exp"):
                profilename = "allProfileNa_withAltitude.npy"
            else:
                error="Param_wfs proftype unknown: got '"+wfs[i].proftype+"', expect one of: \n''\n'Gauss1'\n'Gauss2'\n'Gauss3'\n'Exp'"
                raise ValueError(error)
            prof=np.load(chakra_ao_savepath+profilename)
            wfs[i].set_altna(prof[0,:].astype(np.float32))
            wfs[i].set_profna(np.mean(prof[1:,:],axis=0).astype(np.float32))
            # init sensor gs object with necessary data
            wfs[i].prep_lgs_prof(i,p_tel,wfs[i]._profna,wfs[i]._altna,
                                  wfs[i].beamsize,g_wfs)

    return g_wfs

cdef init_wfs_geom(Param_wfs wfs, Param_wfs wfs0, int n, Param_atmos atmos,
                Param_tel tel, Param_geom geom, Param_target p_target,
                Param_loop loop, int init=0, int comm_size=1, int verbose=0):
    """TODO doc
    """
    if(verbose==0):print "*-----------------------"
    if(verbose==0):print "Doing inits on WFS", n

    cdef long pdiam=0

    if(init==0):
        if(wfs.type_wfs=="sh"):
            pdiam=geom.pupdiam/wfs.nxsub
            if(geom.pupdiam%wfs.nxsub>0):
                    pdiam+=1
        if(wfs.type_wfs=="geo"):
            if(geom.pupdiam==0):
                pdiam=20
            else:
                pdiam=geom.pupdiam/wfs.nxsub
                if(geom.pupdiam%wfs.nxsub>0):
                    pdiam+=1
    else:
        pdiam=-1

    
    cdef int Nfft=0 
    cdef int Ntot=0
    cdef int nrebin=0
    cdef float pixsize=0
    cdef float qpixsize=0
    cdef np.ndarray pyr_focmask, pupvalid, pupreb,istart,jstart#, validsubs
    cdef np.ndarray phase_shift,pshift,cx,cy, phasemap,tmp,sincar, fluxPerSub,halfxy
    cdef np.ndarray validx,validy, isvalid
    cdef np.ndarray binmap,binindices

    cdef np.ndarray x=np.empty((1,1),dtype=np.float32)
    cdef np.ndarray y=np.empty((1,1),dtype=np.float32)

    cdef int i,j,indi#, indj
    cdef long n1,n2,npup
    cdef float coef1,coef2


    #TODO define psize properly
    #defined in yoga_rtc.i /and yoga_dm.i as 
    #   (1): y_geom.pupdiam
    #   (2): y_tel.diam/y_geom.pupdiam
    cdef int psize=0#int(tel.diam/geom.pupdiam)
    init_wfs_size(wfs,wfs0,n,atmos,tel,psize,&pdiam,&Nfft,&Ntot,&nrebin,&pixsize,&qpixsize,verbose)

    if(wfs.type_wfs!="geo"):
        wfs.pixsize   = pixsize
        wfs._Nfft     = Nfft
        wfs._Ntot     = Ntot
        wfs._nrebin   = nrebin
        wfs._qpixsize = qpixsize

    wfs._subapd   = tel.diam/wfs.nxsub
    wfs._pdiam    = pdiam


    if(wfs.type_wfs=="pyr" or wfs.type_wfs=="roof"):
        wfs.npix = pdiam

    if(init==1 or ( wfs.type_wfs=="geo" and n==1)):
        #this is the wfs with largest # of subaps
        #the overall geometry is deduced from it
        geom.geom_init(tel,pdiam*wfs.nxsub,p_target.apod)
    if(wfs.type_wfs=="pyr" or wfs.type_wfs=="roof"):
        padding=2
        npup  =  wfs._Ntot;
        n1    = geom.ssize / 2 - geom.pupdiam / 2 - padding * wfs.npix;
        n2    = n1 + geom.pupdiam + 2 * padding * wfs.npix;

        geom._mpupil = geom._ipupil[n1:n2,n1:n2]
        geom._n1     = n1
        geom._n2     = n2
        geom._n      = npup

        #pup   = pup(ii,ii);
        #phase = phase(ii,ii);
        mod_ampl_pixels = wfs.pyr_ampl / wfs._qpixsize # modulation in pixels
        fsradius_pixels = long(wfs.fssize / wfs._qpixsize / 2.)

        if (wfs.fstop=="round"):
            focmask = mkP.dist(npup,xc=npup/2.+0.5,yc=npup/2.+0.5)<(fsradius_pixels);
            fstop_area = np.pi * (wfs.fssize/2.)**2.;
        elif (wfs.fstop=="square"): 
              x,y = indices(npup) 
              x-=(npup+1.)/2. 
              y-=(npup+1.)/2. 
              focmask = ( np.abs(x) <= (fsradius_pixels) ) *     \
                ( np.abs(y) <= (fsradius_pixels) )
              fstop_area = wfs.fssize**2.
        else:
             msg="wfs "+str(n)+". fstop must be round or square"
             raise ValueError(msg)
    
        pyr_focmask = np.roll(focmask,focmask.shape[0]/2,axis=0)
        pyr_focmask = np.roll(pyr_focmask,focmask.shape[1]/2,axis=1)
        wfs._submask = pyr_focmask
        
        pup = geom._spupil
        pupreb = bin2d(pup*1.,wfs.npix)/wfs.npix**2.
        wsubok = np.where(pupreb>=wfs.fracsub)
        pupvalid = pupreb * 0.
        pupvalid[wsubok] = 1
        wfs._isvalid    = pupvalid.astype(np.int32)

        
        pup = geom._mpupil

        pupreb = bin2d(pup*1.,wfs.npix)/wfs.npix**2.
        wsubok = np.where(pupreb>=wfs.fracsub)
        pupvalid = pupreb * 0.
        pupvalid[wsubok] = 1
        wfs._nvalid     = wsubok[0].size
        validx=np.where(pupvalid)[1].astype(np.int32)
        validy=np.where(pupvalid)[0].astype(np.int32)
        wfs._validsubsx=validx
        wfs._validsubsy=validy

            
        istart =(np.linspace(0.5, geom.pupdiam + 2 * padding * wfs.npix + 0.5,wfs.nxsub+2*padding)+1).astype(np.int32)[:-1]
        jstart=np.copy(istart) 
        wfs._istart= istart.astype(np.int32)
        wfs._jstart= jstart.astype(np.int32)

        x,y = indices(npup)
        x-=(npup+1)/2.
        y-=(npup+1)/2.
        phase_shift = np.roll( np.exp(1j*2*np.pi*(0.5*(x+y))/npup), x.shape[0]/2, axis=0 )
        phase_shift = np.roll( phase_shift, x.shape[1]/2, axis=1 )
        wfs._halfxy = phase_shift

        x,y= indices(Nfft)
        x-=(Nfft+1)/2.
        y-=(Nfft+1)/2.
        if(wfs.nxsub*nrebin%2==1):
            coef1=0.
        else:
            coef1=-0.5
        if(wfs.nxsub%2==1 and wfs.npix*nrebin%2==1):
            coef2=1
        else:
            coef2=0.5

        pshift = np.exp(1j*2*np.pi*(coef1/Nfft+coef2*nrebin/wfs.npix/Nfft)*(x+y))
        wfs._pyr_offsets = pshift

        if(wfs.pyrtype=="Pyramid"):
            if(wfs.pyr_pos == None):
                cx=np.round(mod_ampl_pixels*np.sin((np.arange(wfs.pyr_npts)+1)*2.*np.pi/wfs.pyr_npts))
                cy=np.round(mod_ampl_pixels*np.cos((np.arange(wfs.pyr_npts)+1)*2.*np.pi/wfs.pyr_npts))
                mod_npts = wfs.pyr_npts
            else:
                if(verbose==0):print "Using user-defined positions for the pyramid modulation"
                cx=np.round(wfs.pyr_pos[:,0]/qpixsize)
                cy=np.round(wfs.pyr_pos[:,1]/qpixsize)
                mod_npts=cx.shape[0]
        elif(wfs.pyrtype=="RoofPrism"):
            cx = np.round(2.*mod_ampl_pixels*((np.arange(wfs.pyr_npts)+1)-(wfs.pyr_npts+1)/2.)/wfs.pyr_npts)
            cy = cx
            mod_npts = wfs.pyr_npts
        else:
            if(wfs.pyr_pos==None):
                cx = np.round(mod_ampl_pixels*np.sin((np.arange(wfs.pyr_npts)+1)*2.*np.pi/wfs.pyr_npts))
                cy = np.round(mod_ampl_pixels*np.cos((np.arange(wfs.pyr_npts)+1)*2.*np.pi/wfs.pyr_npts))
                mod_npts = wfs.pyr_npts
            else:
                if(verbose==0):print "Using user-defined positions for the pyramid modulation"
                cx=np.round(wfs.pyr_pos[:,0]/qpixsize)
                cy=np.round(wfs.pyr_pos[:,1]/qpixsize)
                mod_npts=cx.shape[0]

        wfs._pyr_cx=cx.astype(np.int32)
        wfs._pyr_cy=cy.astype(np.int32)

        wfs._nphotons = wfs.zerop*2.51189**(-wfs.gsmag)*loop.ittime*wfs.optthroughput
    
        # spatial filtering by the pixel extent:
        # *2/2 intended. min should be 0.40 = sinc(0.5)^2.
        x=x/(Nfft-1)*2/2
        y=y/(Nfft-1)*2/2
        sincar = np.roll(np.sinc(np.pi*x)*np.sinc(np.pi*y),x.shape[0],axis=0)
        sincar = np.roll(np.pi*x*np.pi*y,x.shape[1],axis=1)
        wfs._sincar = sincar.astype(np.float32)
        #must be set even if unused
        wfs._hrmap=np.array([0],dtype=np.int32)

        # this defines how we cut the phase into subaps
        phasemap=np.zeros((pdiam,pdiam,wfs._nvalid),dtype=np.int32)
        x,y=indices(geom._n) #we need c-like indice
        x-=1
        y-=1
        tmp=x+y*geom._n
        for i in range(wfs._nvalid):
            indi=istart[wfs._validsubsx[i]]+1 #+2-1 (yorick->python
            indj=jstart[wfs._validsubsy[i]]+1
            phasemap[:,:,i]=tmp[indi:indi+pdiam, indj:indj+pdiam]

        
        wfs._phasemap=phasemap

    if(wfs.type_wfs == "sh" or wfs.type_wfs == "geo"):
        # this is the i,j index of lower left pixel of subap
        istart =((np.linspace(0.5, geom.pupdiam + 0.5  ,wfs.nxsub+1)+1)[:-1]).astype(np.int64)

        
        jstart=np.copy(istart)
        wfs._istart = istart.astype(np.int32)
        wfs._jstart = jstart.astype(np.int32)

        # sorting out valid subaps
        fluxPerSub=np.zeros((wfs.nxsub,wfs.nxsub),dtype=np.float32)
        
        for i in range(wfs.nxsub):
            indi=istart[i]+1 #+2-1 (yorick->python
            for j in range(wfs.nxsub):
                indj=jstart[j]+1 #+2-1 (yorick->python
                fluxPerSub[i,j] = np.sum(geom._mpupil[indi:indi+pdiam,indj:indj+pdiam])

        fluxPerSub = fluxPerSub/pdiam**2.

        pupvalid = (fluxPerSub >= wfs.fracsub)*1
        wfs._isvalid= pupvalid.astype(np.int32)
        wfs._nvalid=np.sum(pupvalid)
        wfs._fluxPerSub =fluxPerSub
        if(verbose==0):
            pupvalid2=pupvalid[:,pupvalid.shape[0]/2:]
        else:
            pupvalid2=pupvalid[:,:pupvalid.shape[0]/2]
        validx=np.where(pupvalid)[1].astype(np.int32)
        validy=np.where(pupvalid)[0].astype(np.int32)
        wfs._validsubsx=validx
        wfs._validsubsy=validy

        # this defines how we cut the phase into subaps
        phasemap=np.zeros((pdiam*pdiam,wfs._nvalid),dtype=np.int32)

        x,y=indices(geom._n)
        x-=1
        y-=1
        tmp=x+y*geom._n

        n=wfs._nvalid
        #for i in range(wfs._nvalid):
        for i in range(n):
            indi=istart[wfs._validsubsx[i]]+1 #+2-1 (yorick->python
            indj=jstart[wfs._validsubsy[i]]+1
            phasemap[:,i]=tmp[indi:indi+pdiam, indj:indj+pdiam].flatten("C")
        wfs._phasemap=phasemap

        #this is a phase shift of 1/2 pix in x and y
        if(wfs.type_wfs=="sh"):
            halfxy=np.linspace(0,2*np.pi,wfs._Nfft+1)[0:wfs._pdiam] / 2.
            halfxy=np.tile(halfxy,(wfs._pdiam,1))
            halfxy+=halfxy.T
            #wfs._halfxy = <float*>(halfxy*0.).data #dont work: half*0 is temp python obj

            if(wfs.npix % 2 == 1 and wfs._nrebin %2 == 1):
                #wfs._halfxy = <float*>(halfxy*0.)
                halfxy=np.zeros((wfs._pdiam,wfs._pdiam),dtype=np.float32)
                wfs._halfxy=halfxy.astype(np.float32) 
            else:
                wfs._halfxy = halfxy.astype(np.float32)

        else:
            halfxy=np.linspace(0,2*np.pi,wfs._pdiam+1)[1:wfs._pdiam]/2
            halfxy=np.tile(halfxy,(wfs._pdiam,1))
            halfxy=halfxy*0.
            wfs._halfxy = halfxy.astype(np.float32)

        #must be set even if unused
        wfs._submask=np.array([0],dtype=np.float32)

    if (wfs.type_wfs == "sh"):
        #this defines how we create a larger fov if required
        if(wfs._Ntot!=wfs._Nfft):
            indi=long((wfs._Ntot-wfs._Nfft)/2.) #+1 -1 (yorick>python)
            indj=long(indi+wfs._Nfft)
            x,y=indices(wfs._Nfft)
            #hrpix
            tmp=np.zeros((wfs._Ntot,wfs._Ntot))
            tmp[indi:indj,indi:indj]=np.roll( x+(y-1)*wfs._Nfft, wfs._Nfft/2,axis=0)
            tmp[indi:indj,indi:indj]=np.roll( tmp[indi:indj,indi:indj], wfs._Nfft/2,axis=1)
            #hrmap=roll(hrpix)
            tmp=np.roll(tmp,wfs._Ntot/2,axis=0)
            tmp=np.roll(tmp,wfs._Ntot/2,axis=1)

            tmp=np.where(tmp.flatten())[0]

            wfs._hrmap=np.copy(tmp.astype(np.int32))
            #must be set even if unused
            wfs._sincar=np.array([0],dtype=np.float32)

        else:
            tmp=np.zeros((1))
            wfs._hrmap=np.copy(tmp.astype(np.int32))
            #must be set even if unused
            wfs._sincar=np.array([0],dtype=np.float32)

        if(wfs._nrebin*wfs.npix %2 != wfs._Ntot%2):
            indi=long((wfs._Ntot-wfs._nrebin*wfs.npix)/2.)+1 #+2-1 yorick>python
        else:
            indi=long((wfs._Ntot-wfs._nrebin*wfs.npix)/2.)+0
        indj=long(indi+wfs._nrebin*wfs.npix-1)

        x,y=indices(wfs._nrebin*wfs.npix)
        x=((x-1)/wfs._nrebin).astype(np.int64)
        y=((y-1)/wfs._nrebin).astype(np.int64)

        #binindices
        binindices=np.zeros((wfs._Ntot,wfs._Ntot))
        binindices[indi:indj+1,indi:indj+1]=x+y*wfs.npix+1


        binmap=np.zeros((wfs._nrebin*wfs._nrebin,wfs.npix*wfs.npix))

        x,y=indices(wfs._Ntot)
        x-=1
        y-=1
        tmp=x+y*wfs._Ntot

        if(wfs.gsalt<=0):
            binindices=np.roll(binindices,binindices.shape[0]/2,axis=0)
            binindices=np.roll(binindices,binindices.shape[1]/2,axis=1)

        for i in range(wfs.npix*wfs.npix):
            binmap[:,i]=tmp[np.where(binindices==i+1)]
        #binmap=np.reshape(binmap.flatten("F"),(binmap.shape[0],binmap.shape[1]))
        wfs._binmap=np.copy(binmap.astype(np.int32))


        dr0=tel.diam/atmos.r0*(0.5/wfs.Lambda)**1.2 /np.cos(geom.zenithangle*dtor)**0.6
        fwhmseeing = wfs.Lambda/(tel.diam/np.sqrt(wfs.nxsub**2.+(dr0/1.5)**2.))/4.848
        kernelfwhm = np.sqrt(fwhmseeing**2.+wfs.kernel**2.)


        tmp=makegaussian(wfs._Ntot,kernelfwhm/wfs._qpixsize,wfs._Ntot/2+1,
                            wfs._Ntot/2+1).astype(np.float32)

        tmp=np.roll(tmp,tmp.shape[0]/2,axis=0)
        tmp=np.roll(tmp,tmp.shape[1]/2,axis=1)

        tmp[0,0]=1.     #this insures that even with fwhm=0, the kernel is a dirac
        tmp=tmp/np.sum(tmp)
        tmp=np.fft.fft2(tmp).astype(np.complex64)/(wfs._Ntot*wfs._Ntot)
        wfs._ftkernel = np.copy(tmp).astype(np.complex64)

        # dealing with photometry
        telSurf  = np.pi/4.*(1-tel.cobs**2.)*tel.diam**2.

        # from the guide star
        if(wfs.gsalt == 0):
            if(wfs.zerop == 0):
                wfs.zerop = 1e11
            wfs._nphotons = wfs.zerop*10**(-0.4*wfs.gsmag)*\
                     wfs.optthroughput* \
                    (tel.diam/wfs.nxsub)**2./telSurf* \
                    loop.ittime                        
# include throughput to WFS
# for unobstructed subaperture
# per iteration

        else: # we are dealing with a LGS
            wfs._nphotons = wfs.lgsreturnperwatt* \
                    wfs.laserpower*               \
                    wfs.optthroughput*            \
                    (tel.diam/wfs.nxsub)**2.*1e4* \
                    loop.ittime
# detected by WFS
# ... for given power
# include throughput to WFS
# for unobstructed subaperture
# per iteration

        if(verbose==0):print "nphotons : ", wfs._nphotons


cdef init_wfs_size( Param_wfs wfs, Param_wfs wfs0, int n, Param_atmos atmos,
                Param_tel tel, int psize, long *pdiam, int *Nfft, int *Ntot, int *nrebin, 
                float *pixsize, float *qpixsize,int verbose=0):
    """   init_wfs_size(wfs,pdiam,Nfft,Ntot,nrebin,pixsize,qpixsize)

    computes the quatum pixel sizes and all useful dimensions for a given wfs
    requires 2 externals : y_atmos & y_tel
    y_atmos : a atmos_struct()
    y_tel   : a tel_struct()
    wfs     : wfs_struct as input (input)
    pdiam   : pupil diam for each subap (pixels) (output)
    Nfft    : fft size for a subap (pixels) (output)
    Ntot    : hr image size for a subap (pixels) (output)
    nrebin  : rebin factor for a subap (output)
    pixsize : simulated pixel size for a subap (arcsec) (output)
    qpixsize: quantum pixel size for a subap (arcsec) (output)

    Scheme to determine arrays sizes
    sh :
    k = 6
    p = k * d/r0
    n = long(2*d*v/lambda/RASC)+1
    N = fft_goodsize(k*n/v*lambda/r0*RASC)
    u = k * lambda / r0 * RASC / N
    n = v/u - long(v/u) > 0.5 ? long(v/u)+1 : long(v/u)
    v = n * u
    Nt = v * Npix


    """ 
    cdef float r0 = atmos.r0*(wfs.Lambda*2)**(6./5)
    cdef int k, padding
    cdef long nphase, npix,
    cdef float subapdiam

    cdef np.ndarray w,npix_ok 

    if(r0!=0):
        if(verbose==0):print "r0 for WFS :","%3.2f"%r0," m"
        #seeing = RASC * (wfs.lambda * 1.e-6) / r0
        if(verbose==0):print "seeing for WFS : ","%3.2f"%(RASC * (wfs.Lambda * 1.e-6) / r0),"\""

    if(pdiam[0]<=0):
    # this case is usualy for the wfs with max # of subaps
    # we look for the best compromise between pixsize and fov
        subapdiam = tel.diam / float(wfs.nxsub)  # diam of subap
        k=6
        pdiam[0]=long(k * subapdiam / r0) # number of phase points per subap
        if (pdiam[0] < 16):
            pdiam[0] = 16

        if(wfs.type_wfs=="sh"):
            nrebin[0] = long(2 * subapdiam * wfs.pixsize / (wfs.Lambda*1.e-6) / RASC) + 1
            nrebin[0] = max(2,nrebin[0])
            # first atempt on a rebin factor
    
            # since we clipped pdiam we have to be carreful in nfft computation
            Nfft[0] = fft_goodsize(long(pdiam[0]/ subapdiam * nrebin[0] / wfs.pixsize * RASC * (wfs.Lambda*1.e-6)))
            # size of the support in fourier domain
    
            #qpixsize = k * (wfs.Lambda*1.e-6) / r0 * RASC / Nfft
            qpixsize[0] = (pdiam[0] * (wfs.Lambda*1.e-6) / subapdiam  * RASC) / Nfft[0]

        if(wfs.type_wfs=="pyr" or wfs.type_wfs=="roof"):
            #while (pdiam % wfs.npix != 0) pdiam+=1;
            padding = 2
            nphase  =  pdiam[0] * wfs.nxsub+ 2 * padding * pdiam[0]
            qpixsize[0]   = (pdiam[0] * wfs.nxsub / float(nphase))*wfs.Lambda/tel.diam/4.848
      
            fssize_pixels = long(wfs.fssize / qpixsize[0] / 2.)
            #nrebin  = pdiam / wfs.npix
      
            npix = (wfs.nxsub+2*padding)
            npix_ok = npix*(np.arange(100)+1)
            npix_ok = npix_ok[np.where(npix_ok.flatten()<=(nphase/2.))[0]]

            # now we want the quadrant image dim to be > fssize_pixels:
            w = np.where(npix_ok.flatten()>=fssize_pixels)[0]
            if (w.size==0):
                maxfs = npix_ok[0]*2*psize
                msg="wfs ",n,". ffsize too large (max=",maxfs,")!"
                raise ValueError(msg)
      
            #npix_ok = npix_ok[w[1]]

            nrebin[0]  = npix_ok[w[0]]/npix
            Nfft[0]    = npix_ok[w[0]]
            Ntot[0]    = nphase
            pixsize[0] = qpixsize[0] * nrebin[0]

        #quantum pixel size
    else:
        #this case is for a wfs with fixed # of phase points
        subapdiam = tel.diam / float(wfs.nxsub) # diam of subap
        if (wfs.type_wfs != "geo"):
            Nfft[0] = fft_goodsize(2* pdiam[0]);
            # size of the support in fourier domain
  
            qpixsize[0] = pdiam[0] * (wfs.Lambda*1.e-6) / subapdiam * RASC / Nfft[0];
            # quantum pixel size
    if (wfs.type_wfs == "sh"):
        # actual rebin factor
        if(wfs.pixsize/qpixsize[0] - long(wfs.pixsize/qpixsize[0]) > 0.5):
            nrebin[0] = long(wfs.pixsize/qpixsize[0])+1 
        else:
            nrebin[0] = long(wfs.pixsize/qpixsize[0])

        # actual pixel size
        pixsize[0] = nrebin[0] * qpixsize[0]

        if (pixsize[0] * wfs.npix > qpixsize[0] * Nfft[0]): 
            Ntot[0] = fft_goodsize(long(pixsize[0] * wfs.npix / qpixsize[0]) + 1);
        else:
            Ntot[0] = Nfft[0]

        if (Ntot[0]%2 != Nfft[0]%2):
            Ntot[0]+=1

    if (wfs.type_wfs != "geo"):
        if(verbose==0):print "quantum pixsize : ", "%5.4f"%qpixsize[0],"\""
        if(verbose==0):print "simulated FoV : ","%3.2f"%(Ntot[0] * qpixsize[0]) ,"\" x ","%3.2f"%(Ntot[0] * qpixsize[0]),"\""
        if(verbose==0):print "actual pixsize : ", "%5.4f"%pixsize[0]
        if(verbose==0):print "actual FoV : ","%3.2f"%(pixsize[0] * wfs.npix) ,"\" x ","%3.2f"%(pixsize[0] * wfs.npix),"\""
        if(verbose==0):print "number of phase points : ",pdiam[0]
        if(verbose==0):print "size of fft support : ",Nfft[0]
        if (Ntot[0] > Nfft[0]):
            if(verbose==0):print "size of HR spot support : ",Ntot[0]


    
