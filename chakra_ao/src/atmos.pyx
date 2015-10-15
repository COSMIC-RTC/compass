
#################################################
# P-Class (parametres) Param_atmos
#################################################
cdef class Param_atmos:
    def __cinit__(self):
        self.nscreens=0

    def set_nscreens( self, long n):
        """Set the number of turbulent layers

        :param n: (long) number of screens.
        """
        self.nscreens=n

    def set_r0( self, float r):
        """Set the global r0

        :param r: (float) : global r0
        """
        self.r0=r

    def set_pupixsize(self,float xsize):
        """Set the pupil pixel size

        :param xsize: (float) : pupil pixel size
        """
        self.pupixsize= xsize

    def set_L0( self, l):
        """Set the L0 per layers 

        :param l: (lit of float) : L0 for each layers
        """
        self.L0=np.array(l,dtype=np.float32)

    def set_dim_screens( self, l):
        """Set the size of the phase screens 

        :param l: (lit of float) : phase screens sizes
        """
        self.dim_screens=np.array(l,dtype=np.float32)

    def set_alt( self, l):
        """Set the altitudes of each layer

        :param l: (lit of float) : altitudes
        """
        self.alt=np.array(l,dtype=np.float32)

    def set_winddir( self, l):
        """Set the wind direction for each layer

        :param l: (lit of float) : wind directions
        """
        self.winddir=np.array(l,dtype=np.float32)

    def set_windspeed( self, l):
        """Set the the wind speed for each layer

        :param l: (lit of float) : wind speeds
        """
        self.windspeed=np.array(l,dtype=np.float32)

    def set_frac( self, l):
        """Set the fraction of r0 for each layers

        :param l: (lit of float) : fraction of r0
        """
        self.frac=np.array(l,dtype=np.float32)

    def set_deltax( self, l):
        """Set the translation speed on axis x for each layer

        :param l: (lit of float) : translation speed
        """
        self.deltax=np.array(l,dtype=np.float32)

    def set_deltay( self, l):
        """Set the translation speed on axis y for each layer

        :param l: (lit of float) : translation speed
        """
        self.deltay=np.array(l,dtype=np.float32)

    def set_seeds( self, l):
        """Set the seed for each layer

        :param l: (lit of float) : seed
        """
        self.seeds=np.array(l,dtype=np.float32)



    def atmos_init(self, chakra_context c, Param_tel tel,  Param_geom geom,
                    Param_loop loop, Param_wfs wfs=None, Param_target target=None,
                    int rank=0):
        """Create and initialise an atmos object
        TODO doc

        :parameters:
            c: (chakra_context) : context

            tel: (Param_tel) : telescope settings

            geom: (Param_geom) : geometry settings

            loop: (Param_loop) : loop settings

            wfs: (Param_wfs) : (optional) wfs settings

            target: (Param_target) : (optional) target_settings
        """

        cdef double ittime=loop.ittime
        if(self.r0 is None):
            self.r0=0

        cdef int i

        # ajust layers alt using zenith angle
        self.alt=(self.alt/np.cos(geom.zenithangle*dtor))
        # pixel size in meter
        self.pupixsize=tel.diam/geom.pupdiam


        cdef long max_size
        # compute total fov using targets and wfs gs
        if ((wfs is not None) and (target is not None)):
            max_size = np.max((np.linalg.norm(target.xpos,target.ypos),
                            np.linalg.norm(wfs.xpos,wfs.ypos)))
        elif (target is not None):
                    max_size = np.max(np.linalg.norm(target.xpos,target.ypos))
        elif (wfs is not None):
                max_size = np.max(np.linalg.norm(wfs.xpos,wfs.ypos))
        else:
            max_size = 0


        # compute corresponding meta-pupil diameters at each alt
        cdef np.ndarray patch_diam
        patch_diam = geom._n+2*(max_size*4.84814e-6* self.alt)/ self.pupixsize+4
        self.dim_screens = (patch_diam+patch_diam % 2).astype(np.int64)

        #compute phase screens speed in pixels / iteration
        self.deltax  = geom.pupdiam/tel.diam*self.windspeed*np.cos(dtor*geom.zenithangle)*ittime
        self.deltay  = self.deltax*np.sin(dtor*(self.winddir))
        self.deltax  = self.deltax*np.cos(dtor*(self.winddir))

        if(self.frac.size==1):
            self.frac=np.ones((1),dtype=np.float32)
        else:
            self.frac /= sum(self.frac);
      
        if (self.L0 is None): 
            self.L0 = np.ones(self.nscreens,dtype=np.float32)*(1.e5)# infinite L0
        else: 
            if (self.L0.shape[0] == 1):
                self.L0 = np.ones(self.nscreens,dtype=np.float32)*self.L0[0]
            for i in range(self.nscreens):
                frac_l0=tel.diam/self.L0[i]
                self.L0[i]=geom.pupdiam/frac_l0
        return atmos_create(c,self.nscreens,self.r0,self.L0,self.pupixsize,
            self.dim_screens,self.frac,self.alt,self.windspeed,
            self.winddir,self.deltax,self.deltay,rank)




#################################################
# P-Class atmos
#################################################
cdef class Atmos:
    def __cinit__(self):
        self.context = None
        
    cdef realinit(self,chakra_context ctxt,int nscreens,
                np.ndarray[ndim=1,dtype=np.float32_t] r0,
                np.ndarray[ndim=1,dtype=np.int64_t] size,
                np.ndarray[ndim=1,dtype=np.float32_t] altitude,
                np.ndarray[ndim=1,dtype=np.float32_t] windspeed,
                np.ndarray[ndim=1,dtype=np.float32_t] winddir,
                np.ndarray[ndim=1,dtype=np.float32_t] deltax,
                np.ndarray[ndim=1,dtype=np.float32_t] deltay,
                int device ):
        """Create a sutra_atmos object

        :parameters:
            c: (chakra_context): context

            nscreens: (int): number of turbulent layers

            r0: (np.ndarray[ndim=1,dtype=np.float32_t]): global r0

            size: (np.ndarray[ndim=1,dtype=np.int64_t]): screen size of each layer

            altitude: (np.ndarray[ndim=1,dtype=np.float32_t]): altitude of each layer

            windspeed: (np.ndarray[ndim=1,dtype=np.float32_t]): wind speed of each layer

            winddir: (np.ndarray[ndim=1,dtype=np.float32_t]): wind direction of each layer

            deltax: (np.ndarray[ndim=1,dtype=np.float32_t]): x translation speed

            deltay: (np.ndarray[ndim=1,dtype=np.float32_t]): y translation speed

            device: (int): device index
        """
        
        cdef np.ndarray[ndim=1,dtype=np.int64_t]size2
        size2 = compute_size2(size)

        self.s_a = new sutra_atmos(ctxt.c,nscreens,
                <np.float32_t *>r0.data,
                <long*>size.data,
                <long*>size2.data,
                <np.float32_t*>altitude.data,
                <np.float32_t*>windspeed.data,
                <np.float32_t*>winddir.data,
                <np.float32_t*>deltax.data,
                <np.float32_t*>deltay.data,
                device)
        self.context = ctxt



    def get_screen(self, float alt):
        """Return a numpy array containing the turbulence at a given altitude

        :param alt: (float) :altitude of the screen to get
        """
        cdef carma_obj[float] *screen=self.s_a.d_screens[alt].d_tscreen.d_screen
        self.context.set_activeDevice(screen.getDevice(),1)
        cdef const long *dims
        dims=screen.getDims()
        cdef np.ndarray data_F=np.empty((dims[2],dims[1]),dtype=np.float32)
        cdef np.ndarray data  =np.empty((dims[1],dims[2]),dtype=np.float32)
        screen.device2host(<float*>data_F.data)
        data=np.reshape(data_F.flatten("F"),(dims[1],dims[2]))
        return data

    def disp(self,float alt):
        """Display the screen's phase at a given altitude 
        
        :param alt: (float) : altitude of the screen to display
        """
        cdef carma_obj[float] *c_phase =self.s_a.d_screens[alt].d_tscreen.d_screen
        cdef const long *dims = c_phase.getDims()
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data=np.zeros((dims[2],dims[1]),dtype=np.float32)
        cdef np.ndarray[ndim=2,dtype=np.float32_t] phase=np.ndarray((dims[1],dims[2]),
                                                dtype=np.float32)

        c_phase.device2host(<float*>data.data)
        phase=np.reshape(data.flatten("F"),(dims[1],dims[2]))

        pl.ion()
        pl.clf()
        pl.imshow(phase,cmap="Blues")
        pl.show()


    def add_screen(self,long size, float amplitude, float altitude, 
        float windspeed, float winddir, float deltax, float deltay, int device):
        """Add a screen to the atmos object.
        TODO doc

        :parameters:
            size: (float) :

            amplitude: (float) :

            altitude: (float) :

            windspeed: (float) :

            winddir: (float) :

            deltax: (float) :

            deltay: (float) :

            device: (int) :
        """
        cdef long size2=compute_size2(np.array([size],dtype=np.int64))[0]

        if(self.s_a.d_screens.find(altitude)!=self.s_a.d_screens.end()):
            print "There is already a screen at this altitude"
            print "No screen created"
            return

        cdef sutra_tscreen *screen= new sutra_tscreen(self.s_a.current_context, size,size2,amplitude,altitude,windspeed,winddir,deltax,deltay,device)

        cdef pair[float, sutra_tscreen *] p
        p.first,p.second=altitude,screen
        self.s_a.d_screens.insert(p)
        self.s_a.nscreens+=1
    
    def del_screen(self,float alt):
        """Delete a screen from the atmos object

        :param alt: (float) : altitude of the screen to delete
        """
        if(self.s_a.d_screens.find(alt)==self.s_a.d_screens.end()):
            print "No screen at this altitude"
            print "No screen deleted"
            return
        self.s_a.nscreens-=1
        self.s_a.d_screens.erase(alt)



    def list_alt(self):
        """Display the list of the screens altitude"""

        cdef map[float,sutra_tscreen*].iterator it
        cdef int i=0
        cdef np.ndarray alt=np.zeros(self.s_a.nscreens,dtype=np.float32)
        it=self.s_a.d_screens.begin()

        while it !=self.s_a.d_screens.end():
            alt[i]=deref(it).first
            inc(it)
            i+=1
        print alt


    def move_atmos(self):
        """Generate turbulence"""

        self.s_a.move_atmos()


    def __str__(self):
        cdef map[float,sutra_tscreen*].iterator it = self.s_a.d_screens.begin()
        cdef sutra_tscreen *screen
        cdef int i=1
        info= "Atmos obect:\n"
        info+= "Contains "+str(self.s_a.nscreens)+" turbulent screen(s):\n" 
        info+= "Screen # | alt.(m) | speed (m/s) | dir.(deg) | r0 (pix) | deltax | deltay\n"
        while it !=self.s_a.d_screens.end():
            screen=deref(it).second
            info+= "%8d"%i+ " | "+"%7.4f"%screen.altitude+" | "+"%11.4f"%screen.windspeed+\
                  " | "+"%9.4f"%screen.winddir+ " | "+ "%8.4f"%screen.amplitude**-(6./5.)+\
                  " | "+ "%6.4f"%screen.deltax+ " | "+ "%6.4f"%screen.deltay+"\n"
            i=i+1
            inc(it)

        info+= "--------------------------------------------------------"

        return info

 
cdef atmos_create(chakra_context c, int nscreens,
              float r0,
              np.ndarray[dtype=np.float32_t] L0,
              float pupixsize,
              np.ndarray[ndim=1,dtype=np.int64_t] dim_screens,
              np.ndarray[ndim=1,dtype=np.float32_t] frac,
              np.ndarray[ndim=1,dtype=np.float32_t] alt,
              np.ndarray[ndim=1,dtype=np.float32_t] windspeed,
              np.ndarray[ndim=1,dtype=np.float32_t] winddir,
              np.ndarray[ndim=1,dtype=np.float32_t] deltax,
              np.ndarray[ndim=1,dtype=np.float32_t] deltay,
              int verbose):
    """Create and initialise an atmos object.
    
    TODO doc
    :parameters:
        c: (chakra_context) : context

        nscreens: (float) : number of turbulent layers 

        r0: (float) : global r0

        L0: (np.ndarray[ndim=1, dtype=np.float32_t]) : L0

        pupixsize: (float) :

        dim_screens: (np.ndarray[ndim=1,dtype=np.int64_t]) :

        frac: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        alt: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        windspeed: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        winddir: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        deltax: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        deltay: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        verbose: (int) :
    """

    # get fraction of r0 for corresponding layer
    r0_layers = r0/(frac**(3./5.)*pupixsize)
    # create atmos object on gpu
    atmos_obj = Atmos()

    cdef carma_context *context=carma_context.instance()
    cdef int device = context.get_activeDevice()
    atmos_obj.realinit(chakra_context(),nscreens, r0_layers, dim_screens,alt,
                    windspeed,winddir,deltax,deltay,device)

    cdef int i,j
    cdef np.ndarray[ndim=2,dtype=np.float32_t] A,B,A_F,B_F
    cdef np.ndarray[ndim=1,dtype=np.uint32_t] istx,isty
    cdef sutra_tscreen *tscreen

    cdef str file_A,file_B,file_istx,file_isty


    for i in range(nscreens):
        file_A=os.path.normpath(chakra_ao_savepath+"/turbu/A_"+str(dim_screens[i])+"_L0_"+str(int(L0[i]))+".npy")
        file_B=os.path.normpath(chakra_ao_savepath+"turbu/B_"+str(dim_screens[i])+"_L0_"+str(int(L0[i]))+".npy")
        file_istx=os.path.normpath(chakra_ao_savepath+"/turbu/istx_"+str(dim_screens[i])+"_L0_"+str(int(L0[i]))+".npy")
        file_isty=os.path.normpath(chakra_ao_savepath+"/turbu/isty_"+str(dim_screens[i])+"_L0_"+str(int(L0[i]))+".npy") 

        if(os.path.isfile(file_A) and os.path.isfile(file_B) and
           os.path.isfile(file_istx) and os.path.isfile(file_isty)):
            if(verbose<1):print "reading files"
            A=np.load(file_A)
            B=np.load(file_B)
            A_F=np.reshape(A.flatten("F"),(A.shape[0],A.shape[1]))
            B_F=np.reshape(B.flatten("F"),(B.shape[0],B.shape[1]))
            istx=np.load(file_istx).astype(np.uint32)
            isty=np.load(file_isty).astype(np.uint32)

        else:
            A,B,istx,isty=itK.AB(dim_screens[i],L0[i],verbose)
            if(verbose==0):
                print "writing files"
                np.save(file_A,A)
                np.save(file_B,B)
                np.save(file_istx,istx)
                np.save(file_isty,isty)
            A_F=np.reshape(A.flatten("F"),(A.shape[0],A.shape[1]))
            B_F=np.reshape(B.flatten("F"),(B.shape[0],B.shape[1]))

        tscreen=atmos_obj.s_a.d_screens[alt[i]]
        tscreen.init_screen(<float*>(A_F.data),<float*>(B_F.data),
                <unsigned int*>istx.data,<unsigned int*>isty.data,1234)
        for j in range(2*tscreen.screen_size):
            tscreen.extrude(1)
    return atmos_obj


cdef compute_size2(np.ndarray[ndim=1,dtype=np.int64_t]size):
    """Compute the size of a stencil, given the screen size
    
    :param size: (np.ndarray[ndim=1,dtype=np.int64_t]) :screen size
    """
    
    cdef n=size.shape[0]
    cdef np.ndarray[ndim=1,dtype=np.int64_t] size2=np.zeros(n,dtype=np.int64)

    cdef int i
    for i in range(n):
        size2[i]=itK.stencil_size(size[i])
    return size2


