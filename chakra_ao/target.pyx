
#################################################
# P-Class (parametres) param_target
#################################################
cdef class param_target:
    def __cinit__(self):
        self.ntargets=0


    def set_nTargets(self, int n):
        """Set the attribute to self.n to n
        n -- int : number of targets."""
        self.ntargets=n
    def set_apod(self, int a):
        """Set the attribute apod to a
        a -- int"""
        self.apod=a
    def set_Lambda(self, l):
        """Set the attribute to
        l -- list of float""" 
        self.Lambda=np.array(l,dtype=np.float32)
    def set_xpos(self,l):
        """l -- list of float"""
        self.xpos =np.array(l,dtype=np.float32)
    def set_ypos(self,l):
        """l -- list of float"""
        self.ypos =np.array(l,dtype=np.float32)
    def set_mag(self,l):
        """l -- list of float"""
        self.mag =np.array(l,dtype=np.float32)

    def target_init(self,chakra_context ctxt, param_atmos atm,
                    param_geom geom, param_wfs=None,param_dm=None):
        """Create a cython target from parametres structures
        

        ctxt    -- chakra_context
        atm     -- param_atmos
        geom    -- param_geom
        wfs     -- param_wfs
        dm      --param_dm
        """
        cdef str type_target="atmos"

        cdef target Target
        cdef float xoff
        cdef float yoff
        #if (y_target != []) {
        cdef np.ndarray sizes= np.ones(self.ntargets,dtype=np.int64)*geom.pupdiam
        #sizes = sizes(-::y_target.ntargets-1);
# ATTEMPT AT MAKING IT POSSIBLE TO WORK WITH NON BINARY PUPILS
# switching *y_geom._spupil and *y_geom._apodizer with ceil(*y_geom._spupil) and ceil(*y_geom._apodizer)
        cdef np.ndarray ceiled = np.empty((sizes,sizes),dtype=np.int32)
        ceiled_pupil=np.ceil(geom._spupil)

        ceiled_pupil[np.where(ceiled_pupil>1)]=1;

        if(self.apod==1):
                ceiled_apodizer=np.ceil(geom._apodizer*geom._spupil)
                ceiled_apodizer[np.where(ceiled_apodizer>1)]=1
                Target = target(ctxt, self.ntargets,self.xpos,self.ypos,
                            self.Lambda, self.mag,sizes,ceiled_apodizer)
        else:
                Target= target(ctxt, self.ntargets,self.xpos,self.ypos,
                                self.Lambda, self.mag,sizes,ceiled_pupil)

        cdef int i
        cdef int j
        for i in range(self.ntargets):
            if(atm.nscreens>0):
                for j in range(atm.nscreens):
                    xoff=self.xpos[i]*4.848e-6*atm.alt[j]/atm.pupixsize
                    yoff=self.ypos[i]*4.848e-6*atm.alt[j]/atm.pupixsize
                    Target.add_layer(i,type_target,atm.alt[j],xoff,yoff)
            #if(param_dm is not None):
                """TODO
                cf yorick func target_init,  file yoga_ao/yorick/yoga_ao.i  l 369"""

            Target.init_strehlmeter(i)

        return Target



#################################################
# P-Class target
#################################################
cdef class target:

    def __cinit__(self,chakra_context ctxt, int ntargets, 
                    np.ndarray[dtype=np.float32_t] xpos,
                    np.ndarray[dtype=np.float32_t] ypos,
                    np.ndarray[dtype=np.float32_t] Lambda,
                    np.ndarray[dtype=np.float32_t] mag,
                    np.ndarray[dtype=np.int64_t] size,
                    np.ndarray[ndim=2, dtype=np.float32_t] pupil,
                    int device=-1
                    ):

        if(device<0):
            device = ctxt.get_activeDevice()
        else:
            ctxt.set_activedevice(device)

        self.ntargets=ntargets
        self.xpos=xpos
        self.ypos=ypos
        self.Lambda=Lambda
        self.mag=mag
        self.device=device
        self.context=ctxt


        self.target= new sutra_target(ctxt.c,ntargets,
                    <float*>xpos.data,<float*>ypos.data,
                    <float*>Lambda.data,<float*>mag.data,
                    <long*>size.data,<float*>pupil.data,device)


    def add_layer(self, int n, bytes l_type, float alt, float xoff, float yoff):
        """TODO

        n       -- int   : 
        l_type  -- str   :
        alt     -- float :
        xoff    -- float :
        yoff    -- float :
        """
        self.context.set_activeDevice(self.device)
        self.target.d_targets[n].add_layer(<char*>l_type, alt, xoff, yoff)


    def init_strehlmeter(self, int nTarget):
        """TODO

        nTarget -- int :
        """
        cdef sutra_source *s_s_ptr
        self.target.d_targets[nTarget].init_strehlmeter()

    def atmos_trace(self, int nTarget, atmos atm):
        """ TODO

        int --  nTarget :
        atm --  atmos
        """
        self.context.set_activeDevice(self.device)
        self.target.d_targets[nTarget].raytrace(atm.s_a)

    def get_image(self,int nTarget, bytes type_im, long puponly=0):
        """TODO

        nTarget -- int :
        type_im -- str :
        puponly -- int : (default=0)
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source *src = self.target.d_targets[nTarget]

        cdef const long *dims
        dims=src.d_image.getDims()
        cdef np.ndarray data=np.empty((dims[1],dims[2]),dtype=np.float32)

        src.comp_image(puponly)
        if(type_im=="se"):
            src.d_image.device2host(<float*>data.data)
        elif(type_im=="le"):
            src.d_leimage.device2host(<float*>data.data)

        return data


    def get_phase(self, int nTarget):
        """TODO

        nTarget -- int :
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source *src = self.target.d_targets[nTarget]

        cdef const long *dims
        dims = src.d_phase.d_screen.getDims()
        cdef np.ndarray data=np.empty((dims[1],dims[2]),dtype=np.float32)
        src.d_phase.d_screen.device2host(<float*>data.data)

        return data


    def get_phasetele(self, int nTarget):
        """TODO

        nTarget -- int :
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source *src = self.target.d_targets[nTarget]

        cdef const long *dims
        dims=src.phase_telemetry.getDims()
        cdef np.ndarray data=np.empty((dims[1],dims[2]),dtype=np.float32)

        src.phase_telemetry.fill_into(<float*>data.data)

        return data



    def get_amplipup(self, int nTarget):
        """TODO

        nTarget -- int :
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source *src = self.target.d_targets[nTarget]
        
        cdef const long *dims
        dims=src.d_amplipup.getDims()

        cdef np.ndarray data=np.empty((dims[1],dims[2]),dtype=np.complex64)
        src.d_amplipup.device2host(<cuFloatComplex*>data.data)

        return data


    def get_strehl(self, int nTarget):
        """TODO

        nTarget -- int :
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source *src = self.target.d_targets[nTarget]

        cdef np.ndarray strehl=np.empty(4,dtype=np.float32)
        src.comp_image(0)
        src.comp_strehl()

        strehl[0]=src.strehl_se
        strehl[1]=src.strehl_le
        strehl[2]=src.phase_var
        strehl[3]=src.phase_var_avg/(src.phase_var_count + 0.00001)

        return strehl


