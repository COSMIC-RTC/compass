
#################################################
# P-Class (parametres) Param_target
#################################################
cdef class Param_target:
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

    def target_init(self,chakra_context ctxt, Param_atmos atm,
                    Param_geom geom, Param_tel tel,wfs=None, 
                    Sensors sensors=None,
                    param_dm=None):
        """Create a cython target from parametres structures
        

        ctxt    -- chakra_context
        atm     -- Param_atmos
        geom    -- Param_geom
        wfs     -- Param_wfs
        dm      -- Param_dm
        """
        cdef bytes type_target=bytes("atmos")

        cdef Target target
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
                target = Target(ctxt, self.ntargets,self.xpos,self.ypos,
                            self.Lambda, self.mag,sizes,ceiled_apodizer)
        else:
                target= Target(ctxt, self.ntargets,self.xpos,self.ypos,
                                self.Lambda, self.mag,sizes,ceiled_pupil)

        cdef int i
        cdef int j
        for i in range(self.ntargets):
            if(atm.nscreens>0):
                for j in range(atm.nscreens):
                    xoff=self.xpos[i]*4.848e-6*atm.alt[j]/atm.pupixsize
                    yoff=self.ypos[i]*4.848e-6*atm.alt[j]/atm.pupixsize
                    target.add_layer(i,type_target,atm.alt[j],xoff,yoff)
            #if(param_dm is not None):
                """TODO dm !=None
                cf yorick func target_init,  file yoga_ao/yorick/yoga_ao.i  l 369"""

            target.init_strehlmeter(i)

        if(wfs is not None):
            if(sensors is not None):
                for i in range(len(wfs)):
                    if(wfs[i].gsalt!=0):
                        gsalt=1/wfs[i].gsalt
                    else:
                         gsalt=0
                    
                    if(wfs[i].atmos_seen):
                        for j in range(atm.nscreens):
                            xoff=(gsalt * atm.alt[j] * (tel.diam/2.) + wfs[i].xpos*4.848e-6*atm.alt[j])/atm.pupixsize
                            yoff = (gsalt * atm.alt[j] * (tel.diam/2.) + wfs[i].ypos*4.848e-6**atm.alt[j])/atm.pupixsize
                            xoff = xoff+(atm.dim_screens[j]-geom._n)/2
                            yoff = yoff+(atm.dim_screens[j]-geom._n)/2
                            sensors.sensors.d_wfs[i].d_gs.add_layer(type_target, atm.alt[j],xoff,yoff)
                            print "atmos_seen == true"
        return target

"""
  if (y_wfs != []) {
    if ((y_wfs != []) && (g_wfs != [])) {
      for (cc=1;cc<=numberof(y_wfs);cc++) {
        if (y_wfs(cc).gsalt != 0.)
          gsalt = 1./y_wfs(cc).gsalt;
        else gsalt = 0.;
        if (y_atmos != [] && y_wfs(cc).atmos_seen) {
          for (dd=1;dd<=y_atmos.nscreens;dd++) {
            xoff = (gsalt * (*y_atmos.alt)(dd) * (y_tel.diam/2.) + (y_wfs.xpos)(cc)*4.848e-6*(*y_atmos.alt)(dd))/y_atmos.pupixsize;
            yoff = (gsalt * (*y_atmos.alt)(dd) * (y_tel.diam/2.) + (y_wfs.ypos)(cc)*4.848e-6*(*y_atmos.alt)(dd))/y_atmos.pupixsize;
            xoff = float(xoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
            yoff = float(yoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
            sensors_addlayer,g_wfs,cc-1,type,(*y_atmos.alt)(dd),xoff,yoff;
          }
        }
        if (y_dm != [] && !y_wfs(cc).openloop) {
          for (ddd=1;ddd<=numberof(*y_wfs(cc).dms_seen);ddd++) {
	    dd = (*y_wfs(cc).dms_seen)(ddd);
            dims = y_dm(dd)._n2 - y_dm(dd)._n1 + 1;
            dim  = dimsof(*y_geom._mpupil)(2);
            dim_dm = max([dim,dims]);
            xoff = (y_wfs.xpos)(cc)*4.848e-6*(y_dm.alt)(dd)/(y_tel.diam / y_geom.pupdiam);
            yoff = (y_wfs.ypos)(cc)*4.848e-6*(y_dm.alt)(dd)/(y_tel.diam / y_geom.pupdiam);
            xoff = float(xoff+(dim_dm-y_geom._n)/2);
            yoff = float(yoff+(dim_dm-y_geom._n)/2);
            sensors_addlayer,g_wfs,cc-1,y_dm(dd).type,(y_dm.alt)(dd),xoff,yoff;
          }
        }
      }
    }
  }
"""





#################################################
# P-Class target
#################################################
cdef class Target:

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
        """TODO doc

        n       -- int   : 
        l_type  -- str   :
        alt     -- float :
        xoff    -- float :
        yoff    -- float :
        """
        self.context.set_activeDevice(self.device)
        self.target.d_targets[n].add_layer(<char*>l_type, alt, xoff, yoff)


    def init_strehlmeter(self, int nTarget):
        """TODO doc

        nTarget -- int :
        """
        cdef sutra_source *s_s_ptr
        self.target.d_targets[nTarget].init_strehlmeter()

    def atmos_trace(self, int nTarget, Atmos atm):
        """ TODO doc

        int --  nTarget :
        atm --  atmos
        """
        self.context.set_activeDevice(self.device)
        self.target.d_targets[nTarget].raytrace(atm.s_a)

    def get_image(self,int nTarget, bytes type_im, long puponly=0):
        """TODO doc

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
        """TODO doc

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
        """TODO doc

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
        """TODO doc

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
        """TODO doc

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


