import numpy as np
#################################################
# P-Class target
#################################################
cdef class Target:

    def __cinit__(self,naga_context ctxt, Telescope telescope, int ntargets, 
                    np.ndarray[ndim=1,dtype=np.float32_t] xpos,
                    np.ndarray[ndim=1,dtype=np.float32_t] ypos,
                    np.ndarray[ndim=1,dtype=np.float32_t] Lambda,
                    np.ndarray[ndim=1,dtype=np.float32_t] mag,
                    float zerop,
                    np.ndarray[ndim=1,dtype=np.int64_t] size,
                    int Npts,
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

        self.target= new sutra_target(ctxt.c, telescope.telescope, ntargets,
                    <float*>xpos.data,<float*>ypos.data,
                    <float*>Lambda.data,<float*>mag.data,
                    zerop, <long*>size.data,
                    Npts, device)



    def __dealloc__(self):
        del self.target

    def add_layer(self, int n, bytes l_type, float alt, float xoff, float yoff):
        """Add a phase screen dm or atmos as layers of turbulence

        :parameters:
            n: (int) : index of the target 
            l_type: (str) : "atmos" or "dm"

            alt: (float) : altitude

            xoff: (float) : x-offset

            yoff: (float) : y-offset
        """
        self.context.set_activeDevice(self.device)
        self.target.d_targets[n].add_layer(<char*>l_type, alt, xoff, yoff)


    def init_strehlmeter(self, int nTarget):
        """Initialise target's strehl

        :param nTarget: (int) : index of the target
        """
        cdef sutra_source *s_s_ptr
        self.target.d_targets[nTarget].init_strehlmeter()

    def atmos_trace(self, int nTarget, Atmos atm, Telescope tel):
        """Raytracing of the target through the atmosphere 

        :parameters:
            int: (nTarget)   : index of the target

            atm: (atmos)     : atmos to get through

            tel: (Telescope) : telescope
        """

        self.context.set_activeDevice(self.device)
        self.target.d_targets[nTarget].raytrace(atm.s_a)
        cdef carma_obj[float] *d_screen=self.target.d_targets[nTarget].d_phase.d_screen
        cdef carma_obj[float] *d_tel=tel.telescope.d_phase_ab_M1

        d_screen.axpy(1.0,d_tel,1,1)


    def get_image(self,int nTarget, bytes type_im, long puponly=0):
        """Return the image from the target (or long exposure image according to the requested type)
        :parameters:
            nTarget: (int) : index of the target

            type_im: (str) : type of the image to get ("se" or "le")

            puponly: (int) : if 1, image computed from phase on the pupil only
        """
        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source *src = self.target.d_targets[nTarget]

        cdef const long *dims
        dims=src.d_image.getDims()
        cdef np.ndarray data_F=np.empty((dims[2],dims[1]),dtype=np.float32)
        cdef np.ndarray data=np.empty((dims[1],dims[2]),dtype=np.float32)

        src.comp_image(puponly)
        if(type_im=="se"):
            src.d_image.device2host(<float*>data_F.data)
        elif(type_im=="le"):
            src.d_leimage.device2host(<float*>data_F.data)

        data=np.reshape(data_F.flatten("F"),(dims[1],dims[2]))
        data = np.roll(np.roll(data,data.shape[0]/2,axis=0),data.shape[1]/2,axis=1)
        return data


    def get_phase(self, int nTarget):
        """Return the phase's screen of the target

        :param nTarget: (int) : index of the target
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source *src = self.target.d_targets[nTarget]

        cdef const long *dims
        dims = src.d_phase.d_screen.getDims()
        cdef np.ndarray data_F=np.empty((dims[2],dims[1]),dtype=np.float32)
        cdef np.ndarray data=np.empty((dims[1],dims[2]),dtype=np.float32)
        src.d_phase.d_screen.device2host(<float*>data_F.data)

        data=np.reshape(data_F.flatten("F"),(dims[1],dims[2]))
        return data


    def get_phasetele(self, int nTarget):
        """Return the telemetry phase of the target

        :param nTarget: (int) : index of the target
        :return data: (np.ndarray(ndim=2,np.float32)) : phase screen
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source *src = self.target.d_targets[nTarget]

        cdef const long *dims
        dims=src.phase_telemetry.getDims()
        cdef np.ndarray data_F=np.empty((dims[2],dims[1]),dtype=np.float32)
        cdef np.ndarray data=np.empty((dims[1],dims[2]),dtype=np.float32)

        src.phase_telemetry.fill_into(<float*>data_F.data)

        data=np.reshape(data_F.flatten("F"),(dims[1],dims[2]))
        return data



    def get_amplipup(self, int nTarget):
        """Return the complex amplitude in the pupil plane of the target.

        :param nTarget: (int) : index of the target
        """

        self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source *src = self.target.d_targets[nTarget]
        
        cdef const long *dims
        dims=src.d_amplipup.getDims()

        cdef np.ndarray data_F=np.empty((dims[2],dims[1]),dtype=np.complex64)
        cdef np.ndarray data=np.empty((dims[1],dims[2]),dtype=np.complex64)
        src.d_amplipup.device2host(<cuFloatComplex*>data_F.data)

        data=np.reshape(data_F.flatten("F"),(dims[1],dims[2]))
        return data


    def get_strehl(self, int nTarget):
        """Return the target's strehl

        :param nTarget: (int) : index of the target
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


    cpdef dmtrace(self,int ntar, Dms dms, int reset=0):
        cdef carma_context *context=carma_context.instance()
        context.set_activeDevice(self.target.d_targets[ntar].device,1)

        self.target.d_targets[ntar].raytrace(dms.dms,reset)


    def __str__(self):
        info = "Target objects:\n"
        info+= "Contains "+str(self.target.ntargets)+" Target(s):\n"
        cdef int i
        cdef sutra_source *target
        info+= "Source # | position(\") |  Mag  | Lambda (mic.)\n"
        for i in range(self.target.ntargets):
            target=self.target.d_targets.at(i)
            info+= "%8d"%(i+1)+ " | "+ "%4d"%target.tposx+" , "+"%-4d"%target.tposy+\
                  " | "+"%5.2f"%target.mag+" | "+"%5.3f"%target.Lambda+"\n"
        info+= "--------------------------------------------------------"
        return info


def target_init(naga_context ctxt, Telescope telescope, Param_target p_target, Param_atmos atm,
                Param_geom geom, Param_tel tel,wfs=None, 
                Sensors sensors=None,
                dm=None):
    """Create a cython target from parametres structures
    
    :parameters:
        ctxt: (naga_context) :

        atm: (Param_atmos) : atmos settings

        geom: (Param_geom) : geom settings

        wfs: (Param_wfs) : wfs settings

        dm: (Param_dm) : dm settings
    """
    cdef bytes type_target=bytes("atmos")

    cdef Target target
    cdef float xoff
    cdef float yoff
    cdef int Npts

    cdef int i, j, k

    if(p_target.ntargets>0):
        if(p_target.dms_seen is None):
            for i in range(p_target.ntargets):
                if(dm is not None):
                    p_target.dms_seen=np.arange(len(dm))

    cdef np.ndarray sizes= np.ones(p_target.ntargets,dtype=np.int64)*geom.pupdiam

    #sizes = sizes(-::y_target.ntargets-1);
# ATTEMPT AT MAKING IT POSSIBLE TO WORK WITH NON BINARY PUPILS
# switching *y_geom._spupil and *y_geom._apodizer with ceil(*y_geom._spupil) and ceil(*y_geom._apodizer)
    ceiled_pupil=np.ceil(geom._spupil)

    ceiled_pupil[np.where(ceiled_pupil>1)]=1;

    if(p_target.apod==1):
            Npts=0
            #TODO apodizer, Npts=nb element of apodizer>0
            ceiled_apodizer=np.ceil(geom._apodizer*geom._spupil)
            ceiled_apodizer[np.where(ceiled_apodizer>1)]=1
            target = Target(ctxt, telescope, p_target.ntargets,p_target.xpos,p_target.ypos,
                        p_target.Lambda, p_target.mag,p_target.zerop,sizes,Npts)
            #TODO if brama ...
    else:
            Npts=np.sum(ceiled_pupil)
            target= Target(ctxt, telescope, p_target.ntargets,p_target.xpos,p_target.ypos,
                            p_target.Lambda, p_target.mag,p_target.zerop,sizes,Npts)
            #TODO if brama ...

    cdef long dims,dim, dim_dm

    #cc=i
    for i in range(p_target.ntargets):
        if(atm.nscreens>0):
            for j in range(atm.nscreens):
                xoff=p_target.xpos[i]*4.848e-6*atm.alt[j]/atm.pupixsize
                yoff=p_target.ypos[i]*4.848e-6*atm.alt[j]/atm.pupixsize
                xoff+=float((atm.dim_screens[j]-geom._n)/2)
                yoff+=float((atm.dim_screens[j]-geom._n)/2)
                pupdiff = (geom._n - geom.pupdiam)/2
                xoff += pupdiff
                yoff += pupdiff
                target.add_layer(i,type_target,atm.alt[j],xoff,yoff)

        #if (y_dm != []) {
        if(dm is not None):
            #j=ddd
            #for (ddd=1;ddd<=numberof(*y_target(cc).dms_seen);ddd++) {
            for j in range(p_target.dms_seen.size):
                #k=dd
                #dd = (*y_target(cc).dms_seen)(ddd)
                k=p_target.dms_seen[j]
                dims=dm[k]._n2-dm[k]._n1+1
                dim=geom._mpupil[2].size
                dim_dm = max(dim,dims)
                xoff=target.xpos[i]*4.848e-6*dm[k].alt/tel.diam*geom.pupdiam
                yoff=target.ypos[i]*4.848e-6*dm[k].alt/tel.diam*geom.pupdiam


                xoff+=float((dim_dm-geom._n)/2)
                yoff+=float((dim_dm-geom._n)/2)

                pupdiff=(geom._n-geom.pupdiam)/2
                xoff+=pupdiff
                yoff+=pupdiff

                if(dm[k].type_dm=="kl"):
                    xoff+=2
                    yoff+=2
                target.add_layer(i,dm[k].type_dm,dm[k].alt,xoff,yoff)

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
                        yoff=(gsalt * atm.alt[j] * (tel.diam/2.) + wfs[i].ypos*4.848e-6*atm.alt[j])/atm.pupixsize
                        xoff=xoff+(atm.dim_screens[j]-geom._n)/2
                        yoff=yoff+(atm.dim_screens[j]-geom._n)/2
                        sensors.sensors.d_wfs[i].d_gs.add_layer(type_target, atm.alt[j],xoff,yoff)

                if(dm is not None and not wfs[i].openloop):
                    for j in range(wfs[i].dms_seen.size):
                        k=wfs[i].dms_seen[j]
                        dims=dm[k]._n2-dm[k]._n1+1
                        dim=geom._mpupil.shape[0]
                        if(dim<dims):
                            dim=dims
                        xoff=wfs[i].xpos*4.848e-6*dm[k].alt/tel.diam*geom.pupdiam
                        yoff=wfs[i].ypos*4.848e-6*dm[k].alt/tel.diam*geom.pupdiam
                        xoff=xoff+(dim-geom._n)/2
                        yoff=yoff+(dim-geom._n)/2
                        sensors.sensors.d_wfs[i].d_gs.add_layer(dm[k].type_dm, dm[k].alt,xoff,yoff)



    return target

