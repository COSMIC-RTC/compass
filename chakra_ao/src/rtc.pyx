cdef class Param_rtc:
    def set_nwfs( self, n):
        """Set the number of wfs

        :param n: (int) number of wfs
        """
        self.nwfs=n
    def set_centroiders(self,l):
        """Set the centroiders

        :param l: (list of Param_centroider) : centroiders settings
        """
        self.centroiders=l

    def set_controllers(self,l):
        """Set the controller

        :param l: (list of Param_controller) : controllers settings
        """
        self.controllers=l


cdef class Rtc:
    def __cinit__(self, device=-1):
        cdef carma_context *context =carma_context.instance()
        if(device==-1):
            device=context.get_activeDevice()
        context.set_activeDevice(device,1)
        self.use_brama=0
        self.device=device

        self.rtc=new sutra_rtc(context)


    def add_centroider(self,Sensors sensor, long nwfs, long nvalid,
                        bytes type_centro, float offset, float scale):

        """TODO doc

        :parameters:
            sensor: (Sensors) :

            nwfs : (long) : number of wfs

            nvalid: (long) : number of valid subaps

            type_centro: (str) : centroider's type

            offset: (float) : 

            scale: (float) :

        """
        cdef carma_context *context = carma_context.instance()
        cdef int activeDevice=self.rtc.device
        context.set_activeDeviceForCpy(self.rtc.device,1)
        self.rtc.add_centroider(sensor.sensors,nwfs,nvalid,offset,scale,activeDevice,type_centro)



    cdef add_controller(self, int nactu, float delay, bytes type_control, Dms dms,
                 char **type_dmseen, np.ndarray[ndim=1,dtype=np.float32_t] alt,
                 int ndm, long Nphi=-1):
        """TODO doc

        :parameters:
            nactu: (int) : number of controled actuator

            delay: (float) :

            type_control: (str) : controller's type

            dms: (Dms) : 

            type_dmseen: (char**) : 

            alt: (np.ndarray[ndim=1,dtype=np.float32_t]) : 

            ndm: (int) :

            Nphi: (long) :
        """

        cdef float *ptr_alt    =<float*>alt.data
        cdef char *ptr_dmseen =<char *> type_dmseen
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.device,1)
        if(Nphi>-1):
            self.rtc.add_controller_geo(nactu,Nphi,delay,self.device,dms.dms,&ptr_dmseen,ptr_alt,ndm)
        else:
            self.rtc.add_controller(nactu,delay,self.device,type_control,dms.dms,&ptr_dmseen,ptr_alt,ndm)



    def rmcontrol(self):
        """TODO doc"""
        cdef carma_context *context=carma_context.instance()
        context.set_activeDevice(self.device,1)
        self.dms.rm_controller()


    def setthresh(self, int ncentro, float thresh):
        """set threshold

        :parameters:
            ncentro: (int) : centroider's index

            thresh: (float) : threshold
        """

        cdef carma_context *context=carma_context.instance()
        cdef sutra_centroider_tcog *tcog=NULL
        cdef sutra_centroider *centro=NULL
        context.set_activeDeviceForCpy(self.device,1)
        if(self.dms.d_centro.at(ncentro).is_type("tcog")):
            centro=self.rtc.d_centro.at(ncentro)
            tcog=dynamic_cast_centroider_tcog_ptr(centro)
            tcog.set_threshold(thresh)


    def setnmax(self,int ncentro, int nmax):
        """TODO doc

        :parameters:
            ncentro: (int) : centroider's index

            nmax: (int) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDevice(self.device,1)

        #TODO add centroider_bpcog, cast to centroiderbpcog
        #cdef sutra_centroider_bpcog *bpcog=NULL
        #cdef sutra_centroider *centro=NULL

        #if(self.d_centro.at(ncentro).is_type("bpcog")):
        #    centro=self.rtc.d_centro.at(ncentro)
        #    bpcog=dynamic_cast_centroider_tcog_ptr(centro)
        #    bpcog.setnmax(nmax)

    def sensors_initbcube(self,int ncentro):
        """TODO doc

        :param ncentro: (int) : centroider's index
        """
            
        cdef carma_context *context=carma_context.instance()
        context.set_activeDevice(self.device,1)
        cdef sutra_centroider *centro=NULL
        centro=self.rtc.d_centro.at(ncentro)
        cdef sutra_centroider_corr *corr
        if(self.rtc.d_centro[ncentro].is_type(<bytes>"corr")):
            corr= dynamic_cast_centroider_corr_ptr( centro)
            corr.init_bincube()
        return 1


    def sensors_initweight(self,int ncentro, np.ndarray[ndim=1, dtype=np.float32_t] w):
        """Set weigth of a centroider

        :parameters:
            ncentro: (int) : centroider's index

            w: (np.ndarray[ndim=1, dtype=np.float32_t]) : weight
        """
        cdef carma_context *context = carma_context.instance()
        context.set_activeDevice(self.rtc.device,1)
        cdef sutra_centroider_wcog *centro=\
            dynamic_cast_centroider_wcog_ptr(self.rtc.d_centro[ncentro])
        centro.init_weights()
        cdef np.ndarray w_F=w.flatten("F")
        centro.load_weights(<float*>w_F.data,int(w.ndim))


    def sensors_initcorr(self,int ncentro, np.ndarray[ndim=1,dtype=np.float32_t] w,
                        np.ndarray[ndim=2,dtype=np.float32_t] corr_norm,
                        int sizex, int sizey,
                        np.ndarray[ndim=2,dtype=np.float32_t] interpmat):
        """TODO doc

        :parameters:
            ncentro: (int) : centroider's index

            w: (np.ndarray[ndim=1,dtype=np.float32_t]) : weight

            corr_norm: (np.ndarray[ndim=2,dtype=np.float32_t]) :

            sizex: (int) :

            sizey: (int) :

            interpmat: ([ndim=2,dtype=np.float32_t]) :

        """

        cdef carma_context *context = carma_context.instance()
        context.set_activeDevice(self.rtc.device,1)

        cdef sutra_centroider_corr * centro_corr =dynamic_cast_centroider_corr_ptr(\
                self.rtc.d_centro[ncentro]) 
        cdef np.ndarray w_F=w.flatten("F")
        cdef np.ndarray[dtype=np.float32_t] corr_norm_F=corr_norm.flatten("F")
        cdef np.ndarray[dtype=np.float32_t] interpmat_F=interpmat.flatten("F")
        centro_corr.init_corr(sizex,sizey,<float*>interpmat_F.data)
        centro_corr.load_corr(<float*>w_F.data,<float*>corr_norm_F.data,int(w.ndim))

    #TODO possible error -> check it
    cpdef getcentroids(self,int ncontrol, Sensors g_wfs=None, int nwfs=0):
        """TODO doc

        :parameters:
            ncontrol: (int) : controller's index

            g_wfs: (Sensors) : (optional)

            nwfs: (int) : (optional) number of wfs
        """

        cdef const long *dims
        cdef sutra_wfs *wfs
        cdef sutra_rtc *rtc = self.rtc
        cdef carma_obj[float] *d_tmp
        cdef carma_obj[float] *d_data
        cdef carma_context *context=carma_context.instance()
        cdef np.ndarray[ndim=1,dtype=np.float32_t] data


        if(ncontrol>=abs(rtc.d_control.size())):
            if(g_wfs is None):
                raise ValueError("Controller not initialized on the GPU, you have to specify the WFS")
            wfs=g_wfs.sensors.d_wfs[nwfs]
            dims=wfs.d_slopes.getDims()
            data=np.zeros((dims[1]),dtype=np.float32)
            d_tmp= new carma_obj[float](context,wfs.d_subsum.getDims())
            d_data=new carma_obj[float](context,dims)
            rtc.d_centro[ncontrol].get_cog(d_tmp.getData(),d_data.getData())
            d_data.device2host(<float*>data.data)

        else:
            dims=rtc.d_control[ncontrol].d_centroids.getDims()
            data=np.zeros((dims[1]),dtype=np.float32)
            rtc.d_control[ncontrol].d_centroids.device2host(<float*>data.data)
        return data


    cpdef docentroids(self,int ncontrol=-1):
        """TODO doc

        :param ncontrol: (optional)
        """
        if(ncontrol>-1):
            self.rtc.do_centroids(ncontrol)
        else:
            self.rtc.do_centroids()
        

    cdef init_proj(self,int ncontro,Dms dms,np.ndarray[ndim=1,dtype=np.int32_t] indx_dm,
            np.ndarray[ndim=1,dtype=np.float32_t] unitpervolt, np.ndarray[ndim=1,dtype=np.int32_t] indx_pup):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            dms: (Dms) :

            indx_dm: (np.ndarray[ndim=1,dtype=np.int32_t]) :

            unitpervolt: (np.ndarray[ndim=1,dtype=np.float32_t]) :

            indx_pup: (np.ndarray[ndim=1,dtype=np.int32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        cdef sutra_controller_geo *controller_geo=\
            dynamic_cast_controller_geo_ptr(self.rtc.d_control[ncontro])

        context.set_activeDeviceForCpy(self.rtc.device,1)
        controller_geo.init_proj_sparse(dms.dms,<int*>indx_dm.data,
            <float*>unitpervolt.data, <int*>indx_pup.data)


    cdef init_modalOpti(self,int ncontro,int nmodes,int nrec, np.ndarray[ndim=2,dtype=np.float32_t] M2V,
            float gmin, float gmax, int ngain, float Fs):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            nmodes: (int) : number of modes

            nrec: (int) : number of recorded open slopes measurements

            M2V: (np.ndarray[ndim=2,dtype=np.float32_t]) :

            gmin: (float) : minimum gain for modal optimization

            gmax: (float) : maximum gain for modal optimization

            ngain: (int) : Number of tested gains

            Fs: (float) : sampling frequency
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls * controller_ls

        cdef np.ndarray[dtype=np.float32_t] M2V_F=M2V.flatten("F")

        if(<bytes>self.d_control[ncontro].get_type()=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.init_modalOpti(nmodes, nrec, <float*>M2V_F.data, gmin, gmax, ngain, Fs)
        else:
            raise TypeError("**** ERROR : Modal Optimization only for controller type ls ****")


    cdef loadOpenLoop(self,int ncontro, np.ndarray[ndim=2, dtype=np.float32_t] ol_slopes):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            pl_slopes: (np.ndarray[ndim=2, dtype=np.float32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls * controller_ls

        cdef np.ndarray[dtype=np.float32_t] slopes_F = ol_slopes.flatten("F")

        if(<bytes>self.rtc.d_control[ncontro].get_type()=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.loadOpenLoopSlp(<float*>slopes_F.data)
        else:
            raise TypeError("Controller type must be ls")

    cdef modalControlOptimization(self,int ncontro):
        """TODO doc

        :param ncontro: controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls * controller_ls

        if(<bytes>self.rtc.d_control[ncontro].get_type()=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            if(controller_ls.is_modopti):
                controller_ls.modalControlOptimization()
            else:
                raise ValueError("**** ERROR : Modal Optimization not initialized ****")
        else:
            raise TypeError("**** ERROR : Modal Optimization only for controller type ls ***")

    cdef set_gain(self, int ncontro, float gain):
        """TODO doc
        
        :parameters:
            ncontro: (int) : controller index

            gain: (float) 
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef bytes type_control = <bytes>self.rtc.d_control[ncontro].get_type()
        cdef sutra_controller_ls        *controller_ls
        cdef sutra_controller_mv        *controller_mv
        cdef sutra_controller_cured     *controller_cured
        cdef sutra_controller_geo       *controller_geo
        cdef sutra_controller_kalman    *controller_kl

        if(type_control=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.set_gain(gain)

        elif(type_control=="mv"):
            controller_mv =dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            controller_mv.set_gain(gain)

        elif(type_control=="cured"):
            controller_cured =dynamic_cast_controller_cured_ptr(self.rtc.d_control[ncontro])
            controller_cured.set_gain(gain)

        elif(type_control=="geo"):
            controller_geo =dynamic_cast_controller_geo_ptr(self.rtc.d_control[ncontro])
            controller_geo.set_gain(gain)
        elif(type_control=="kalman_CPU" or type_control=="kalman_GPU" or 
             type_control=="kalman_uninitialized"):
            controller_kl =dynamic_cast_controller_kl_ptr(self.rtc.d_control[ncontro])
            controller_kl.set_gain(gain)


    cdef set_mgain(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] mgain):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            mgain: (np.ndarray[ndim=1,dtype=np.float32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef sutra_controller_ls *controller_ls
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.set_mgain(<float*>mgain.data)
        elif(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            controller_mv.set_mgain(<float*>mgain.data)
        else:
            raise TypeError("Controller needs to be ls or mv")


    cdef set_imat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            data: (np.ndarray[ndim=2,dtype=np.float32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef np.ndarray[dtype=np.float32_t] data_F=data.flatten("F")

        cdef sutra_controller_ls *controller_ls
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.d_imat.host2device(<float*>data_F.data)
        elif(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            controller_mv.d_imat.host2device(<float*>data_F.data)
        else:
            raise TypeError("Controller needs to be ls or mv")


    cpdef get_imat(self, int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        
        cdef sutra_controller_ls *controller_ls
        cdef sutra_controller_mv *controller_mv
        cdef sutra_controller_cured *controller_cured
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()

        cdef const long *dims=NULL
        cdef np.ndarray[ndim=2,dtype=np.float32_t] imat_F
        cdef np.ndarray[ndim=2,dtype=np.float32_t] imat

        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            dims=controller_ls.d_imat.getDims()
            imat_F=np.zeros((dims[2],dims[1]),dtype=np.float32)
            controller_ls.d_imat.device2host(<float*>imat_F.data)

        elif(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            dims=controller_mv.d_imat.getDims()
            imat_F=np.zeros((dims[2],dims[1]),dtype=np.float32)
            controller_mv.d_imat.device2host(<float*>imat_F.data)

        elif(type_contro=="cured"):
            controller_cured=dynamic_cast_controller_cured_ptr(self.rtc.d_control[ncontro])
            dims=controller_cured.d_imat.getDims()
            imat_F=np.zeros((dims[2],dims[1]),dtype=np.float32)
            controller_cured.d_imat.device2host(<float*>imat_F.data)

        imat=np.reshape(imat_F.flatten("F"),(dims[1],dims[2]))
        return imat


    cdef set_cmat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            data: (np.ndarray[ndim=2,dtype=np.float32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef np.ndarray[dtype=np.float32_t] data_F=data.flatten("F")

        cdef sutra_controller_ls *controller_ls
        cdef sutra_controller_mv *controller_mv
        cdef sutra_controller_generic *controller_generic
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.d_cmat.host2device(<float*>data_F.data)
        elif(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            controller_mv.d_cmat.host2device(<float*>data_F.data)
        elif(type_contro=="generic"):
            controller_generic=dynamic_cast_controller_generic_ptr(self.rtc.d_control[ncontro])
            controller_generic.d_cmat.host2device(<float*>data_F.data)
        else:
            raise TypeError("Controller needs to be ls, mv or generic")


    cpdef get_cmat(self,int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls *controller_ls
        cdef sutra_controller_mv *controller_mv
        cdef sutra_controller_generic *controller_generic
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data_F ,cmat
        cdef const long *cdims
        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            cdims=controller_ls.d_cmat.getDims()
            data_F=np.zeros((cdims[2],cdims[1]),dtype=np.float32)
            controller_ls.d_cmat.device2host(<float*>data_F.data)
            cmat=np.reshape(data_F.flatten("F"),(cdims[1],cdims[2]))
            return cmat
        elif(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            cdims=controller_mv.d_cmat.getDims()
            data_F=np.zeros((cdims[2],cdims[1]),dtype=np.float32)
            controller_mv.d_cmat.device2host(<float*>data_F.data)
            cmat=np.reshape(data_F.flatten("F"),(cdims[1],cdims[2]))
            return cmat
        elif(type_contro=="generic"):
            controller_generic=dynamic_cast_controller_generic_ptr(self.rtc.d_control[ncontro])
            cdims=controller_generic.d_cmat.getDims()
            data_F=np.zeros((cdims[2],cdims[1]),dtype=np.float32)
            controller_generic.d_cmat.device2host(<float*>data_F.data)
            cmat=np.reshape(data_F.flatten("F"),(cdims[1],cdims[2]))
            return cmat
        else:
            raise TypeError("Controller needs to be ls, mv or generic")



    cdef set_decayFactor(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] decay):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            decay: (np.ndarray[ndim=1,dtype=np.float32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_generic *controller_generic
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="generic"):
            controller_generic=dynamic_cast_controller_generic_ptr(self.rtc.d_control[ncontro])
            controller_generic.set_decayFactor(<float*>decay.data)
        else:
            raise TypeError("Controller needs to be generic")


    cdef set_matE(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] matE):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            matE: (np.ndarray[ndim=2,dtype=np.float32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_generic *controller_generic
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()

        cdef np.ndarray[dtype=np.float32_t] matE_F=matE.flatten("F")

        if(type_contro=="generic"):
            controller_generic=dynamic_cast_controller_generic_ptr(self.rtc.d_control[ncontro])
            controller_generic.set_matE(<float*>matE_F.data)
        else:
            raise TypeError("Controller needs to be generic")


    cdef doimat_geom(self, int ncontro, Dms g_dms,int geom):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            g_dms: (Dms) :

            geom: (int) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        self.rtc.do_imat_geom(ncontro,g_dms.dms,geom)
        print "TODO call imat_geom"
        print "TODO set_imat"

    cdef doimat(self, int ncontro, Dms g_dms):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            g_dms: (Dms) :
        """

        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)


        cdef sutra_controller *control=self.rtc.d_control[ncontro]
        cdef carma_obj[float] *d_imat=NULL

        if(<bytes>control.get_type()==<bytes>"ls"):
            d_imat=dynamic_cast_controller_ls_ptr(control).d_imat
        elif(<bytes>control.get_type()==<bytes>"mv"):
            d_imat=dynamic_cast_controller_mv_ptr(control).d_imat

        cdef int inds1,j,idx_cntr, device
        cdef float tmp_noise
        inds1=0
        cdef sutra_dm *dm
        cdef sutra_wfs *wfs
        cdef carma_obj[float] *screen
        cdef vector[sutra_dm *].iterator it_dm

        cdef float *d_centroids
        cdef np.ndarray[ndim=1,dtype=np.float32_t] h_centroids

        cdef int rank
        IF USE_MPI==1:
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD,&rank)
        ELSE:
            rank=0

        it_dm=control.d_dmseen.begin()
        while(it_dm!=control.d_dmseen.end()):
            dm=deref(it_dm)
            for j in range(dm.ninflu):
                dm.comp_oneactu(j,dm.push4imat)
                for idx_cntr in range(<int>self.rtc.d_centro.size()):
                    wfs=self.rtc.d_centro[idx_cntr].wfs
                    screen=wfs.d_gs.d_phase.d_screen
                    tmp_noise=wfs.noise
                    wfs.noise=-1
                    wfs.kernconv=True
                    wfs.sensor_trace(g_dms.dms,1)
                    IF USE_MPI==1:
                        Bcast(screen,0)
                    wfs.comp_image()
                    wfs.noise=tmp_noise
                    wfs.kernconv=False

                self.rtc.do_centroids(ncontro,True)

                h_centroids=self.getCentroids(ncontro)
                control.d_centroids.host2device(<float*>h_centroids.data)

                device=control.d_centroids.getDevice()
                d_centroids=control.d_centroids.getData()

                convert_centro(d_centroids,
                        d_centroids,0,
                        (0.5/dm.push4imat),
                        control.d_centroids.getNbElem(),
                        context.get_device(device))

                control.d_centroids.copyInto(
                        &d_imat.getData()[inds1],
                        control.nslope())
                dm.reset_shape()

                dm.comp_oneactu(j,-1.*dm.push4imat)

                for idx_cntr in range(<int>self.rtc.d_centro.size()):
                    wfs=self.rtc.d_centro[idx_cntr].wfs
                    tmp_noise=wfs.noise
                    wfs.noise=-1
                    wfs.kernconv=True
                    wfs.sensor_trace(g_dms.dms,1)
                    wfs.comp_image()
                    wfs.noise=tmp_noise
                    wfs.kernconv=False

                self.rtc.do_centroids(ncontro,True)

                h_centroids=self.getCentroids(ncontro)
                control.d_centroids.host2device(<float*>h_centroids.data)

                device=control.d_centroids.getDevice()
                d_centroids=control.d_centroids.getData()
                convert_centro(d_centroids,
                        d_centroids,0,
                        (0.5/dm.push4imat),
                        control.d_centroids.getNbElem(),
                        context.get_device(device))

                carma_axpy[float](context.get_cublasHandle(),
                    control.d_centroids.getNbElem(), <float>-1.0,
                    d_centroids,1,
                    &d_imat.getData()[inds1],1)

                dm.reset_shape()
                inds1+=control.nslope()
            inc(it_dm)


    cdef sensors_compslopes(self, int ncentro, int nmax=-1, float thresh=-1):
        """TODO doc

        :parameters:
            ncentro: (int) : centroider index

            nmax: (int) : (optional)

            thresh: (float) : (optional)
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        self.rtc.d_centro[ncentro].get_cog()


    cdef imat_svd(self,int ncontro):
        """TODO doc

        :param ncontro: controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()
        cdef sutra_controller_ls *controller_ls
        if(type_contro=="ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])

            if(controller_ls.svdec_imat()==1):
                raise RuntimeError("sutra controller has no SVD implementation")

        else:
            raise TypeError("Controller needs to be ls")


    cdef setU(self,int ncontro,np.ndarray[ndim=2,dtype=np.float32_t] U):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            U: (np.ndarray[ndim=2,dtype=np.float32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef sutra_controller_ls *controller_ls
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()

        cdef np.ndarray[dtype=np.float32_t] data_F= U.flatten("F")
        if(type_contro=="ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.d_U.host2device(<float*>data_F.data)


    cdef setEigenvals(self, int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] eigenvals):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            eigenvals: (np.ndarray[ndim=1,dtype=np.float32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef sutra_controller_ls *controller_ls
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.h_eigenvals.fill_from(<float*>eigenvals.data)
        else:
            raise TypeError("Controller needs to be ls")


    cdef getU(self, int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls *controller_ls
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()
        cdef np.ndarray[ndim=2, dtype=np.float32_t] data_F
        cdef np.ndarray[ndim=2, dtype=np.float32_t] data
        cdef const long *dims=NULL

        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            dims=controller_ls.d_U.getDims()
            data_F=np.zeros((dims[2],dims[1]),dtype=np.float32)
            controller_ls.d_U.device2host(<float*>data_F.data)

        data=np.reshape(data_F.flatten("F"),(dims[1],dims[2]))
        return data


    cdef getEigenvals(self,int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls *controller_ls
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()
        cdef np.ndarray[ndim=1, dtype=np.float32_t] data
        cdef const long *dims
        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            dims=controller_ls.h_eigenvals.getDims()
            data=np.zeros((dims[1]),dtype=np.float32)
            controller_ls.h_eigenvals.fill_into(<float*>data.data)

        return data


    cpdef getCenbuff(self, int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls *controller_ls
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()
        cdef np.ndarray[ndim=2, dtype=np.float32_t] data_F
        cdef np.ndarray[ndim=2, dtype=np.float32_t] data
        cdef const long *dims=NULL

        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            dims=controller_ls.d_cenbuff.getDims()
            data_F=np.zeros((dims[2],dims[1]),dtype=np.float32)
            controller_ls.d_cenbuff.device2host(<float*>data_F.data)

        data=np.reshape(data_F.flatten("F"),(dims[1],dims[2]))
        return data


    cdef getErr(self,int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls *controller_ls
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()
        cdef np.ndarray[ndim=1, dtype=np.float32_t] data
        cdef const long *dims
        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            dims=controller_ls.d_err.getDims()
            data=np.zeros((dims[1]),dtype=np.float32)
            controller_ls.d_err.device2host(<float*>data.data)

        return data

    cpdef getCom(self,int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef np.ndarray[ndim=1, dtype=np.float32_t] data
        cdef const long *dims
        dims=self.rtc.d_control[ncontro].d_com.getDims()
        data=np.zeros((dims[1]),dtype=np.float32)
        self.rtc.d_control[ncontro].d_com.device2host(<float*>data.data)

        return data


    cpdef getVoltage(self,int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef np.ndarray[ndim=1, dtype=np.float32_t] data
        cdef const long *dims
        dims=self.rtc.d_control[ncontro].d_voltage.getDims()
        data=np.zeros((dims[1]),dtype=np.float32)
        self.rtc.d_control[ncontro].d_voltage.device2host(<float*>data.data)

        return data


    cpdef getCentroids(self,int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef np.ndarray[ndim=1, dtype=np.float32_t] data
        cdef const long *dims
        dims=self.rtc.d_control[ncontro].d_centroids.getDims()
        data=np.zeros((dims[1]),dtype=np.float32)
        self.rtc.d_control[ncontro].d_centroids.device2host(<float*>data.data)

        IF USE_MPI==1:
            cdef np.ndarray[ndim=1, dtype=np.float32_t] all_centroids
            cdef int comm_size, rank
            mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD,&comm_size)
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD,&rank)

            cdef int *count=<int*>malloc(comm_size*sizeof(int))
            cdef int *disp=<int*>malloc((comm_size+1)*sizeof(int))
            cdef int d=dims[1]/(comm_size*2)
            cdef int i

            disp[0]=0
            for i in range(comm_size):
                if(i<(dims[1]/2)%comm_size):
                    count[i]=d+1
                else:
                    count[i]=d

                disp[i+1]=disp[i]+count[i]

            all_centroids=np.zeros(disp[comm_size]*2,dtype=np.float32)

            cdef float *send=<float*>data.data
            cdef float *recv=<float*>all_centroids.data

            # gather centroids X axis
            mpi.MPI_Allgatherv(send,count[rank],mpi.MPI_FLOAT,
                                recv,count,disp, mpi.MPI_FLOAT,
                                mpi.MPI_COMM_WORLD)

            # gather centroids Y axis
            mpi.MPI_Allgatherv(&send[disp[comm_size]],count[rank],mpi.MPI_FLOAT,
                                &recv[disp[comm_size]],count,disp,
                                mpi.MPI_FLOAT, mpi.MPI_COMM_WORLD)

            free(count)
            free(disp)
            return all_centroids

        ELSE:
            return data


    cdef buildcmat(self,int ncontro,int nfilt, int filt_tt=0):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            nfilt: (int) : 

            filt_tt: (int) : (optional)
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls *controller_ls
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            if(filt_tt>0):
                controller_ls.build_cmat(nfilt,True)
            else:
                controller_ls.build_cmat(nfilt)


    cdef buildcmatmv(self,int ncontro,float cond):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            cond: (float) : 
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            controller_mv.build_cmat(cond)



    cdef loadnoisemat(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] N):
        """TODO doc

        :parameters:
            ncontro: (int) : controller index

            N: (np.ndarray[ndim=1,dtype=np.float32_t]) :
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            controller_mv.load_noisemat(<float*>N.data)



    cpdef docontrol(self,int ncontro):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDevice(self.rtc.device,1)
        self.rtc.do_control(ncontro)


    cpdef applycontrol(self,int ncontro,Dms dms):
        """TODO doc

        :param ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDevice(self.rtc.device,1)
        self.rtc.apply_control(ncontro,dms.dms)




    def __str__(self):
        print "RTC object:"

        cdef sutra_centroider *centro
        cdef sutra_controller *contro
        cdef int i

        info = "Contains "+str(self.rtc.d_centro.size())+" Centroider(s)\n"
        info+= "Centro # | Type  | nwfs | Nvalid\n"

        for i in range(<int>self.rtc.d_centro.size()):
            centro = self.rtc.d_centro[i]
            info+= "%8d"%(i+1)+" | "+"%5s"%centro.get_type()+" | "+"%4d"%(centro.nwfs+1)+\
                    " | "+str(centro.nvalid)+"\n"

        info+= "Contains "+str(self.rtc.d_control.size())+" Controller(s):\n"
        info+= "Control # | Type  | Nslope | Nactu\n"

        for i in range(<int>self.rtc.d_control.size()):
            contro=self.rtc.d_control[i]
            info+= "%9d"%(i+1)+" | "+"%5s"%contro.get_type()+" | "+"%6d"%contro.nslope()+\
                    " | "+str(contro.nactu())+"\n"

        info+= "--------------------------------------------------------"
        return info



def rtc_init(Sensors g_wfs, p_wfs, Dms g_dms, p_dms, Param_geom p_geom, Param_rtc p_rtc,
            Param_atmos p_atmos, Atmos g_atmos, Param_tel p_tel, Param_loop p_loop, 
            Param_target p_tar=None, clean=None, brama=None, doimat=None,simul_name=""):
    """TODO doc

    :parameters:
        g_wfs: (Sensors) :

        p_wfs: (list of Param_wfs) : wfs settings

        g_dms: (Dms) :

        p_dms: (list of Param_dms) : dms settings

        p_geom: (Param_geom) : geom settings

        p_atmos: (Param_atmos) : atmos settings

        g_atmos: (Atmos) :

        p_tel: (Param_tel) : telescope settings

        p_loop: (Param_loop) : loop settings

        p_tar: (Param_target) : (optional) target settings

        clean: (int) : (optional) clean datafiles (imat, U, eigenv)

        brama: (int) : (optional) not implemented yet

        doimat: (int) : (optional) force imat computation

        simul_name: (str) : (optional) simulation's name, use for path to save data (imat, U...)

    """

    g_rtc=Rtc()
    print "TODO brama" 
    """if(brama==1) 
     g_rtc = yoga_rtc_brama(g_rtc);
    """
    if(doimat==None):
        doimat=1

    cdef carma_context *context=carma_context.instance()
    cdef int device = g_rtc.rtc.device
    cdef Param_wfs wfs
    cdef Param_centroider centroider

    cdef int i,j,offset,ncentro,ncontro
    cdef int nwfs =0
    cdef float s_offset=0.
    cdef float s_scale=0.
    cdef Param_controller controller
    cdef np.ndarray[ndim=2,dtype=np.float32_t] imat,corrnorm, KL2V

    cdef int nactu,Nphi
    cdef np.ndarray[ndim=1,dtype=np.float32_t] alt
    cdef char **type_dmseen

    cdef sutra_controller_geo *controller_geo
    cdef np.ndarray[ndim=1,dtype=np.int32_t] indx_pup, indx_dm
    cdef np.ndarray[ndim=1,dtype=np.float32_t] unitpervolt


    cdef np.ndarray[ndim=3,dtype=np.float32_t] tmp,tmp3
    cdef np.ndarray[ndim=2,dtype=np.float32_t] tmp2

    ncentro = len(p_rtc.centroiders)
    ncontro = len(p_rtc.controllers)
    if(p_rtc is not None):
        if(p_wfs is not None):
            for i in range(ncentro):
                centroider=p_rtc.centroiders[i]
                nwfs=centroider.nwfs
                wfs = p_wfs[nwfs]

                if(wfs.type_wfs=="sh"):
                    if( centroider.type_centro!="corr"):
                        s_offset=wfs.npix/2.+0.5
                    else:
                        if(centroider.type_fct=="model"):
                            if(wfs.npix%2==0):
                                s_offset=wfs.npix/2+0.5
                            else:
                                s_offset=wfs.npix/2
                        else:
                            s_offset=wfs.npix/2+0.5

                    s_scale = wfs.pixsize

                elif(wfs.type_wfs=="pyr" or wfs.type_wfs=="roof"):
                    s_offset=0
                    s_scale=0

                g_rtc.add_centroider(g_wfs, nwfs,wfs._nvalid,centroider.type_centro,s_offset,s_scale)
                g_rtc.sensors_initbcube(i)

                if(wfs.type_wfs=="sh"):
                    if(centroider.type_centro=="tcog"):
                        g_rtc.setthresh(i,centroider.thresh)
                    elif(centroider.type_centro=="bpcog"):
                        g_rtc.setnmax(i,centroider.nmax)
                    if(centroider.type_fct=="model"):
                        r0=get_r0(p_atmos.r0, 0.5, wfs.Lambda)
                        seeing= RASC*(wfs.Lambda*1.e-6)/r0
                        npix=seeing/wfs.pixsize
                        if(wfs.gsalt>0):
                            if(wfs.proftype is None or wfs.proftype==""):
                                wfs.proftype="Gauss1"
                            if(wfs.proftype=="Gauss1"):
                                profilename = "allProfileNa_withAltitude_1Gaussian.npy"
                            elif(wfs.proftype=="Gauss2"):
                                profilename = "allProfileNa_withAltitude_2Gaussian.npy"
                            elif(wfs.proftype=="Gauss3"):
                                profilename = "allProfileNa_withAltitude_3Gaussian.npy"
                            elif(wfs.proftype=="Exp"):
                                profilename = "allProfileNa_withAltitude.npy"
                            else:
                                error="Param_wfs proftype unknown: got '"+wfs.proftype+"', expect one of: \n''\n'Gauss1'\n'Gauss2'\n'Gauss3'\n'Exp'"
                                raise ValueError(error)
                            prof=np.load(chakra_ao_savepath+profilename).astype(np.float32)

                            wfs.make_lgs_prof1d(p_tel,np.mean(prof[1:,:],axis=0),prof[0,:],
                                wfs.beamsize,center=<bytes>"image")
                            tmp=wfs._lgskern
                            tmp2=makegaussian(tmp.shape[1],npix*wfs._nrebin).astype(np.float32)
                            tmp3=np.zeros((tmp.shape[1],tmp.shape[1],wfs._nvalid),dtype=np.float32)

                            for j in range(wfs._nvalid):
                                tmp3[:,:,j]= np.fft.ifft2(np.fft.fft2(tmp[:,:,j])*np.fft.fft2(tmp2.T)).real 
                                tmp3[:,:,j]*=tmp3.shape[0]*tmp3.shape[1]

                                tmp3[:,:,j]=np.roll(tmp3[:,:,j],tmp3.shape[0]/2,axis=0)
                                tmp3[:,:,j]=np.roll(tmp3[:,:,j],tmp3.shape[1]/2,axis=1)
                            offset= (wfs._Ntot-wfs._nrebin*wfs.npix)/2
                            j=offset+wfs._nrebin*wfs.npix
                            tmp=np.zeros((j-offset+1,j-offset+1,tmp3.shape[2]),dtype=np.float32)
                            tmp3=np.cumsum(tmp3[offset:j,offset:j,:],axis=0)
                            tmp[1:,1:,:]=np.cumsum(tmp3,axis=1)

                            tmp=np.diff(tmp[::wfs._nrebin,::wfs._nrebin,:],axis=0)
                            tmp=np.diff(tmp,axis=1)

                            centroider.set_weights(tmp)

                        elif(centroider.width==0):
                            centroider.width=npix
                    elif(centroider.type_fct=="gauss"):
                        if(wfs.npix%2==1):
                            centroider.weights = makegaussian(wfs.npix,
                                centroider.width,int(wfs.npix/2+1),int(wfs.npix/2+1)).astype(np.float32)
                        elif(centroider.type_centro=="corr"):
                            centroider.weights = makegaussian(wfs.npix,
                                centroider.width,int(wfs.npix/2),int(wfs.npix/2)).astype(np.float32) 
                        else:
                            centroider.weights = makegaussian(wfs.npix,
                                centroider.width,int(wfs.npix/2)+0.5,int(wfs.npix/2)+0.5).astype(np.float32)

                    if(centroider.weights is None):
                        centroider.weights=np.zeros(0,dtype=np.float32)

                    if(centroider.type_centro=="wcog"):
                        g_rtc.sensors_initweights(i,centroider.weights)
                    elif(centroider.type_centro=="corr"):
                        aa=np.zeros((2*wfs.npix,2*wfs.npix),dtype=np.float32)
                        aa[:wfs.npix,:wfs.npix]=1.0

                        #TODO replace with corrnorm=np.zeros(( dims ))
                        corrnorm=np.fft.ifft2(np.abs(np.fft.fft2(aa))**2).real.astype(np.float32)
                        corrnorm=np.roll(corrnorm*corrnorm.shape[0],corrnorm.shape[0]/2,axis=0)
                        corrnorm=np.roll(corrnorm*corrnorm.shape[1],corrnorm.shape[1]/2,axis=1)[1:,1:]
                        corrnorm=corrnorm*0.+1.

                        centroider.sizex=3
                        centroider.sizey=3
                        centroider.interpmat=create_interp_mat(centroider.sizex,
                            centroider.sizey).astype(np.float32)

                        if(centroider.weights is None):
                            raise ValueError("centroider.weights is None")
                        g_rtc.sensors_initcorr(i,centroider.weights,corrnorm,
                            centroider.sizex,centroider.sizey,centroider.interpmat)


            if(p_wfs is not None and p_dms is not None):
                for i in range(ncontro):
                    controller=p_rtc.controllers[i]
                    print "filtering unseen actuators... "
                    if(p_wfs[0].type_wfs=="sh"):
                        imat=imat_geom(g_wfs,p_wfs,controller,g_dms,p_dms,meth=0)
                    else:
                        imat=manual_imat(g_rtc,g_wfs,p_wfs,g_dms,p_dms)
                    correct_dm(p_dms,g_dms,controller,p_geom, imat,simul_name)
                    #TODO add timer
                    if(controller.type_control !="geo"):
                        nwfs=controller.nwfs
                        if(len(p_wfs)==1):
                            nwfs=p_rtc.controllers[0].nwfs-1
                            # TODO fixing a bug ... still not understood
                        controller.set_nvalid(p_wfs[nwfs]._nvalid)
                    #parameter for add_controller(_geo)
                    ndms=controller.ndm.tolist()
                    controller.set_nactu([p_dms[n-1]._ntotact for n in ndms])
                    nactu=sum([p_dms[j-1]._ntotact for j in ndms])
                    alt=np.array([p_dms[j-1].alt for j in controller.ndm],
                                dtype=np.float32)

                    list_dmseen=[p_dms[j-1].type_dm for j in controller.ndm]

                    type_dmseen=<char**>malloc(controller.ndm.size*sizeof(char*))

                    for j in range(controller.ndm.size):
                        type_dmseen[j]= p_dms[controller.ndm[j]-1].type_dm

                    context.set_activeDeviceForCpy(device,1)
                    if(controller.type_control=="geo"):
                        Nphi=np.where(p_geom._spupil)[0].size
                        list_dmseen,alt,controller.ndm.size
                        g_rtc.rtc.add_controller_geo( nactu, Nphi, controller.delay,
                            device,g_dms.dms, type_dmseen,<float*>alt.data,
                            controller.ndm.size)
                    else:
                        g_rtc.rtc.add_controller( nactu, controller.delay,device,
                            controller.type_control,g_dms.dms, type_dmseen,
                            <float*>alt.data,controller.ndm.size)


                    if(controller.type_control=="geo"):
                        indx_pup = np.where(p_geom._spupil.flatten())[0].astype(np.int32)
                        cpt = 0
                        indx_dm = np.zeros((controller.ndm.size*indx_pup.size),dtype=np.int32)
                        for dmn in range(controller.ndm.size):
                            tmp_s=(p_geom._ipupil.shape[0]-(p_dms[dmn]._n2-p_dms[dmn]._n1+1))/2
                            tmp_e0=p_geom._ipupil.shape[0]-tmp_s
                            tmp_e1=p_geom._ipupil.shape[1]-tmp_s
                            pup_dm=p_geom._ipupil[tmp_s:tmp_e0,tmp_s:tmp_e1]
                            cpt+=np.where(pup_dm)[0].size
                        #convert unitpervolt list to a np.ndarray
                        unitpervolt=np.array([p_dms[j].unitpervolt for j in range(len(p_dms))],
                                    dtype=np.float32)

                        g_rtc.init_proj(i, g_dms, indx_dm, unitpervolt, indx_pup)

                    free(type_dmseen)

                    if(controller.type_control=="ls"):
                        if(doimat):
                            imat_init(i,g_rtc,p_rtc,g_dms,g_wfs,p_wfs,p_tel,clean=1,simul_name=simul_name)
                            if(controller.modopti == 1):
                                print "Initializing Modal Optimization : "
                                if(controller.nrec==0):
                                    controller.nrec=2048
                                else:
                                    #next power of 2 (for fft)
                                    controller.nrec=int(2**np.ceil(np.log2(controller.nrec)))
                                if(controller.nmodes==0):
                                    controller.nmodes=np.sum(p_dms[ndms]._ntotact)
                                if(controller.gmax==0):
                                    controller.gmax=1.0
                                if(controller.ngain==0):
                                    controller.ngain=15
                                KL2V = compute_KL2V(p_dms,controller)
                                g_rtc.init_modalOpti(i,controller.nmodes,controller.nrec,KL2V,
                                    controller.gmin,controller.gmax,controller.ngain,1./p_loop.ittime)
                                ol_slopes=openLoopSlp(g_atmos,g_rtc,controller.nrec,i,g_wfs,p_wfs,p_tar)
                                g_rtc.loadOpenLoop(i,ol_slopes)
                                g_rtc.modalControlOptimization(i)
                            else:
                                cmat_init(i,g_rtc,p_rtc,p_wfs,clean=1,simul_name=simul_name)
                                g_rtc.set_gain(0,controller.gain)
                                mgain=np.ones(sum([p_dms[j]._ntotact for j in range(len(p_dms))]),dtype=np.float32)
                                #filtering tilt ...
                                #mgain(-1:0) = 0.0f;
                                g_rtc.set_mgain(0,mgain)
                        else:
                            nactu=np.sum(controller.nactu)
                            nvalid=np.sum(controller.nvalid)
                            imat=np.zeros((nactu*nvalid*2),dtype=np.float32)
                            g_rtc.set_imat(i,imat)
                            g_rtc.set_cmat(i,imat)
                            g_rtc.set_gain(0,controller.gain)
                            mgain=np.ones(sum([p_dms[j]._ntotact for j in range(len(p_dms))]),dtype=np.float32)
                            # filtering tilt ...
                            g_rtc.set_mgain(0,mgain)

                    if(controller.type_control=="cured"):
                        print "initializing cured controller"
                        if( "tt" in [p_dms[j].type_dm for j in range(len(p_dms))]):
                            tt_flag=1
                        else:
                            tt_flag=0
                        print "TODO controller_initcured"
                        #controller_initcured,g_rtc,0,int(y_wfs(1).nxsub),int(*y_wfs(1)._isvalid),
                                #int(controllers(i).cured_ndivs),int(tt_flag);
                        g_rtc.set_gain(0,controller.gain)
                    if(controller.type_control=="kalman_CPU" or 
                       controller.type_control=="kalman_GPU"):
                        env_var=os.environ.get("COMPILATION_LAM")
                        #TODO found_str = strfind("standalone",env_var);
                        #TODO if (found_str(2) != -1)
                        print "\nWARNING : Environment variable COMPILATION_LAM contains the word \"standalone\". Make sure that this variable did not contain \"standalone\" when compiling, which would mean that Kalman filter was compiled for standalone version (in lam/kalman_CPU_GPU/test), which is not compatible with COMPASS\n."
                        if(controller.type_control=="kalman_GPU"):
                            print "initializing kalman_GPU controller"
                        else:
                            print "initializing kalman_CPU controller"

                        g_rtc.set_gain(0,controller.gain)

                        print "TODO :rtc_init kalman case"
                        #TODO D_Mo = create_dmo(1,1);

                        # creation de N_Act (en um/V et non normalise)
                        #TODO N_Act  = create_nact(1); //N_Act = N_Act-min(N_Act); N_Act = N_Act/max(N_Act);
                        #TODO PROJ   = LUsolve(N_Act); inverse N_act

                        #Creation de SigmaTur puis conversion de rad^2 en um^2
                        #TODO SigmaTur = create_sigmaTur(1) * (y_wfs(nwfs).lambda/2/pi)^2;

                        #SigmaTur = create_sigmaTur(1);//pli,SigmaTur;
                        #TODO atur = (0.985, sum(y_dm(ndms)._ntotact));
                        #TODO btur = array(0.0f, sum(y_dm(ndms)._ntotact));
                        #TODO ordreAR = anyof(btur)+1;
                        #TODO isZonal=1 ; isSparse=1;
                        #TODO SigmaV = create_sigmav(SigmaTur, isZonal, ordreAR, atur, btur) ;

                        
               
                        if (controller.type_control  == "kalman_CPU"):
                            print "TODO rtc_initkalman, g_rtc, 0, avg(noise_cov(1)), D_Mo, N_Act, PROJ, SigmaV, atur, btur, isZonal, isSparse, 0;"
                        else:
                            print  "TODO rtc_initkalman, g_rtc, 0, avg(noise_cov(1)), D_Mo, N_Act, PROJ, SigmaV, atur, btur, isZonal, isSparse, 1;"

                    elif(controller.type_control=="generic"):
                        size=sum([p_dms[j]._ntotact for j in range(len(p_dms))])
                        decayFactor=np.zeros(size ,dtype=np.float32)
                        mgain=np.zeros(size ,dtype=np.float32)
                        matE=np.zeros((size,size) ,dtype=np.float32)
                        cmat=np.zeros((size,np.sum(controller.nvalid)*2),
                                       dtype=np.float32)
                        #TODO ? law = "integrator";

                        g_rtc.set_decayFactor(0,decayFactor)
                        g_rtc.set_mgain(0,mgain)
                        g_rtc.set_cmat(0,cmat)
                        g_rtc.set_matE(0,matE)
            
    return g_rtc
"""
          // Florian features
          if (controllers(i).type == "mv"){   
            write,format="%s", "doing imat_geom... ";
            tic;
            //imat_init,i,clean=clean;
            imat = imat_geom(i,meth=0);
            rtc_setimat,g_rtc,i-1,imat;
            write,format = "done in : %f s\n",tac();  
            rtc_setgain,g_rtc,0,controllers(i).gain;
            mgain = array(1.0f,(y_dm._ntotact)(sum));
            rtc_loadmgain,g_rtc,0,mgain;
        
            Cphim = mat_cphim_gpu(i-1);

            cmat_init,i,clean=clean;
           
          }
"""






cdef correct_dm(p_dms, Dms g_dms, Param_controller p_control, Param_geom p_geom, np.ndarray imat, bytes simul_name):
    """TODO doc

    :parameters:
        p_dms: (list of Param_dm) : dms settings

        g_dms: (Dms) : 

        p_control: (Param_controller) : controller settings

        p_geom: (Param_geom) : geom settings

        imat: (np.ndarray) : interaction matrix

        simul_name: (str) : simulation's name, use for data files' path
    """
    #cdef carma_context *context = carma_context.instance() #UNUSED
    
    cdef int i, nm, nmc, inds,ndm,nactu_nm
    cdef np.ndarray[ndim=1,dtype=np.float32_t] resp

    cdef bytes filename
    cdef bytes dirsave = chakra_ao_savepath+"mat/"

    cdef long dims,ninflu,influsize,NR,NP

    if (simul_name=="" ): 
        imat_clean = 1
    else:
        imat_clean=0

    ndm=p_control.ndm.size
    for i in range(ndm):
        nm=p_control.ndm[i]-1
        g_dms.remove_dm(p_dms[nm].type_dm,p_dms[nm].alt)

    resp=np.sqrt(np.sum(imat**2,axis=0))

    inds=0

    cdef int rank
    IF USE_MPI==1:
        mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD,&rank)
    ELSE:
        rank=0

    for nmc in range(ndm):
        nm=p_control.ndm[nmc]-1
        nactu_nm=p_dms[nm]._ntotact
        # filter actuators only in stackarray mirrors:

        if(p_dms[nm].type_dm=="pzt"):
            filename=dirsave+"pztok-"+str(nm)+"-"+simul_name+".npy"

            if(imat_clean<1 and os.path.isfile(filename)):
                ok=np.load(filename)
                filename=dirsave+"pztnok-"+str(nm)+"-"+simul_name+".npy"
                nok=np.load(filename)
            else:
                tmp=resp[inds:inds+p_dms[nm]._ntotact]
                ok=np.where(tmp.flatten("F")>p_dms[nm].thresh*np.max(tmp))[0]
                nok=np.where(tmp.flatten("F")<=p_dms[nm].thresh*np.max(tmp))[0]
                if(simul_name !="" and rank==0):
                    np.save(filename,ok)
                    filename=dirsave+"pztnok-"+str(nm)+"-"+simul_name+".npy"
                    np.save(filename,nok)

            p_dms[nm].set_xpos(p_dms[nm]._xpos[ok])
            p_dms[nm].set_ypos(p_dms[nm]._ypos[ok])
            p_dms[nm].set_i1(p_dms[nm]._i1[ok])
            p_dms[nm].set_j1(p_dms[nm]._j1[ok])
            p_dms[nm].set_influ(p_dms[nm]._influ[:,:,ok.tolist()])
            p_dms[nm].set_ntotact(p_dms[nm]._influ.shape[2])

            comp_dmgeom(p_dms[nm],p_geom)

            ninflu=long(p_dms[nm]._ntotact)
            influsize=long(p_dms[nm]._influsize)
            ninflupos=long(p_dms[nm]._influpos.size)
            n_npts=long(p_dms[nm]._ninflu.size)

            dims=long(p_dms[nm]._n2-p_dms[nm]._n1+1)
            dims=max(dims,p_geom._mpupil.shape[1])
            g_dms.add_dm(<bytes>"pzt",p_dms[nm].alt, dims, ninflu, influsize,
                            ninflupos, n_npts,p_dms[nm].push4imat)
            g_dms.load_pzt(p_dms[nm].alt,p_dms[nm]._influ,p_dms[nm]._influpos.astype(np.int32),
                p_dms[nm]._ninflu, p_dms[nm]._influstart, p_dms[nm]._i1,p_dms[nm]._j1,p_dms[nm]._influkernel)

        elif(p_dms[nm].type_dm=="tt"):
            dim=long(p_dms[nm]._n2-p_dms[nm]._n1+1)
            g_dms.add_dm(<bytes>"tt",p_dms[nm].alt,dim,2,dim,1,1,p_dms[nm].push4imat)
            g_dms.load_tt(p_dms[nm].alt,p_dms[nm]._influ)


        elif(p_dms[nm].type_dm=="kl"):
            dim=long(p_dms[nm]._n2-p_dms[nm]._n1+1)
            ninflu=long(p_dms[nm].nkl)
            influsize=long(p_dms[nm]._klbas.ncp)
            _nr=long(p_dms[nm]._klbas.nr)
            _np=long(p_dms[nm]._klbas.np)
            g_dms.add_dm(<bytes>"kl",p_dms[nm].alt, dim, ninflu, influsize,
                            _nr, _np,p_dms[nm].push4imat)
            g_dms.load_kl(p_dms[nm].alt,p_dms[nm]._klbas.rabas,p_dms[nm]._klbasazbas,
                            p_dms[nm]._klbas.ord,p_dms[nm]._klbas.cr,p_dms[nm]._klbas.cp)

        inds += nactu_nm



cdef imat_geom(Sensors g_wfs, p_wfs, Param_controller p_control,Dms g_dms, p_dms, int meth=0):
    """TODO doc

    :parameters:
        g_wfs: (Sensors) :

        p_wfs: (list of Param_wfs) : wfs settings

        p_control: (Param_controller) : controller settings

        g_dms: (Dms) :

        p_dms: (list of Param_dm) : dms settings

        meth: (int) : (optional)
    """

    cdef int nwfs=p_control.nwfs.size
    cdef int ndm =p_control.ndm.size
    cdef int nslps=0
    cdef int nmc,nm,nw,i,ind
    cdef int wfs
    cdef int imat_size1=0
    cdef int imat_size2=0 

    cdef np.ndarray[ndim=2,dtype=np.float32_t] imat_cpu

    for nw in range(nwfs):
        nm=p_control.ndm[nw]-1
        imat_size1+=p_wfs[nw]._nvalid*2

    for nmc in range(ndm):
        imat_size2 +=p_dms[nmc]._ntotact

    imat_cpu = np.zeros((imat_size1,imat_size2),dtype=np.float32)
    ind=0

    for nmc in range(ndm):
        nm=p_control.ndm[nmc]-1
        g_dms.resetdm(p_dms[nm].type_dm,p_dms[nm].alt)
        for i in range(p_dms[nm]._ntotact):
            g_dms.oneactu(p_dms[nm].type_dm,p_dms[nm].alt,i,p_dms[nm].push4imat)
            nslps=0
            for nw in range(nwfs):
                wfs=p_control.nwfs[nw]-1
                g_wfs.sensors_trace(wfs,"dm",None,dms=g_dms,rst=1)
                g_wfs.slopes_geom(wfs,meth)
                imat_cpu[nslps:nslps+p_wfs[wfs]._nvalid*2,ind]=g_wfs._get_slopes(wfs)
                nslps+=p_wfs[wfs]._nvalid*2
            imat_cpu[:,ind]=imat_cpu[:,ind]/float(p_dms[nm].push4imat)
            ind=ind+1
            g_dms.resetdm(p_dms[nm].type_dm,p_dms[nm].alt)
    return imat_cpu


cdef manual_imat(Rtc g_rtc,Sensors g_wfs, p_wfs, Dms g_dms, p_dms):
    """TODO doc

    :parameters:
        g_rtc: (Rtc) :

        g_wfs: (Sensors) :

        g_dms: (Dms) :

        p_dms: (list of Param_dm) : dm settings
    """

    cdef int nm,i,ind

    #cdef np.ndarray[ndim=1, dtype=np.float32_t] slps=g_rtc.getcentroids(0, g_wfs, 0) #UNUSED
    cdef int nslps = g_wfs._get_slopesDims(0)
    cdef int imat_size2=0

    cdef np.ndarray[ndim=2,dtype=np.float32_t] imat_cpu
    cdef np.ndarray[ndim=1,dtype=np.float32_t] com

    for nm in range(len(p_dms)):
        g_dms.resetdm(p_dms[nm].type_dm,p_dms[nm].alt)
        imat_size2+=p_dms[nm]._ntotact

    imat_cpu=np.zeros((nslps,imat_size2),dtype=np.float32)    


    ind=0

    for nm in range(len(p_dms)):
        for i in range(p_dms[nm]._ntotact):
            com=np.zeros((p_dms[nm]._ntotact),dtype=np.float32)
            com[i]=float(p_dms[nm].push4imat)
            g_dms.set_comm(p_dms[nm].type_dm,p_dms[nm].alt,com)
            g_dms.shape_dm(p_dms[nm].type_dm,p_dms[nm].alt)
            g_wfs.sensors_trace(0,"dm", None,dms=g_dms,rst=1)
            #equivalent to Bcast(g_wfs.sensors.d_wfs[0].d_gs.d_phase.d_screen)
            #g_wfs.Bcast_dscreen()
            IF USE_MPI==1:
                Bcast(g_wfs.sensors.d_wfs[0].d_gs.d_phase.d_screen,0)
            g_wfs.sensors_compimg(0)
            g_rtc.docentroids()


            if(p_wfs[0].type_wfs!="pyr"):
                imat_cpu[:,ind]=g_rtc.getcentroids(0,g_wfs,0)/float(p_dms[0].push4imat)
            else:
                imat_cpu[:,ind]=g_wfs._get_slopes(0)/float(p_dms[0].push4imat)
            g_dms.resetdm(p_dms[nm].type_dm,p_dms[nm].alt)
            ind+=1
   
    return imat_cpu


cdef get_r0(float r0_at_lambda1, float lambda1, float lambda2):
    """Compute r0

    :parameters:
        r0_at_lambda1: (float) :

        lambda1: (float) :

        lambda2: (float) :
    """
    return (lambda2/lambda1)**(6./5.)*r0_at_lambda1



def create_interp_mat(int dimx, int dimy):
    """TODO doc

    :parameters:
        dimx: (int) :

        dimy: (int) :
    """
    n=max(dimx,dimy)
    tmp1,tmp2=indices(n)
    tmp1=tmp1[:dimx,:dimy]-(dimx/2+1)
    tmp2=tmp2[:dimx,:dimy]-(dimy/2+1)

    tmp=np.zeros((tmp1.size,6),np.int32)
    tmp[:,0]=(tmp1**2).flatten()
    tmp[:,1]=(tmp2**2).flatten()
    tmp[:,2]=(tmp1*tmp2).flatten()
    tmp[:,3]=tmp1.flatten()
    tmp[:,4]=tmp2.flatten()
    tmp[:,5]=1

    return np.dot(np.linalg.inv(np.dot(tmp.T,tmp)),tmp.T).T



    

cdef compute_KL2V(p_dms, Param_controller controller):
    """TODO doc

    :parameters:
        p_dms: (list of Param_dm) : dms settings

        controller: (Param_controller) : controller settings
    """
    cdef int i,nTT,indx_act, ndm
    cdef np.ndarray[ndim=1,dtype=np.int64_t] ntotact=\
        np.array([p_dms[i]._ntotact for i in range(len(p_dms))],dtype=np.int64)
    cdef np.ndarray[ndim=2,dtype=np.float32_t] KL2V=\
        np.zeros((np.sum(ntotact),np.sum(ntotact)),dtype=np.float32)

    indx_act=0
    nTT=0

    for i in range(controller.ndm.size):
        ndm=controller.ndm[i]
        if(p_dms[ndm].type_dm=="pzt"):
            print "TODO KL2V"
            #KL2V[indx_act:indx_act+ntotact[ndm],indx_act:indx_act+ntotact[ndm]]=\
                #TODO compute_klbasis(ndm)
        elif(p_dms[ndm].type_dm=="tt"):
            nTT+=1

    KL2V=KL2V[:,controller.nmodes]
    if(nTT!=0):
        KL2V[:,:controller.nmodes-2]= KL2V[:,2:]
        KL2V[:,controller.nmodes-2:]= np.zeros((np.sum(ntotact,2)),dtype=np.float32)
        KL2V[np.sum(ntotact)-2:,controller.nmodes-2:]=np.identity(2,dtype=np.float32)

    return KL2V

cdef openLoopSlp(Atmos g_atm, Rtc g_rtc,int nrec, int ncontro, Sensors g_wfs=None,  
        p_wfs=None, Param_target p_tar=None,Target g_tar=None):
    """TODO doc

    :parameters:
        g_atm: (Atmos) :

        g_rtc: (Rtc) :

        nrec: (int) :

        ncontro: (int) :

        g_wfs: (Sensors) :

        p_wfs: (list of Param_wfs) : (optional) wfs settings

        p_tar: (Param_target) : (optional) target settings

        g_tar: (Target) : (optional)
    """
    #TEST IT
    cdef int i,j
    cdef np.ndarray[ndim=2,dtype=np.float32_t] ol_slopes=\
        np.zeros((sum([p_wfs[i]._nvalid for i in range(len(p_wfs))]),nrec),
        dtype=np.float32)

    for i in range(nrec):
        print "Reconring"+str(nrec)+"open-loop slpoes :"+str(int(i/nrec)*100)
        if(g_tar is not None):
            g_atm.move_atmos()
            if(p_tar is not None):
                for j in range(p_tar.ntarget):
                    g_tar.atmos_trace(j,g_atm)

        if(p_wfs is not None and g_wfs is not None):
            for j in range(len(p_wfs)):
                g_wfs.sensors_trace(j,"atmos",g_atm)
                g_wfs.sensors_compimg(j)
                g_rtc.sensors_compslopes(ncontro)
                ol_slopes[j*p_wfs[j]._nvalid*2:j*p_wfs[j]._nvalid*2,i]=g_wfs.get_slopes(j)

    return ol_slopes

cdef imat_init(int ncontro, Rtc g_rtc, Param_rtc p_rtc, Dms g_dms, Sensors g_wfs,
        p_wfs, Param_tel p_tel, int clean=1, bytes simul_name=<bytes>""):
    """TODO doc

    :parameters:
        ncontro: (int) :

        g_rtc: (Rtc) :

        p_rtc: (Param_rtc) : rtc settings

        g_dms: (Dms) :

        g_wfs: (Sensors) :

        p_wfs: (list of Param_wfs) : wfs settings

        p_tel: (Param_tel) : telescope settings

        clean: (int) : (optional) : clean datafiles (imat, U, eigenv)

        simul_name: (str) : (optional) simulation's name, use for data files' path
    """
    cdef bytes dirsave=chakra_ao_savepath+"mat/"
    cdef bytes filename=dirsave+"imat-"+str(ncontro)+"-"+simul_name+".npy"
    cdef bytes profilename=chakra_ao_savepath+<bytes>"allProfileNa_withAltitude_1Gaussian.npy"
    cdef int imat_clean=1
    cdef int i
    cdef double t0

    cdef int rank
    cdef int world

    IF USE_MPI==1:
        mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD,&rank)
        mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD,&world)
    ELSE:
        rank=0
        world=1

    cdef sutra_wfs *wfs
    cdef carma_obj[float] *screen


    if(simul_name!=""):
        imat_clean=int(not os.path.isfile(filename) or clean)

    if(imat_clean):
        # first check if wfs is using lgs
        # if so, load new lgs spot, just for imat
        for i in range(len(p_wfs)):
            if(p_wfs[i].gsalt>0):
                prof=np.load(profilename)
                h=prof[0,:]
                prof=prof[1:,:]
                prof=np.mean(prof,axis=0)
                p_wfs[i].prep_lgs_prof(i,p_tel,prof,h,
                                        p_wfs[i].beamsize,g_wfs,imat=1)

        print "doing imat..."
        t0=time.time()
        g_rtc.doimat(ncontro,g_dms)
        print "done in ",time.time()-t0
        p_rtc.controllers[ncontro].set_imat(g_rtc.get_imat(ncontro))
        if(simul_name!="" and rank==0):
            np.save(filename,p_rtc.controllers[ncontro].imat)

    else:
        p_rtc.controllers[ncontro].set_imat(np.load(filename))
        g_rtc.set_imat(ncontro, p_rtc.controllers[ncontro].imat)

    #now restore original profile in lgs spots
    for i in range(len(p_wfs)):
        if(p_wfs[i].gsalt>0):
            p_wfs[i].prep_lgs_prof(i,p_tel,p_wfs[i]._profna,p_wfs[i]._altna,
                            p_wfs[i].beamsize,g_wfs)




cdef cmat_init(int ncontro, Rtc g_rtc, Param_rtc p_rtc, list p_wfs,
                clean=1, bytes simul_name=<bytes>""):
    """TODO doc

    :parameters:
        ncontro: (int) :
        g_rtc: (Rtc) :
        p_rtc: (Param_rtc) : rtc settings
        p_wfs: (list of Param_wfs) : wfs settings
        clean: (int) : (optional) clean datafiles (imat, U, eigenv)
        simul_name: (str) : (optional) simulation's name, use for data files' path
    """

    cdef bytes dirsave=chakra_ao_savepath+"mat/"
    cdef bytes filename

    cdef int cmat_clean
    cdef double t0
    cdef np.ndarray[ndim=1,dtype=np.float32_t] eigenv, N
    cdef np.ndarray[ndim=2,dtype=np.float32_t] U, imat
    cdef np.ndarray[ndim=1,dtype=np.int64_t] mfilt

    cdef float maxcond
    cdef int nfilt,ind,i


    if(simul_name==""):
        cmat_clean=1
    else:
        filename=dirsave+"imat-"+str(ncontro)+"-"+simul_name+".npy"
        cmat_clean=int(os.path.isfile(filename) or clean)

    cdef int rank
    IF USE_MPI:
        mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD,&rank)
    ELSE:
        rank=0

    if(p_rtc.controllers[ncontro].type_control == "ls"):
        if(cmat_clean):
            print "doing svd"
            t0=time.time()
            g_rtc.imat_svd(ncontro)
            print "svd time",time.time()-t0
            eigenv = g_rtc.getEigenvals(ncontro)
            if(simul_name!="" and rank==0):
                U=g_rtc.getU(ncontro)
                filename=dirsave+"eigenv-"+str(ncontro)+"-"+simul_name
                np.save(filename,eigenv)
                filename=dirsave+"U-"+str(ncontro)+"-"+simul_name
                np.save(filename,U)
        else:
            filename=dirsave+"eigenv-"+str(ncontro)+"-"+simul_name+".npy"
            eigenv=np.load(filename)
            filename=dirsave+"U-"+str(ncontro)+"-"+simul_name+".npy"
            U=np.load(filename)
            g_rtc.seteigenvals(ncontro,eigenv)
            g_rtc.setU(ncontro,U)

        imat = g_rtc.get_imat(ncontro)
        maxcond=p_rtc.controllers[ncontro].maxcond
        if(eigenv[0]<eigenv[eigenv.shape[0]-1]):
            mfilt=np.where((eigenv/eigenv[eigenv.shape[0]-3]) < 1./maxcond)[0]
        else:
            mfilt=np.where( (1./(eigenv/eigenv[2]))>maxcond)[0]
        nfilt=mfilt.shape[0]

        #print "TODO wfs_disp"
        """
        if ( (wfs_disp!=[]) && (numberof(*wfs_disp._winits) > 0)) {
            if ((*wfs_disp._winits)(5)) {
                window,(*wfs_disp._wins)(5);fma;logxy,0,1;
                if (eigenv(1) < eigenv(0)) {
                    plg, eigenv(::-1), marks=0;
                    plmk, eigenv(::-1), msize = 0.3, marker=4;
                } else {
                    plg, eigenv, marks=0;
                    plmk, eigenv, msize = 0.3, marker=4;
                }
            x0 = dimsof(imat)(3) - nfilt + 0.5;
            pldj, x0 ,min(eigenv), x0, max(eigenv), color="red";
            }
        }

        """

        print "building cmat"
        print "filtering ",nfilt," modes"
        t0=time.time()
        g_rtc.buildcmat(ncontro,nfilt)
        print "cmat time ",time.time()-t0

    if(p_rtc.controllers[ncontro].type_control=="mv"):
        N=np.zeros((2*np.sum(p_wfs[p_rtc.controllers[ncontro].nwfs]._nvalid)),dtype=np.int32)
        ind=0
        for i in range(p_rtc.controllers[ncontro].nwfs.size):
            k=p_rtc.controllers[ncontro].nwfs[i]
            #N[]=noise_cov(k)
            ind+=2*p_wfs[k]._nvalid

        g_rtc.loadnoisemat(ncontro,N)
        print "Building cmat..."
        g_rtc.build_cmatmv(ncontro,p_rtc.controllers[ncontro].maxcond)

    '''
    TODO what is y_controllers(ncontrol-1).TTcond

    if (((*y_rtc.controllers(ncontrol)).type)(1) == "mv"){
        rtc_buildcmatmv,g_rtc,ncontrol-1,y_controllers(ncontrol-1).maxcond;
        if(y_controllers(ncontrol-1).TTcond == 0) 
            y_controllers(ncontrol-1).TTcond = y_controllers(ncontrol-1).maxcond;
        if(anyof(y_dm.type == "tt"))
            rtc_filtercmatmv,g_rtc,ncontrol-1,y_controllers(ncontrol-1).TTcond;
    }
    '''
    p_rtc.controllers[ncontro].set_cmat(g_rtc.get_cmat(ncontro))


