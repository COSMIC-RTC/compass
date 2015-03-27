cdef class Sensors:
    def __cinit__(self, int nsensors, bytes type_data,
                    np.ndarray[dtype=np.int64_t] npup,
                    np.ndarray[dtype=np.int64_t] nxsub,
                    np.ndarray[dtype=np.int64_t] nvalid,
                    np.ndarray[dtype=np.int64_t] nphase,
                    np.ndarray[dtype=np.float32_t] pdiam,
                    np.ndarray[dtype=np.int64_t] npix=None,
                    np.ndarray[dtype=np.int64_t] nrebin=None,
                    np.ndarray[dtype=np.int64_t] nfft=None,
                    np.ndarray[dtype=np.int64_t] ntota=None,
                    np.ndarray[dtype=np.float32_t] nphot=None,
                    np.ndarray[dtype=np.int32_t] lgs=None,
                    int odevice=-1,
                    int comm_size=1,
                    int rank=0
                ):

        cdef char *type
        cdef carma_context *context= carma_context.instance()
        if odevice <0:
            odevice=context.get_activeDevice()

        if(type_data=="geo"):
            self.sensors= new sutra_sensors(context, nsensors,
                        <long*>nxsub.data,
                        <long*>nvalid.data,
                        <long*>nphase.data,
                        npup[0],
                        <float*>pdiam.data,
                        odevice)
        else:
            type=<char*>type_data
            self.sensors= new sutra_sensors(context, &type, nsensors,
                        <long*>nxsub.data,
                        <long*>nvalid.data,
                        <long*>npix.data,
                        <long*>nphase.data,
                        <long*>nrebin.data,
                        <long*>nfft.data,
                        <long*>ntota.data,
                        <long*>npup.data,
                        <float*>pdiam.data,
                        <float*>nphot.data,
                        <int*>lgs.data,
                        odevice)
        self.sensors.define_mpi_rank(rank,comm_size)
        self.sensors.allocate_buffers()
        self.sensors.device=odevice


    cdef sensors_initgs(self,np.ndarray[dtype=np.float32_t] xpos,
                             np.ndarray[dtype=np.float32_t] ypos,
                             np.ndarray[dtype=np.float32_t] Lambda,
                             np.ndarray[dtype=np.float32_t] mag,
                             np.ndarray[dtype=np.int64_t  ] size,
                             np.ndarray[dtype=np.float32_t] noise,
                             np.ndarray[dtype=np.int64_t  ] seed):
       
        if(noise.size==0):
            self.sensors.sensors_initgs(<float*>xpos.data, <float*>ypos.data,
                        <float*>Lambda.data, <float*>mag.data, <long*>size.data)

        elif(seed.size ==0):
            self.sensors.sensors_initgs(<float*>xpos.data, <float*>ypos.data,
                        <float*>Lambda.data, <float*>mag.data, <long*>size.data,
                        <float*>noise.data)
        else:
            self.sensors.sensors_initgs(<float*>xpos.data, <float*>ypos.data,
                        <float*>Lambda.data, <float*>mag.data, <long*>size.data,
                        <float*>noise.data, <long*>seed.data)


#            int *istart, int *jstart, float *ftkern_sinc, int *hrmap, int *binmap,
#            float *halfxy):

    cdef sensors_initarr(self,int n, Param_wfs wfs, Param_geom geom):

        cdef np.ndarray tmp,tmp2,tmp_istart,tmp_jstart
        cdef string type_sh="sh"
        cdef string type_pyr="pyr"
        cdef string type_roof="roof"
        cdef string type_geo="geo"

        cdef sutra_wfs_geom *wfs_geom=NULL
        cdef sutra_wfs_sh *wfs_sh = NULL
        cdef sutra_wfs_pyr_roof *wfs_roof = NULL
        cdef sutra_wfs_pyr_pyr4 *wfs_pyr = NULL

        cdef int* phasemap=<int*>wfs._phasemap.data
        cdef int* hrmap=<int*>wfs._hrmap.data
        cdef int* binmap=<int*>wfs._binmap.data
        cdef int* validx=<int*>wfs._validsubsx.data
        cdef int* validy=<int*>wfs._validsubsy.data
        tmp_istart=np.copy(wfs._istart+1)
        cdef int* istart=<int*>tmp_istart.data
        tmp_jstart=np.copy(wfs._jstart+1)
        cdef int* jstart=<int*>tmp_jstart.data
        cdef int* cx = <int*>wfs._pyr_cx.data
        cdef int* cy = <int*>wfs._pyr_cy.data

        cdef float* halfxy=<float*>wfs._halfxy.data
        cdef float* pupil=<float*>geom._mpupil.data
        cdef float* ftkernel=<float*>wfs._ftkernel.data
        cdef float* submask=<float*>wfs._submask.data
        cdef float* sincar=<float*>wfs._sincar.data


        tmp=wfs._fluxPerSub.astype(np.float32)
        tmp2=wfs._isvalid.astype(np.int32)
        tmp=tmp[np.where(tmp2>0)]
        cdef float* fluxPerSub=<float*>tmp.data


        if(self.sensors.d_wfs[n].type==type_geo):
            wfs_geom=dynamic_cast_wfs_geom_ptr(self.sensors.d_wfs[n])
            # TODO tmp
            #wfs_geom.wfs_initarrays(phasemap,offset, pupil,flux_foc,isvalid,
            #    validsubsx, validsubsy)
            #tmp.data=wfs._isvalid
            #tmp2=wfs._fluxPerSub.flatten()[tmp]#wfs._isvalid]
            #wfs_geom.wfs_initarrays(wfs._phasemap, wfs._halfxy, pupil, <float*>tmp.data,
            #        wfs._isvalid, wfs._validsubsx,wfs._validsubsy)
            """
  sensors_initarr,g_wfs,i-1,int(*y_wfs(i)._phasemap),float(*y_wfs(i)._halfxy),float(*y_geom._mpupil),
    (*y_wfs(i)._fluxPerSub)(where(*y_wfs(i)._isvalid)),int(*y_wfs(i)._isvalid),
    int((*y_wfs(i)._validsubs)(1,)-1),int((*y_wfs(i)._validsubs)(2,)-1),int(*y_wfs(i)._istart+1),
    int(*y_wfs(i)._jstart+1);

            """


        elif(self.sensors.d_wfs[n].type==type_pyr):
            wfs_pyr = dynamic_cast_wfs_pyr_pyr4_ptr(self.sensors.d_wfs[n])
            """ TODO init tmp arrays
  tmp = array(float,2,y_wfs(i)._Ntot,y_wfs(i)._Ntot);
  tmp(1,,) = (*y_wfs(i)._halfxy).re;
  tmp(2,,) = (*y_wfs(i)._halfxy).im;
  
  tmp2 = array(float,2,y_wfs(i)._Nfft,y_wfs(i)._Nfft);
  tmp2(1,,) = (*y_wfs(i)._pyr_offsets).re;
  tmp2(2,,) = (*y_wfs(i)._pyr_offsets).im;

  sensors_initarr,g_wfs,i-1,float(tmp),float(tmp2),float(*y_wfs(i)._submask),float(*y_geom._mpupil),
    int(*y_wfs(i)._isvalid),int(*y_wfs(i)._pyr_cx),int(*y_wfs(i)._pyr_cy),float(*y_wfs(i)._hrmap),
    int(*y_wfs(i)._phasemap),int((*y_wfs(i)._validsubs)(1,)-1),int((*y_wfs(i)._validsubs)(2,)-1);
            """
            tmp=wfs._halfxy.astype(np.complex64)
            tmp2=wfs._pyr_offsets.astype(np.complex64)
            wfs_pyr.wfs_initarrays(<cuFloatComplex*>tmp.data,<cuFloatComplex*>tmp2.data,
                    submask, pupil,cx, cy, 
                    sincar, phasemap, validx,validy)

        elif(self.sensors.d_wfs[n].type==type_roof):
            wfs_roof = dynamic_cast_wfs_pyr_roof_ptr(self.sensors.d_wfs[n])
            """ TODO init tmp arrays
  tmp = array(float,2,y_wfs(i)._Ntot,y_wfs(i)._Ntot);
  tmp(1,,) = (*y_wfs(i)._halfxy).re;
  tmp(2,,) = (*y_wfs(i)._halfxy).im;
  
  tmp2 = array(float,2,y_wfs(i)._Nfft,y_wfs(i)._Nfft);
  tmp2(1,,) = (*y_wfs(i)._pyr_offsets).re;
  tmp2(2,,) = (*y_wfs(i)._pyr_offsets).im;

  sensors_initarr,g_wfs,i-1,float(tmp),float(tmp2),float(*y_wfs(i)._submask),float(*y_geom._mpupil),
    int(*y_wfs(i)._isvalid),int(*y_wfs(i)._pyr_cx),int(*y_wfs(i)._pyr_cy),float(*y_wfs(i)._hrmap),
    int(*y_wfs(i)._phasemap),int((*y_wfs(i)._validsubs)(1,)-1),int((*y_wfs(i)._validsubs)(2,)-1);
            """
            tmp=wfs._halfxy.astype(np.complex64)
            tmp2=wfs._pyr_offsets.astype(np.complex64)
            wfs_roof.wfs_initarrays(<cuFloatComplex*>tmp.data,<cuFloatComplex*>tmp2.data,
                    submask, pupil,cx, cy, 
                    sincar, phasemap, validx,validy)


        elif(self.sensors.d_wfs[n].type==type_sh):
            wfs_sh=dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[n])
            wfs_sh.wfs_initarrays(phasemap,hrmap,binmap,halfxy,
                    pupil,fluxPerSub,validy, validx,
                    istart, jstart,<cuFloatComplex*>ftkernel)
            

            """
            sensors_initarr,g_wfs,i-1,int(*y_wfs(i)._phasemap),int(*y_wfs(i)._hrmap),
    int(*y_wfs(i)._binmap),float(*y_wfs(i)._halfxy),float(*y_geom._mpupil),
    (*y_wfs(i)._fluxPerSub)(where(*y_wfs(i)._isvalid)),int(*y_wfs(i)._isvalid),
    int((*y_wfs(i)._validsubs)(1,)-1),int((*y_wfs(i)._validsubs)(2,)-1),int(*y_wfs(i)._istart+1),
    int(*y_wfs(i)._jstart+1),float(*y_wfs(i)._ftkernel);

            """

    cdef sensors_addlayer(self,int i, bytes type, float alt,
        float xoff, float yoff):
       
       cdef carma_context * context = carma_context.instance()
       context.set_activeDevice(self.sensors.device,1)

       self.sensors.d_wfs[i].d_gs.add_layer(<char*>type, alt, xoff, yoff)
         

    def sensors_compimg(self, int n):
        cdef carma_context *context = carma_context.instance()
        context.set_activeDeviceForCpy(self.sensors.device,1)
        self.sensors.d_wfs[n].comp_image()


    def get_offsets(self, int n):
        cdef carma_obj[float] *img
        cdef const long *cdims
        cdef np.ndarray[ndim=1,dtype=np.float32_t] data
        img=self.sensors.d_wfs[n].d_offsets
        cdims=img.getDims() 
        #data=np.zeros((cdims[1],cdims[2]),dtype=np.float32)
        data=np.empty((cdims[1]),dtype=np.float32)
        img.device2host(<float*>data.data)
        return data



    def get_imgtele(self, int n):
        cdef carma_host_obj[float] *img
        cdef const long *cdims
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data
        img=self.sensors.d_wfs[n].image_telemetry
        cdims=img.getDims() 
        #data=np.zeros((cdims[1],cdims[2]),dtype=np.float32)
        data=np.empty((cdims[1],cdims[2]),dtype=np.float32)
        
        #img.fill_into(<float*>data.data)
        return data
    
    cdef _get_binimg(self, int n):
        cdef carma_obj[float] *img
        cdef const long *cdims
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data
        img=self.sensors.d_wfs[n].d_binimg
        cdims=img.getDims() 
        data=np.empty((cdims[1],cdims[2]),dtype=np.float32)
        dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[n]).fill_binimage(0)
        img.device2host(<float*>data.data)
        return data

    def get_binimg(self,int n):
        return self._get_binimg(n)

    cdef _get_bincube(self, int n):
        cdef carma_obj[float] *cube
        cdef const long *cdims
        cdef np.ndarray[ndim=3,dtype=np.float32_t] data
        cube=self.sensors.d_wfs[n].d_bincube
        cdims=cube.getDims() 
        data=np.empty((cdims[1],cdims[2],cdims[3]),dtype=np.float32)
        cube.device2host(<float*>data.data)
        return data

    def get_bincube(self,int n):
        return self._get_bincube(n)

    def get_phase(self, int n):
        cdef carma_obj[float] *cube
        cdef const long *cdims
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data
        cube=self.sensors.d_wfs[n].d_gs.d_phase.d_screen
        cdims=cube.getDims() 
        data=np.empty((cdims[1],cdims[2]),dtype=np.float32)
        cube.device2host(<float*>data.data)
        return data


    def get_camplipup(self, int n):
        cdef carma_obj[cuFloatComplex] *cube
        cdef const long *cdims
        cdef np.ndarray[ndim=3,dtype=np.complex64_t] data
        cube=self.sensors.d_wfs[n].d_camplipup
        cdims=cube.getDims() 
        data=np.empty((cdims[3],cdims[2],cdims[1]),dtype=np.complex64)
        cube.device2host(<cuFloatComplex*>data.data)
        return data

    def get_amplifoc(self, int n):
        cdef carma_obj[cuFloatComplex] *cube
        cdef const long *cdims
        cdef np.ndarray[ndim=3,dtype=np.complex64_t] data
        cube=self.sensors.d_wfs[n].d_camplifoc
        cdims=cube.getDims() 
        data=np.empty((cdims[1],cdims[2],cdims[3]),dtype=np.complex64)
        cube.device2host(<cuFloatComplex*>data.data)
        return data


    def sensors_trace(self,int n, str type_trace, Atmos atmos, int rst): 
        """
        TODO dm 
        """

        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForce(self.sensors.device,1)
        if(type_trace=="all"):
            #TODO add dm
            self.sensors.d_wfs[n].d_gs.raytrace(atmos.s_a)
        elif(type_trace=="atmos"):
            self.sensors.d_wfs[n].d_gs.raytrace(atmos.s_a)
        #elif(type_trace=="dm"): # TODO dm
        #    self.sensors.d_wfs[n].sensors_trace(dm)



    cdef Bcast_dscreen(self):
        cdef carma_obj[float] *screen
        cdef float *ptr
        cdef int i,nsensors, size_i
        cdef long size

        nsensors=self.sensors.nsensors()
        for i in range(nsensors):
            screen=self.sensors.d_wfs[i].d_gs.d_phase.d_screen
            size=screen.getNbElem()
            size_i=size
            ptr=<float*>malloc(size*sizeof(float))
            screen.device2host(ptr)
            mpi.MPI_Barrier(mpi.MPI_COMM_WORLD)
            mpi.MPI_Bcast(ptr,size_i,mpi.MPI_FLOAT,0,mpi.MPI_COMM_WORLD)
            screen.host2device(ptr)
            free(ptr)
            mpi.MPI_Barrier(mpi.MPI_COMM_WORLD)



    cdef gather_bincube(self,MPI.Intracomm comm,int n):
        cdef np.ndarray[ndim=3,dtype=np.float32_t] tmp,bincube
        cdef np.ndarray[ndim=1,dtype=np.float32_t] tmp_flat,bincube_flat
        tmp=self._get_bincube(n).astype(np.float32)
        cdef int nx=self.sensors.d_wfs[n].npix
        cdef int ny=nx
        cdef int nz=self.sensors.d_wfs[n].nvalid
        cdef int nz_t=self.sensors.d_wfs[n].nvalid_tot
        cdef int size=nx*ny*nz
        cdef int recv_size=nx*ny*nz_t
        cdef int *count_bincube=self.sensors.d_wfs[n].count_bincube
        cdef int *displ_bincube=self.sensors.d_wfs[n].displ_bincube

        bincube=np.empty((nx,ny,nz_t),dtype=np.float32,order="F")
        bincube_flat=np.empty((nx*ny*nz_t),dtype=np.float32)
        tmp_flat=np.copy(tmp.flatten()[:nx*ny*nz])

        mpi.MPI_Gatherv(<float*>tmp_flat.data,size,mpi.MPI_FLOAT,
                    <float*>bincube_flat.data,
                    count_bincube  ,
                    displ_bincube ,
                    mpi.MPI_FLOAT,0,mpi.MPI_COMM_WORLD)

        if(self._get_rank(n)==0):
            bincube=np.reshape(bincube_flat,(nx,ny,nz_t),order="F")
            self.sensors.d_wfs[n].d_bincube.host2device(<float*>bincube.data)

    cdef _get_rank(self,int n):
        return self.sensors.d_wfs[n].rank

    def get_rank(self,int n):
        return self.sensors.d_wfs[n].rank
        
