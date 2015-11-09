include "../par.pxi"

import numpy as np
import time

import os

from cython.operator cimport dereference as deref, preincrement as inc

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

        """Add a centroider in the sutra_centroiders vector of the RTC on the GPU

        :parameters:
            sensor: (Sensors) : sutra_sensors object (GPU)

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


    def add_Controller(self, int nactu, float delay, bytes type_control, Dms dms,
                 list dmseen, np.ndarray[ndim=1,dtype=np.float32_t] alt,
                 int ndm, long Nphi=-1):
        """Add a controller in the sutra_controller vector of the RTC on the GPU

        :parameters:
            nactu: (int) : number of actuators

            delay: (float) : loop delay

            type_control: (str) : controller's type

            dms: (Dms) : sutra_dms object (GPU)

            type_dmseen: (char**) : dms indices controled by the controller

            alt: (np.ndarray[ndim=1,dtype=np.float32_t]) : altitudes of the dms seen

            ndm: (int) : number of dms controled

            Nphi: (long) : number of pixels in the pupil (used in geo controler case only)
        """
        type_dmseen=<char**>malloc(len(dmseen)*sizeof(char*))
        for j in range(controller.ndm.size):
                        type_dmseen[j]= dmseen[j]
                        
        self.add_controller( nactu, delay,
                            type_control, dms, type_dmseen,
                            alt,ndm,Nphi)
    
    cdef add_controller(self, int nactu, float delay, bytes type_control, Dms dms,
                 char **type_dmseen, np.ndarray[ndim=1,dtype=np.float32_t] alt,
                 int ndm, long Nphi=-1):
        """Add a controller in the sutra_controller vector of the RTC on the GPU

        :parameters:
            nactu: (int) : number of actuators

            delay: (float) : loop delay

            type_control: (str) : controller's type

            dms: (Dms) : sutra_dms object (GPU)

            type_dmseen: (char**) : dms indices controled by the controller

            alt: (np.ndarray[ndim=1,dtype=np.float32_t]) : altitudes of the dms seen

            ndm: (int) : number of dms controled

            Nphi: (long) : number of pixels in the pupil (used in geo controler case only)
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
        """Remove a controller"""
        cdef carma_context *context=carma_context.instance()
        context.set_activeDevice(self.device,1)
        self.dms.rm_controller()


    def setthresh(self, int ncentro, float thresh):
        """set threshold for the centroider #ncentro

        :parameters:
            ncentro: (int) : centroider's index

            thresh: (float) : threshold
        """

        cdef carma_context *context=carma_context.instance()
        cdef sutra_centroider_tcog *tcog=NULL
        cdef sutra_centroider *centro=NULL
        context.set_activeDeviceForCpy(self.device,1)
        if(self.rtc.d_centro[ncentro].is_type("tcog")):
            centro=self.rtc.d_centro.at(ncentro)
            tcog=dynamic_cast_centroider_tcog_ptr(centro)
            tcog.set_threshold(thresh)


    def setnmax(self,int ncentro, int nmax):
        """set the number of brightest pixels to consider for bpcog centroider

        :parameters:
            ncentro: (int) : centroider's index

            nmax: (int) : number of brightest pixels
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
        """Initialize npix in the sutra_centroider_corr object (useless ?)

        :parameters:
            ncentro: (int) : centroider's index
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


    def sensors_initweights(self,int ncentro, np.ndarray[ndim=2, dtype=np.float32_t] w):
        """Load the weight array in sutra_centroider_wcog object

        :parameters:
            ncentro: (int) : centroider's index

            w: (np.ndarray[ndim=2, dtype=np.float32_t]) : weight
        """
        cdef carma_context *context = carma_context.instance()
        context.set_activeDevice(self.rtc.device,1)
        cdef sutra_centroider_wcog *centro=\
            dynamic_cast_centroider_wcog_ptr(self.rtc.d_centro[ncentro])
        centro.init_weights()
        cdef np.ndarray w_F=w.flatten("F")
        centro.load_weights(<float*>w_F.data,int(w.ndim))


    def sensors_initcorr(self,int ncentro, np.ndarray[ndim=2,dtype=np.float32_t] w,
                        np.ndarray[ndim=2,dtype=np.float32_t] corr_norm,
                        int sizex, int sizey,
                        np.ndarray[ndim=2,dtype=np.float32_t] interpmat):
        """Initialize sutra_centroider_corr oblect

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
        """Return the centroids computed by the sutra_rtc object
        If ncontrol <= d_control.size, return rtc.d_centroids 
        Else, compute centroids from wfs[nwfs] with centroider[ncontrol]

        :parameters:
            ncontrol: (int) : controller's index

            g_wfs: (Sensors) : (optional) sutra_sensors object

            nwfs: (int) : (optional) number of wfs
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : centroids (arcsec)
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
        """Compute the centroids with sutra_controller #ncontrol object

        :parameters:
            ncontrol: (optional) controller's index
        """
        if(ncontrol>-1):
            self.rtc.do_centroids(ncontrol)
        else:
            self.rtc.do_centroids()

    cpdef docentroids_geom(self,int ncontrol=-1):
        """Compute the geometric centroids with sutra_controller #ncontrol object

        :parameters: 
            ncontrol: (optional) controller's index
        """
        if(ncontrol>-1):
            self.rtc.do_centroids_geom(ncontrol)
        else:
            self.rtc.do_centroids(0)


    cpdef init_proj(self,int ncontro,Dms dms,np.ndarray[ndim=1,dtype=np.int32_t] indx_dm,
            np.ndarray[ndim=1,dtype=np.float32_t] unitpervolt, np.ndarray[ndim=1,dtype=np.int32_t] indx_pup):
        """Initialize the projection matrix for sutra_controller_geo object. 
        The projection matrix is (IFt.IF)**(-1) * IFt where IF is the DMs influence functions matrix

        :parameters:
            ncontro: (int) : controller index

            dms: (Dms) : Dms object

            indx_dm: (np.ndarray[ndim=1,dtype=np.int32_t]) : indices of where(pup) on DM screen

            unitpervolt: (np.ndarray[ndim=1,dtype=np.float32_t]) : unitpervolt DM parameter

            indx_pup: (np.ndarray[ndim=1,dtype=np.int32_t]) : indices of where(pup) on ipupil screen
        """
        cdef carma_context *context=carma_context.instance()
        cdef sutra_controller_geo *controller_geo=\
            dynamic_cast_controller_geo_ptr(self.rtc.d_control[ncontro])

        context.set_activeDeviceForCpy(self.rtc.device,1)
        controller_geo.init_proj_sparse(dms.dms,<int*>indx_dm.data,
            <float*>unitpervolt.data, <int*>indx_pup.data)


    cpdef init_modalOpti(self,int ncontro,int nmodes,int nrec, np.ndarray[ndim=2,dtype=np.float32_t] M2V,
            float gmin, float gmax, int ngain, float Fs):
        """Initialize the modal optimization controller : compute the slopes-to-modes matrix 
        and the transfer functions

        :parameters:
            ncontro: (int) : controller index

            nmodes: (int) : number of modes

            nrec: (int) : number of recorded open slopes measurements

            M2V: (np.ndarray[ndim=2,dtype=np.float32_t]) : modes-to-volt matrix

            gmin: (float) : minimum gain for modal optimization

            gmax: (float) : maximum gain for modal optimization

            ngain: (int) : Number of tested gains

            Fs: (float) : sampling frequency
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls * controller_ls

        cdef np.ndarray[dtype=np.float32_t] M2V_F=M2V.flatten("F")

        if(<bytes>self.rtc.d_control[ncontro].get_type()=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.init_modalOpti(nmodes, nrec, <float*>M2V_F.data, gmin, gmax, ngain, Fs)
        else:
            raise TypeError("**** ERROR : Modal Optimization only for controller type ls ****")


    cpdef loadOpenLoop(self,int ncontro, np.ndarray[ndim=2, dtype=np.float32_t] ol_slopes):
        """Load an array of recoded open-loop measurements for modal optimization

        :parameters:
            ncontro: (int) : controller index

            ol_slopes: (np.ndarray[ndim=2, dtype=np.float32_t]) : open-loop slopes
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

    cpdef modalControlOptimization(self,int ncontro):
        """Compute the command matrix with modal control optimization

        :parameter:
            ncontro: controller index
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

    cpdef set_gain(self, int ncontro, float gain):
        """Set the loop gain in sutra_controller object
        
        :parameters:
            ncontro: (int) : controller index

            gain: (float) : loop gain
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


    cpdef set_mgain(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] mgain):
        """Set modal gains in sutra_controller object

        :parameters:
            ncontro: (int) : controller index

            mgain: (np.ndarray[ndim=1,dtype=np.float32_t]) : modal gains
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
    
    cpdef setCom(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] comvec):
        """Set the command vector of a sutra_controller object to comvec
        
        :parameters:
            ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef sutra_controller_ls *controller_ls
        cdef sutra_controller_mv *controller_mv
        cdef sutra_controller_geo *controller_geo
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()
        
        self.rtc.d_control[ncontro].d_com.host2device(<float*>comvec.data)

    cpdef get_mgain(self,int ncontro):
        """Return modal gains from sutra_controller

        :parameters:
            ncontro: (int) : controller index
        :return:
            mgain : (np.ndarray[ndim=1,dtype=np.float32_t]) : modal gains
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef sutra_controller_ls *controller_ls
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()
        cdef int size
        cdef np.ndarray[ndim=1,dtype=np.float32_t] mgain

        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            size=controller_ls.d_gain.getNbElem()
            mgain=np.zeros(size,dtype=np.float32)
            controller_ls.d_gain.device2host(<float*>mgain.data)
            return mgain
        elif(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            size=controller_mv.d_gain.getNbElem()
            mgain=np.zeros(size,dtype=np.float32)
            controller_mv.d_gain.device2host(<float*>mgain.data)
            return mgain
        else:
            raise TypeError("Controller needs to be ls or mv")


    cpdef set_imat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data):
        """Set the interaction matrix on a sutra_controller object

        :parameters:
            ncontro: (int) : controller index

            data: (np.ndarray[ndim=2,dtype=np.float32_t]) : interaction matrix to use
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
        """Return the interaction matrix of a sutra_controller object

        :parameters:
            ncontro: (int) : controller index
        :return:
            imat : (np.ndarray[ndim=2,dtype=np.float32_t]) : interaction matrix
        
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


    cpdef set_cmat(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data):
        """Set the command matrix on a sutra_controller object

        :parameters:
            ncontro: (int) : controller index

            data: (np.ndarray[ndim=2,dtype=np.float32_t]) : command matrix to use
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


            

    cpdef set_cmm(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] data):

        """Set the Cmm matrix on a sutra_controller_mv object

        :parameters:

            ncontro: (int) : controller index

            data: (np.ndarray[ndim=2,dtype=np.float32_t]) : Cmm matrix

        """

        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef np.ndarray[dtype=np.float32_t] data_F=data.flatten("F")

        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            controller_mv.d_Cmm.host2device(<float*>data_F.data)

        else:

            raise TypeError("Controller needs to be mv")


    cpdef get_cmm(self, int ncontro):
        """Return the Cmm matrix from a sutra_controller_mv object
        
        :parameters:
            ncontro: (int) : controller index
        :return:
            cmm : (np.ndarray[ndim=2,dtype=np.float32_t]) : Cmm matrix
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data_F, cmm
        cdef const long *cdims
        
        if(type_contro == "mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            cdims=controller_mv.d_Cmm.getDims()
            data_F=np.zeros((cdims[2],cdims[1]),dtype=np.float32)
            controller_mv.d_Cmm.device2host(<float*>data_F.data)
            cmm=np.reshape(data_F.flatten("F"),(cdims[1],cdims[2]))
            return cmm
        else:
            raise TypeError("Controller needs to be mv")

    cpdef get_cphim(self, int ncontro):
        """Return the Cphim matrix from a sutra_controller_mv object
        
        :parameters;
            ncontro: (int) : controller index
        :return:
            cphim : (np.ndarray[ndim=2,dtype=np.float32_t]) : Cphim matrix
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro = <bytes>self.rtc.d_control[ncontro].get_type()
        cdef np.ndarray[ndim=2,dtype=np.float32_t] data_F, cphim
        cdef const long *cdims
        
        if(type_contro == "mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            cdims=controller_mv.d_Cphim.getDims()
            data_F=np.zeros((cdims[2],cdims[1]),dtype=np.float32)
            controller_mv.d_Cphim.device2host(<float*>data_F.data)
            cphim=np.reshape(data_F.flatten("F"),(cdims[1],cdims[2]))
            return cphim
        else:
            raise TypeError("Controller needs to be mv")
            
    cpdef get_cmat(self,int ncontro):
        """Return the command matrix from a sutra_controller object

        :parameters:
            ncontro: (int) : controller index
        :return:
            cmat : (np.ndarray[ndim=2,dtype=np.float32_t]) : command matrix
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



    cpdef set_decayFactor(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] decay):
        """Set the decay factor on a sutra_controller_generic object

        :parameters:
            ncontro: (int) : controller index

            decay: (np.ndarray[ndim=1,dtype=np.float32_t]) : ask to Rico
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


    cpdef set_matE(self,int ncontro, np.ndarray[ndim=2,dtype=np.float32_t] matE):
        """Set the matrix E on a sutra_controller_generic object

        :parameters:
            ncontro: (int) : controller index

            matE: (np.ndarray[ndim=2,dtype=np.float32_t]) : ask to Rico
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


    cpdef doimat_geom(self, int ncontro, Dms g_dms,int geom):
        """Compute the interaction matrix by using a geometric centroiding method

        :parameters:
            ncontro: (int) : controller index

            g_dms: (Dms) : Dms object

            geom: (int) : type of geometric method (0 or 1)
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        self.rtc.do_imat_geom(ncontro,g_dms.dms,geom)
        print "TODO call imat_geom"
        print "TODO set_imat"

    cpdef doimat(self, int ncontro, Dms g_dms):
        """Compute the interaction matrix

        :parameters:
            ncontro: (int) : controller index

            g_dms: (Dms) : Dms object
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
        cdef int nactu

        cdef int rank
        IF USE_MPI==1:
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD,&rank)
        ELSE:
            rank=0
        nactu = d_imat.getDims(2)
        cc = 0
        it_dm=control.d_dmseen.begin()
        print "Doing imat...%d%%"%cc,
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
                cc = cc + 1
                print "\rDoing imat...%d%%"%((cc*100)/nactu),
            inc(it_dm)
        print "\n"


    cpdef sensors_compslopes(self, int ncentro, int nmax=-1, float thresh=-1):
        """Compute the slopes in a sutra_wfs object. This function is equivalent to
        docentroids() but the centroids are stored in the sutra_wfs object instead of
        the sutra_rtc object

        :parameters:
            ncentro: (int) : centroider index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        self.rtc.d_centro[ncentro].get_cog()


    cpdef imat_svd(self,int ncontro):
        """Compute the singular value decomposition of the interaction matrix

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


    cpdef setU(self,int ncontro,np.ndarray[ndim=2,dtype=np.float32_t] U):
        """Set the eigen modes matrix of the imat decomposition in a sutra_controller_ls object

        :parameters:
            ncontro: (int) : controller index

            U: (np.ndarray[ndim=2,dtype=np.float32_t]) : eigen modes matrix
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)

        cdef sutra_controller_ls *controller_ls
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()

        cdef np.ndarray[dtype=np.float32_t] data_F= U.flatten("F")
        if(type_contro=="ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            controller_ls.d_U.host2device(<float*>data_F.data)


    cpdef setEigenvals(self, int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] eigenvals):
        """Set the eigen values of the imat decomposition in a sutra_controller_ls object

        :parameters:
            ncontro: (int) : controller index

            eigenvals: (np.ndarray[ndim=1,dtype=np.float32_t]) : eigen values
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


    cpdef getU(self, int ncontro):
        """Return the eigen modes matrix of the imat decomposition from a sutra_controller_ls object

        :parameters:
            ncontro: (int) : controller index
        :return:
            U : (np.ndarray[ndim=2,dtype=np.float32_t]) : eigen modes matrix
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


    cpdef getEigenvals(self,int ncontro):
        """Return the eigen values of the imat decomposition in a sutra_controller object

        :parameters:
            ncontro: (int) : controller index
        :return:
            eigenvals : (np.ndarray[ndim=1,dtype=np.float32_t]) : eigenvalues
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_ls *controller_ls
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()
        cdef np.ndarray[ndim=1, dtype=np.float32_t] data
        cdef const long *dims
        if(type_contro=="ls"):
            controller_ls=dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontro])
            dims=controller_ls.h_eigenvals.getDims()
            data=np.zeros((dims[1]),dtype=np.float32)
            controller_ls.h_eigenvals.fill_into(<float*>data.data)
        if(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            dims=controller_mv.h_Cmmeigenvals.getDims()
            data=np.zeros((dims[1]),dtype=np.float32)
            controller_mv.h_eigenvals.fill_into(<float*>data.data)


        return data


    cpdef getCmmEigenvals(self,int ncontro):
        """Return the eigen values of the Cmm decomposition in a sutra_controller_mv object
        :parameters:
            ncontro: (int) : controller index
        :return:
            eigenvals : (np.ndarray[ndim=1,dtype=np.float32_t]) : eigenvalues
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()
        cdef np.ndarray[ndim=1, dtype=np.float32_t] data
        cdef const long *dims
        if(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            dims=controller_mv.h_Cmmeigenvals.getDims()
            data=np.zeros((dims[1]),dtype=np.float32)
            controller_mv.h_Cmmeigenvals.fill_into(<float*>data.data)

        return data


    cpdef getCenbuff(self, int ncontro):
        """Return the centroids buffer from a sutra_controller_ls object. This 
        buffer contains centroids from iteration i-delay to current iteration
        :parameters:
            ncontro: (int) : controller index
        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : centroids buffer
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


    cpdef getErr(self,int ncontro):
        """Return the command increment (cmat*slopes) from a sutra_controller_ls object

        :parameters:
            ncontro: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : command increment
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
        """Return the command vector from a sutra_controller object
        
        :parameters:
            ncontro: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : command vector
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef np.ndarray[ndim=1, dtype=np.float32_t] data
        cdef const long *dims
        dims=self.rtc.d_control[ncontro].d_com.getDims()
        data=np.zeros((dims[1]),dtype=np.float32)
        self.rtc.d_control[ncontro].d_com.device2host(<float*>data.data)

        return data

    cpdef getolmeas(self,int ncontro):
        """Return the reconstructed open-loop measurement from a sutra_controller_mv object
        
        :parameters:
            ncontro: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : reconstructed open-loop
        """
        cdef carma_context *context=carma_context.instance()
        cdef sutra_controller_mv *controller_mv
        
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef np.ndarray[ndim=1, dtype=np.float32_t] data
        cdef const long *dims
        
        controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
        dims=controller_mv.d_olmeas.getDims()
        data=np.zeros((dims[1]),dtype=np.float32)
        controller_mv.d_olmeas.device2host(<float*>data.data)

        return data


    cpdef getVoltage(self,int ncontro):
        """Return the voltage vector that will be effectively applied to the DMs
        
        :parameters:
            ncontro: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : voltage vector
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
        """Return the centroids computed by the sutra_rtc object

        :parameters:
            ncontrol: (int) : controller's index

        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : centroids (arcsec)
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


    cpdef buildcmat(self,int ncontro,int nfilt, int filt_tt=0):
        """Compute the command matrix in a sutra_controller_ls object

        :parameters:
            ncontro: (int) : controller index

            nfilt: (int) : number of modes to filter

            filt_tt: (int) : (optional) flag to filter TT
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


    cpdef buildcmatmv(self,int ncontro,float cond):
        """Compute the command matrix in a sutra_controller_mv object

        :parameters:
            ncontro: (int) : controller index

            cond: (float) : conditioning factor for the Cmm inversion
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            controller_mv.build_cmat(cond)



    cpdef loadnoisemat(self,int ncontro, np.ndarray[ndim=1,dtype=np.float32_t] N):
        """Load the noise vector on a sutra_controller_mv object

        :parameters:
            ncontro: (int) : controller index

            N: (np.ndarray[ndim=1,dtype=np.float32_t]) : noise vector
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device,1)
        cdef sutra_controller_mv *controller_mv
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()

        if(type_contro=="mv"):
            controller_mv=dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontro])
            controller_mv.load_noisemat(<float*>N.data)



    cpdef docontrol(self,int ncontro):
        """Compute the command to apply on the DMs on a sutra_controller object

        :parameters:
            ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        context.set_activeDevice(self.rtc.device,1)
        self.rtc.do_control(ncontro)
        
    cpdef docontrol_geo(self, int ncontro, Dms dms, Target target, int ntarget):
        """Compute the command to apply on the DMs on a sutra_controller_geo object

        :parameters:
            ncontro: (int) : controller index
        """
        cdef carma_context *context=carma_context.instance()
        cdef sutra_controller_geo *controller_geo
        cdef bytes type_contro=<bytes>self.rtc.d_control[ncontro].get_type()
        context.set_activeDevice(self.rtc.device,1)
        
        if(type_contro == "geo"):
            controller_geo = dynamic_cast_controller_geo_ptr(self.rtc.d_control[ncontro])
            controller_geo.comp_dphi(target.target.d_targets[ntarget])
            self.rtc.do_control(ncontro)
        else:
            raise TypeError("Controller needs to be geo")
        

    cpdef applycontrol(self,int ncontro,Dms dms):
        """Compute the DMs shapes from the commands computed in a sutra_controller_object.
        From the command vector, it computes the voltage command (adding pertrubation voltages,
        taking delay into account) and then apply it to the dms
        
        :parameters:
            ncontro: (int) : controller index
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
            Target g_tar, Param_target p_tar, clean=1, brama=None, doimat=None,simul_name="", int overwrite=1):
    """Initialize all the sutra_rtc objects : centroiders and controllers

    :parameters:
        g_wfs: (Sensors) : Sensors object

        p_wfs: (list of Param_wfs) : wfs settings

        g_dms: (Dms) : Dms object

        p_dms: (list of Param_dms) : dms settings

        p_geom: (Param_geom) : geom settings

        p_atmos: (Param_atmos) : atmos settings

        g_atmos: (Atmos) : Atmos object

        p_tel: (Param_tel) : telescope settings

        p_loop: (Param_loop) : loop settings

        p_tar: (Param_target) : (optional) target settings

        clean: (int) : (optional) clean datafiles (imat, U, eigenv)

        brama: (int) : (optional) not implemented yet

        doimat: (int) : (optional) force imat computation

        simul_name: (str) : (optional) simulation's name, use for path to save data (imat, U...)
       
        overwrite: (int) : (optional) overwrite data files if overwite=1 (default 1)

    :return:
        Rtc : (Rtc) : Rtc object

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
    #cdef int nwfs =0
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
                            print "reading Na profile from",shesha_savepath+profilename
                            prof=np.load(shesha_savepath+profilename).astype(np.float32)

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
                    correct_dm(p_dms,g_dms,controller,p_geom, imat,simul_name,overwrite)
                    #TODO add timer
                    if(controller.type_control !="geo"):
                        nwfs=controller.nwfs
                        if(len(p_wfs)==1):
                            nwfs=p_rtc.controllers[i].nwfs
                            # TODO fixing a bug ... still not understood
                        nvalid = sum([p_wfs[k]._nvalid for k in nwfs])
                        controller.set_nvalid(nvalid)
                    #parameter for add_controller(_geo)
                    ndms=controller.ndm.tolist()
                    controller.set_nactu([p_dms[n]._ntotact for n in ndms])
                    nactu=np.sum([p_dms[j]._ntotact for j in ndms])
                    alt=np.array([p_dms[j].alt for j in controller.ndm],
                                dtype=np.float32)

                    list_dmseen=[p_dms[j].type_dm for j in controller.ndm]

                    type_dmseen=<char**>malloc(controller.ndm.size*sizeof(char*))
                    #TODO : j is not always a range, it has to be index set in controller.ndm
                    for j in range(controller.ndm.size):
                        type_dmseen[j]= p_dms[controller.ndm[j]].type_dm

                    context.set_activeDeviceForCpy(device,1)
                    if(controller.type_control=="geo"):
                        Nphi=np.where(p_geom._spupil)[0].size
                        #list_dmseen,alt,controller.ndm.size
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
                            indx_dm[cpt:cpt+np.where(pup_dm)[0].size] = np.where(pup_dm.flatten('F'))[0] 
                            cpt+=np.where(pup_dm)[0].size
                        #convert unitpervolt list to a np.ndarray
                        unitpervolt=np.array([p_dms[j].unitpervolt for j in range(len(p_dms))],
                                    dtype=np.float32)
                        
                        g_rtc.init_proj(i, g_dms, indx_dm, unitpervolt, indx_pup)

                    free(type_dmseen)

                    if(controller.type_control=="ls"):
                        if(doimat):
                            imat_init(i,g_rtc,p_rtc,g_dms,g_wfs,p_wfs,p_tel,clean=clean,simul_name=simul_name,overwrite=overwrite)
                            if(controller.modopti == 1):
                                print "Initializing Modal Optimization : "
                                if(controller.nrec==0):
                                    controller.nrec=2048
                                else:
                                    #next power of 2 (for fft)
                                    controller.nrec=int(2**np.ceil(np.log2(controller.nrec)))
                                if(controller.nmodes==0):
                                    controller.nmodes=sum(p_dms[j]._ntotact for k in ndms)
                                if(controller.gmax==0):
                                    controller.gmax=1.0
                                if(controller.ngain==0):
                                    controller.ngain=15
                                KL2V = compute_KL2V(controller,g_dms,p_dms,p_geom,p_atmos,p_tel)
                                g_rtc.init_modalOpti(i,controller.nmodes,controller.nrec,KL2V,
                                    controller.gmin,controller.gmax,controller.ngain,1./p_loop.ittime)
                                ol_slopes=openLoopSlp(g_atmos,g_rtc,controller.nrec,i,g_wfs,p_wfs,p_tar,g_tar)
                                g_rtc.loadOpenLoop(i,ol_slopes)
                                g_rtc.modalControlOptimization(i)
                            else:
                                cmat_init(i,g_rtc,p_rtc,p_wfs,p_atmos,p_tel, p_dms, clean=clean,simul_name=simul_name,overwrite=overwrite)
                                g_rtc.set_gain(i,controller.gain)
                                mgain=np.ones(sum([p_dms[j]._ntotact for j in range(len(p_dms))]),dtype=np.float32)
                                #filtering tilt ...
                                #mgain(-1:0) = 0.0f;
                                g_rtc.set_mgain(i,mgain)
                        else:
                            nactu=np.sum(controller.nactu)
                            nvalid=np.sum(controller.nvalid)
                            imat=np.zeros((nactu*nvalid*2),dtype=np.float32)
                            g_rtc.set_imat(i,imat)
                            g_rtc.set_cmat(i,imat)
                            g_rtc.set_gain(i,controller.gain)
                            mgain=np.ones(sum([p_dms[j]._ntotact for j in range(len(p_dms))]),dtype=np.float32)
                            # filtering tilt ...
                            g_rtc.set_mgain(i,mgain)

                    if(controller.type_control=="cured"):
                        print "initializing cured controller"
                        if( "tt" in [p_dms[j].type_dm for j in range(len(p_dms))]):
                            tt_flag=1
                        else:
                            tt_flag=0
                        print "TODO controller_initcured"
                        #controller_initcured,g_rtc,0,int(y_wfs(1).nxsub),int(*y_wfs(1)._isvalid),
                                #int(controllers(i).cured_ndivs),int(tt_flag);
                        g_rtc.set_gain(i,controller.gain)
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

                        g_rtc.set_gain(i,controller.gain)

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

                    if(controller.type_control == "mv"):
                        imat = imat_geom(g_wfs,p_wfs,controller,g_dms,p_dms,meth=0)
                        #imat_init(i,g_rtc,p_rtc,g_dms,g_wfs,p_wfs,p_tel,clean=1,simul_name=simul_name)
                        g_rtc.set_imat(i,imat)
                        g_rtc.set_gain(i,controller.gain)
                        size=sum([p_dms[j]._ntotact for j in range(len(p_dms))])
                        mgain = np.ones(size,dtype=np.float32)
                        g_rtc.set_mgain(i,mgain)
                        doTomoMatrices(i,g_rtc,p_wfs,g_dms,g_atmos,g_wfs,p_rtc,p_geom,p_dms,p_tel,p_atmos)
                        cmat_init(i,g_rtc,p_rtc,p_wfs,p_atmos,p_tel,p_dms,clean=1,simul_name=simul_name,overwrite=overwrite)
                        
                    elif(controller.type_control=="generic"):
                        size=sum([p_dms[j]._ntotact for j in range(len(p_dms))])
                        decayFactor=np.zeros(size ,dtype=np.float32)
                        mgain=np.zeros(size ,dtype=np.float32)
                        matE=np.zeros((size,size) ,dtype=np.float32)
                        cmat=np.zeros((size,np.sum(controller.nvalid)*2),
                                       dtype=np.float32)
                        #TODO ? law = "integrator";

                        g_rtc.set_decayFactor(i,decayFactor)
                        g_rtc.set_mgain(i,mgain)
                        g_rtc.set_cmat(i,cmat)
                        g_rtc.set_matE(i,matE)
            
    return g_rtc
    

cpdef correct_dm(p_dms, Dms g_dms, Param_controller p_control, Param_geom p_geom, np.ndarray imat, bytes simul_name, int overwrite=1):
    """Correct the geometry of the DMs using the imat (filter unseen actuators)

    :parameters:
        p_dms: (list of Param_dm) : dms settings

        g_dms: (Dms) : Dms object

        p_control: (Param_controller) : controller settings

        p_geom: (Param_geom) : geom settings

        imat: (np.ndarray) : interaction matrix

        simul_name: (str) : simulation's name, use for data files' path

        overwrite: (int) : (optional) overwrite data files if overwite=1 (default 1)
    """
    #cdef carma_context *context = carma_context.instance() #UNUSED
    
    cdef int i, nm, nmc, inds,ndm,nactu_nm
    cdef np.ndarray[ndim=1,dtype=np.float32_t] resp

    cdef bytes filename
    cdef bytes dirsave = shesha_savepath+<bytes>"mat/"

    cdef long dims,ninflu,influsize,NR,NP

    if (simul_name=="" ): 
        imat_clean = 1
    else:
        imat_clean=0

    ndm=p_control.ndm.size
    for i in range(ndm):
        nm=p_control.ndm[i]
        g_dms.remove_dm(p_dms[nm].type_dm,p_dms[nm].alt)

    resp=np.sqrt(np.sum(imat**2,axis=0))

    inds=0

    cdef int rank
    IF USE_MPI==1:
        mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD,&rank)
    ELSE:
        rank=0

    for nmc in range(ndm):
        nm=p_control.ndm[nmc]
        nactu_nm=p_dms[nm]._ntotact
        # filter actuators only in stackarray mirrors:

        if(p_dms[nm].type_dm=="pzt"):
            filename=dirsave+"pztok-"+str(nm)+"-"+simul_name+".npy"

            if(imat_clean<1 and os.path.isfile(filename)):
                print "reading valid actuator from:",filename
                ok=np.load(filename)
                filename=dirsave+"pztnok-"+str(nm)+"-"+simul_name+".npy"
                print "reading unvalid actuator from:",filename
                nok=np.load(filename)
            else:
                tmp=resp[inds:inds+p_dms[nm]._ntotact]
                ok=np.where(tmp.flatten("F")>p_dms[nm].thresh*np.max(tmp))[0]
                nok=np.where(tmp.flatten("F")<=p_dms[nm].thresh*np.max(tmp))[0]
                if(simul_name !="" and rank==0 and overwrite==1):
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



cpdef imat_geom(Sensors g_wfs, p_wfs, Param_controller p_control,Dms g_dms, p_dms, int meth=0):
    """Compute the interaction matrix with a geometric method

    :parameters:
        g_wfs: (Sensors) : Sensors object

        p_wfs: (list of Param_wfs) : wfs settings

        p_control: (Param_controller) : controller settings

        g_dms: (Dms) : Dms object

        p_dms: (list of Param_dm) : dms settings

        meth: (int) : (optional) method type (0 or 1)
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
        nm=p_control.nwfs[nw]
        imat_size1+=p_wfs[nw]._nvalid*2

    for nmc in range(ndm):
        imat_size2 +=p_dms[nmc]._ntotact

    imat_cpu = np.zeros((imat_size1,imat_size2),dtype=np.float32)
    ind=0
    cc = 0
    print "Doing imat geom...%d%%"%cc,
    for nmc in range(ndm):
        nm=p_control.ndm[nmc]
        g_dms.resetdm(p_dms[nm].type_dm,p_dms[nm].alt)
        for i in range(p_dms[nm]._ntotact):
            g_dms.oneactu(p_dms[nm].type_dm,p_dms[nm].alt,i,p_dms[nm].push4imat)
            nslps=0
            for nw in range(nwfs):
                wfs=p_control.nwfs[nw]
                g_wfs.sensors_trace(wfs,"dm",None,dms=g_dms,rst=1)
                g_wfs.slopes_geom(wfs,meth)
                imat_cpu[nslps:nslps+p_wfs[wfs]._nvalid*2,ind]=g_wfs._get_slopes(wfs)
                nslps+=p_wfs[wfs]._nvalid*2
            imat_cpu[:,ind]=imat_cpu[:,ind]/float(p_dms[nm].push4imat)
            ind=ind+1
            cc = cc+1
            g_dms.resetdm(p_dms[nm].type_dm,p_dms[nm].alt)
            
            print "\r Doing imat geom...%d%%"%((cc*100/imat_size2)),
    print "\n"
    return imat_cpu


cpdef manual_imat(Rtc g_rtc,Sensors g_wfs, p_wfs, Dms g_dms, p_dms):
    """Compute the interaction matrix 'manually', ie without sutra_rtc doimat method

    :parameters:
        g_rtc: (Rtc) : Rtc object

        g_wfs: (Sensors) : Sensors object

        g_dms: (Dms) : Dms object

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


cpdef get_r0(float r0_at_lambda1, float lambda1, float lambda2):
    """Compute r0 at lambda2 from r0 value at lambda1 

    :parameters:
        r0_at_lambda1: (float) : r0 value at lambda1

        lambda1: (float) : lambda1

        lambda2: (float) : lambda2
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



    

cpdef compute_KL2V(Param_controller controller, Dms dms, p_dms, Param_geom p_geom, Param_atmos p_atmos, Param_tel p_tel):
    """Compute the Karhunen-Loeve to Volt matrix 
    (transfer matrix between the KL space and volt space for a pzt dm)

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
            KL2V[indx_act:indx_act+ntotact[ndm],indx_act:indx_act+ntotact[ndm]]=\
                compute_klbasis( dms, p_dms[ndm],p_geom,p_atmos,p_tel)

            print "compute klbasis done"
        elif(p_dms[ndm].type_dm=="tt"):
            nTT+=1

    KL2V=KL2V[:,:controller.nmodes]
    
    if(nTT!=0):
        KL2V[:,:controller.nmodes-2]= KL2V[:,2:]
        KL2V[:,controller.nmodes-2:]= np.zeros((np.sum(ntotact),2),dtype=np.float32)
        KL2V[np.sum(ntotact)-2:,controller.nmodes-2:]=np.identity(2,dtype=np.float32)

    return KL2V

cpdef openLoopSlp(Atmos g_atm, Rtc g_rtc,int nrec, int ncontro, Sensors g_wfs,  
        p_wfs, Param_target p_tar,Target g_tar):
    """Return a set of recorded open-loop slopes, usefull for modal control optimization

    :parameters:
        g_atm: (Atmos) : Atmos object

        g_rtc: (Rtc) : Rtc object

        nrec: (int) : number of samples to record

        ncontro: (int) : controller's index

        g_wfs: (Sensors) : Sensors object

        p_wfs: (list of Param_wfs) : wfs settings

        p_tar: (Param_target) : target settings

        g_tar: (Target) : Target object
    """
    #TEST IT
    cdef int i,j
    cdef np.ndarray[ndim=2,dtype=np.float32_t] ol_slopes=\
        np.zeros((sum([2*p_wfs[i]._nvalid for i in range(len(p_wfs))]),nrec),
        dtype=np.float32)
     
    print "Recording "+str(nrec)+" open-loop slopes..."
    for i in range(nrec):
        if(g_tar is not None):
            g_atm.move_atmos()
            if(p_tar is not None):
                for j in range(p_tar.ntargets):
                    g_tar.atmos_trace(j,g_atm)

        if(p_wfs is not None and g_wfs is not None):
            for j in range(len(p_wfs)):
                g_wfs.sensors_trace(j,"atmos",g_atm)
                g_wfs.sensors_compimg(j)
                g_rtc.sensors_compslopes(ncontro)
                ol_slopes[j*p_wfs[j]._nvalid*2:(j+1)*p_wfs[j]._nvalid*2,i]=g_wfs._get_slopes(j)
    print "done"
    return ol_slopes

cpdef imat_init(int ncontro, Rtc g_rtc, Param_rtc p_rtc, Dms g_dms, Sensors g_wfs,
        p_wfs, Param_tel p_tel, int clean=1, bytes simul_name=<bytes>"", overwrite=1):
    """Initialize and compute the interaction matrix on the GPU

    :parameters:
        ncontro: (int) : controller's index

        g_rtc: (Rtc) : Rtc object

        p_rtc: (Param_rtc) : rtc settings

        g_dms: (Dms) : Dms object

        g_wfs: (Sensors) : Sensors object

        p_wfs: (list of Param_wfs) : wfs settings

        p_tel: (Param_tel) : telescope settings

        clean: (int) : (optional) : clean datafiles (imat, U, eigenv)

        simul_name: (str) : (optional) simulation's name, use for data files' path

        overwrite: (int) : (optional
    """
    cdef bytes dirsave=shesha_savepath+<bytes>"mat/"
    cdef bytes filename=dirsave+"imat-"+str(ncontro)+"-"+simul_name+".npy"
    cdef bytes profilename=shesha_savepath+<bytes>"allProfileNa_withAltitude_1Gaussian.npy"
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
                print "reading Na profile from:", profilename
                prof=np.load(profilename)
                h=prof[0,:]
                prof=prof[1:,:]
                prof=np.mean(prof,axis=0)
                prep_lgs_prof(p_wfs[i],i,p_tel,prof,h,
                                        p_wfs[i].beamsize,g_wfs,<bytes>"",imat=1)

        t0=time.time()
        g_rtc.doimat(ncontro,g_dms)
        print "done in ",time.time()-t0
        p_rtc.controllers[ncontro].set_imat(g_rtc.get_imat(ncontro))
        if(simul_name!="" and rank==0 and overwrite==1):
            np.save(filename,p_rtc.controllers[ncontro].imat)

    else:
        print "reading imat from:",filename
        p_rtc.controllers[ncontro].set_imat(np.load(filename))
        g_rtc.set_imat(ncontro, p_rtc.controllers[ncontro].imat)

    #now restore original profile in lgs spots
    for i in range(len(p_wfs)):
        if(p_wfs[i].gsalt>0):
            prep_lgs_prof(p_wfs[i],i,p_tel,p_wfs[i]._profna,p_wfs[i]._altna,
                            p_wfs[i].beamsize,g_wfs)




cpdef cmat_init(int ncontro, Rtc g_rtc, Param_rtc p_rtc, list p_wfs, Param_atmos p_atmos,
               Param_tel p_tel, p_dms, clean=1, bytes simul_name=<bytes>"", int overwrite=1):
    """ Compute the command matrix on the GPU

    :parameters:
        ncontro: (int) :

        g_rtc: (Rtc) :

        p_rtc: (Param_rtc) : rtc settings

        p_wfs: (list of Param_wfs) : wfs settings

        p_tel : (Param_tel) : telescope settings

        clean: (int) : (optional) clean datafiles (imat, U, eigenv)

        simul_name: (str) : (optional) simulation's name, use for data files' path

        overwrite: (int) : (optional) overwrite data files if overwite=1 (default 1)
        
    """

    cdef bytes dirsave=shesha_savepath+<bytes>"mat/"
    cdef bytes filename

    cdef int cmat_clean
    cdef double t0
    cdef np.ndarray[ndim=1,dtype=np.float32_t] eigenv, N
    cdef np.ndarray[ndim=2,dtype=np.float32_t] U, imat
    cdef np.ndarray[ndim=1,dtype=np.int64_t] mfilt
    
    cdef sutra_controller_mv *controller_mv
    

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
            if(simul_name!="" and rank==0 and overwrite==1):
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
        
        controller_mv = dynamic_cast_controller_mv_ptr(g_rtc.rtc.d_control[ncontro])
        nvalidperwfs = np.array([o._nvalid   for o in p_wfs],dtype=np.int64)
        N=np.zeros(2*sum([nvalidperwfs[j]   for j in p_rtc.controllers[ncontro].nwfs]),
                    dtype=np.float32)
        ind=0
        for k in p_rtc.controllers[ncontro].nwfs:
            N[ind:ind+2*p_wfs[k]._nvalid] = noise_cov(k,p_wfs[k],p_atmos,p_tel)
            ind+=2*p_wfs[k]._nvalid

        g_rtc.loadnoisemat(ncontro,N)
        print "Building cmat..."
        g_rtc.buildcmatmv(ncontro,p_rtc.controllers[ncontro].maxcond)
        
        if(p_rtc.controllers[ncontro].TTcond==0):
            p_rtc.controllers[ncontro].set_TTcond(p_rtc.controllers[ncontro].maxcond)

        if("tt" in [dm.type_dm for dm in p_dms]):
            controller_mv.filter_cmat(p_rtc.controllers[ncontro].TTcond)
        

    p_rtc.controllers[ncontro].set_cmat(g_rtc.get_cmat(ncontro))

cpdef doTomoMatrices(int ncontro, Rtc g_rtc, list wfs, Dms g_dm, Atmos g_atmos, Sensors g_wfs,  Param_rtc p_rtc, Param_geom geom, list p_dm, Param_tel p_tel, Param_atmos p_atmos):
    """Compute Cmm and Cphim matrices for the MV controller on GPU

    :parameters:
        g_wfs: (Sensors) :

        p_wfs: (list of Param_wfs) : wfs settings

        g_dms: (Dms) :

        p_dms: (list of Param_dms) : dms settings

        p_geom: (Param_geom) : geom settings

        p_atmos: (Param_atmos) : atmos settings

        g_atmos: (Atmos) :

        p_tel: (Param_tel) : telescope settings

    """
    cdef np.ndarray[ndim=2,dtype=np.float32_t] ipup, spup, F, Nact
    cdef np.ndarray nvalidperwfs = np.array([o._nvalid   for o in wfs],dtype=np.int64)
    cdef np.ndarray[ndim=1, dtype=np.float64_t] frac_d, L0_d, X, Y, Xactu, Yactu, k2, pitch, alphaX, alphaY, alt_DM
    cdef np.ndarray[ndim=1, dtype=np.int64_t] indlayersDM, NlayersDM
    cdef sutra_controller_mv *controller_mv
    
    # Bring bottom left corner of valid subapertures in ipupil
    ipup = geom.get_ipupil()
    spup = geom.get_spupil()
    s2ipup = (ipup.shape[0] - spup.shape[0])/2.
    nvalid = sum([nvalidperwfs[j]   for j in p_rtc.controllers[ncontro].nwfs]) # Total number of subapertures
    ind = 0
    X = np.zeros(nvalid,dtype=np.float64) # X-position of the bottom left corner of each valid subaperture
    Y = np.zeros(nvalid,dtype=np.float64) # Y-position of the bottom left corner of each subaperture

    for k in p_rtc.controllers[ncontro].nwfs:
        posx = wfs[k]._istart + s2ipup 
        posx = posx * wfs[k]._isvalid # X-position of the bottom left corner of each valid subaperture
        posx = posx[np.where(posx > 0)] - ipup.shape[0]/2 - 1 #Select the valid ones, bring the origin in the center of ipupil and 0-index it
        posy = wfs[k]._jstart + s2ipup
        posy = posy * wfs[k]._isvalid
        posy = posy.T[np.where(posy > 0)] - ipup.shape[0]/2 - 1
        sspDiam = posx[1] - posx[0] # Diameter of one ssp in pixels
        p2m = (p_tel.diam / wfs[k].nxsub) / sspDiam # Size of one pixel in meters
        posx *= p2m # Positions in meters
        posy *= p2m
        X[ind:ind+wfs[k]._nvalid] = posx
        Y[ind:ind+wfs[k]._nvalid] = posy
        ind += wfs[k]._nvalid
    
    # Get the total number of pzt DM and actuators to control
    nactu = 0 
    npzt = 0
    for k in p_rtc.controllers[ncontro].ndm:
        if(p_dm[k].type_dm == "pzt"):
            nactu += p_dm[k]._ntotact
            npzt += 1
        
    
    Xactu = np.zeros(nactu,dtype=np.float64) # X-position actuators in ipupil
    Yactu = np.zeros(nactu,dtype=np.float64) # Y-position actuators in ipupil
    k2 = np.zeros(npzt,dtype=np.float64) # k2 factors for computation
    pitch = np.zeros(npzt,dtype=np.float64)
    alt_DM = np.zeros(npzt,dtype=np.float64)
    ind = 0
    indk = 0
    for k in p_rtc.controllers[ncontro].ndm:
        if(p_dm[k].type_dm == "pzt"):
            p2m = p_tel.diam / geom.pupdiam
            actu_x = (p_dm[k]._xpos - ipup.shape[0]/2)*p2m # Conversion in meters in the center of ipupil
            actu_y = (p_dm[k]._ypos - ipup.shape[0]/2)*p2m
            pitch[indk] = actu_x[1] - actu_x[0]
            k2[indk] = wfs[0].Lambda / 2. / np.pi / p_dm[k].unitpervolt
            alt_DM[indk] = <double>p_dm[k].alt
            Xactu[ind:ind+p_dm[k]._ntotact] = actu_x
            Yactu[ind:ind+p_dm[k]._ntotact] = actu_y
            
            ind += p_dm[k]._ntotact
            indk += 1
    
    # Select a DM for each layer of atmos
    NlayersDM = np.zeros(npzt,dtype=np.int64) #Useless for now
    indlayersDM = selectDMforLayers(ncontro, p_atmos, p_rtc, p_dm)
   # print "indlayer = ",indlayersDM
    
    # Get FoV
    #RASC = 180.0/np.pi * 3600.
    wfs_distance = np.zeros(len(p_rtc.controllers[ncontro].nwfs),dtype = np.float64)
    ind = 0
    for k in p_rtc.controllers[ncontro].nwfs:
        wfs_distance[ind] = np.sqrt(wfs[k].xpos**2 + wfs[k].ypos**2)
        ind += 1
    FoV = np.max(wfs_distance) / RASC
    
    # WFS postions in rad
    alphaX = np.zeros(len(p_rtc.controllers[ncontro].nwfs))
    alphaY = np.zeros(len(p_rtc.controllers[ncontro].nwfs))
    
    ind = 0
    for k in p_rtc.controllers[ncontro].nwfs:
        alphaX[ind] = wfs[k].xpos / RASC
        alphaY[ind] = wfs[k].ypos / RASC
        ind +=1
    
    controller_mv = dynamic_cast_controller_mv_ptr(g_rtc.rtc.d_control[ncontro])
    
    L0_d = np.copy(p_atmos.L0).astype(np.float64)
    frac_d = np.copy(p_atmos.frac * (p_atmos.r0 ** (-5.0/3.0))).astype(np.float64)
    
    print "Computing Cphim..."
    controller_mv.compute_Cphim(g_atmos.s_a, g_wfs.sensors, g_dm.dms, <double*>L0_d.data, 
                                <double*>frac_d.data, <double*>alphaX.data, 
                                <double*>alphaY.data,<double*>X.data,<double*>Y.data,
                                <double*>Xactu.data,<double*>Yactu.data,<double>p_tel.diam,
                                <double*>k2.data, <long*>NlayersDM.data,
                                <long*>indlayersDM.data,<double>FoV,<double*>pitch.data,
                                <double*>alt_DM.data)
    print "Done"
    
    print "Computing Cmm..."
    controller_mv.compute_Cmm(g_atmos.s_a, g_wfs.sensors, <double*>L0_d.data, 
                              <double*>frac_d.data, <double*>alphaX.data, <double*>alphaY.data,
                              <double>p_tel.diam, <double>p_tel.cobs)
    print "Done"
    
    Nact = np.zeros([nactu,nactu],dtype=np.float32)
    F = np.zeros([nactu,nactu],dtype=np.float32)
    ind = 0
    for k in p_rtc.controllers[ncontro].ndm:
        if(p_dm[k].type_dm == "pzt"):
            Nact[ind:ind+p_dm[k]._ntotact,ind:ind+p_dm[k]._ntotact] = create_nact_geom(p_dm,k)
            F[ind:ind+p_dm[k]._ntotact,ind:ind+p_dm[k]._ntotact] = create_piston_filter(p_dm,k)
            ind += p_dm[k]._ntotact
    
    controller_mv.filter_cphim(<float*>F.data,<float*> Nact.data)
    
cpdef selectDMforLayers(int ncontro, Param_atmos p_atmos, Param_rtc p_rtc, list p_dms):
    """ For each atmos layer, select the DM which have to handle it in the Cphim computation for MV controller
    :parameters:
        ncontro : (int) : controller number
        p_atmos : (Param_atmos) : atmos parameters
        p_rtc : (Param_rtc) : rtc parameters
        p_dms :(list of Param_dm) : dms parameters
        
    :return:
        indlayersDM : (np.array(dtype=np.int32)) : for each atmos layer, the Dm number corresponding to it
    """
    indlayersDM = np.zeros(p_atmos.nscreens,dtype=np.int64)
    for i in range(p_atmos.nscreens):
        mindif = 1e6
        for j in p_rtc.controllers[ncontro].ndm:
            alt_diff = np.abs(p_dms[j].alt - p_atmos.alt[i])
            if(alt_diff < mindif):
                indlayersDM[i] = j
                mindif = alt_diff
    
    return indlayersDM

cpdef create_nact_geom(list p_dms, int ndm):
    """ Compute the DM coupling matrix
    :param:
        p_dms : (list of Param_dm) : dms parameters
        ndm : (int) : dm number
    :return:
        Nact : (np.array(dtype=np.float64)) : the DM coupling matrix
    """
    nactu = p_dms[ndm]._ntotact
    Nact = np.zeros([nactu,nactu],dtype=np.float32)
    coupling = p_dms[ndm].coupling
    dim = p_dms[ndm]._n2 - p_dms[ndm]._n1 + 1
    
    tmpx = p_dms[ndm]._i1
    tmpy = p_dms[ndm]._j1
    offs = ((p_dms[ndm]._n2 - p_dms[ndm]._n1 + 1) - (np.max(tmpx) - np.min(tmpx)))/2 - np.min(tmpx)
    tmpx = tmpx + offs + 1
    tmpy = tmpy + offs + 1
    mask = np.zeros([dim,dim],dtype=np.float32)
    shape = np.zeros([dim,dim],dtype=np.float32)
    for i in range(len(tmpx)):
        mask[tmpx[i]][tmpy[i]] = 1
    
    mask_act = np.where(mask)

    pitch = mask_act[1][1] - mask_act[1][0]
    
    for i in range(nactu):
        shape *= 0
        #Diagonal
        shape[tmpx[i]][tmpy[i]] = 1
        # Left, right, above and under the current actuator
        shape[tmpx[i]][tmpy[i]-pitch] = coupling
        shape[tmpx[i]- pitch][tmpy[i]] = coupling
        shape[tmpx[i]][tmpy[i]+pitch] = coupling
        shape[tmpx[i]+pitch][tmpy[i]] = coupling
        # Diagonals of the current actuators 
        shape[tmpx[i]-pitch][tmpy[i]-pitch] = coupling**2
        shape[tmpx[i]- pitch][tmpy[i]+pitch] = coupling**2
        shape[tmpx[i]+pitch][tmpy[i]+pitch] = coupling**2
        shape[tmpx[i]+pitch][tmpy[i]-pitch] = coupling **2

        Nact[:,i] = shape.T[mask_act]
    
    return Nact
    
cpdef create_piston_filter(list p_dms, int ndm):
    nactu = p_dms[ndm]._ntotact    
    F = np.ones([nactu,nactu],dtype=np.float32)
    F = F * (-1.0/nactu)
    for i in range(nactu):
        F[i][i] = 1 - 1.0/nactu
    return F

    
