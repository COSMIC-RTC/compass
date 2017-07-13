include "../par.pxi"

import numpy as np
cimport numpy as np
# np.import_array()

import time

import os

import h5py
import hdf5_utils as h5u
import pandas
from subprocess import check_output
import copy as copy

from cython.operator cimport dereference as deref, preincrement as inc

from sys import stdout

import shesha

IF USE_MPI == 1:
    # from mpi4py import MPI
    from mpi4py cimport MPI
    # C-level cdef, typed, Python objects
    # from mpi4py cimport mpi_c as mpi
    from mpi4py cimport libmpi as mpi

IF USE_MPI == 2:
    cimport mpi4py.MPI as MPI
    cimport mpi4py.libmpi as mpi

#################################################
# P-Class Rtc
#################################################
cdef class Rtc:
    # child constructor must have the same prototype (same number of
    # non-optional arguments)

    def __cinit__(self, Sensors sensor=None, Target target=None, device=-1):
        cdef carma_context * context = &carma_context.instance()
        if(device == -1):
            device = context.get_activeDevice()
        context.set_activeDevice(device, 1)
        self.device = device

        self.rtc = new sutra_rtc(context)

    def __dealloc__(self):
        del self.rtc

    def add_centroider(self, Sensors sensor, long nwfs, long nvalid,
                       bytes type_centro, float offset, float scale):
        """
            Add a centroider in the sutra_centroiders vector of the RTC on the GPU

        :parameters:
            sensor: (Sensors) : sutra_sensors object (GPU)
            nwfs : (long) : number of wfs
            nvalid: (long) : number of valid subaps
            type_centro: (str) : centroider's type
            offset: (float) :
            scale: (float) :
        """

        cdef carma_context * context = &carma_context.instance()
        cdef int activeDevice = self.rtc.device
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        self.rtc.add_centroider(
            sensor.sensors,
            nwfs,
            nvalid,
            offset,
            scale,
            activeDevice,
            type_centro)

    def add_Controller(self, int nactu, float delay, bytes type_control, Dms dms,
                       list dmseen, np.ndarray[ndim=1, dtype=np.float32_t] alt,
                       int ndm, long Nphi=-1):
        """
            Add a controller in the sutra_controller vector of the RTC on the GPU

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

        type_dmseen = < char ** > malloc(len(dmseen) * sizeof(char * ))
        for j in range(controller.ndm.size):
            type_dmseen[j] = dmseen[j]

        self._add_controller(nactu, delay,
                             type_control, dms, type_dmseen,
                             alt, ndm, Nphi)

    cdef _add_controller(self, int nactu, float delay, bytes type_control, Dms dms,
                         char ** type_dmseen, np.ndarray[ndim=1, dtype=np.float32_t] alt,
                         int ndm, long Nphi=-1, bool wfs_direction=False):
        """
            Add a controller in the sutra_controller vector of the RTC on the GPU

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

        cdef float * ptr_alt = < float * > alt.data
        cdef char * ptr_dmseen = < char * > type_dmseen
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.device, 1)
        if(Nphi > -1):
            self.rtc.add_controller_geo(nactu, Nphi, delay, self.device, dms.dms, & ptr_dmseen, ptr_alt, ndm, wfs_direction)
        else:
            self.rtc.add_controller(nactu, delay, self.device, type_control, dms.dms, & ptr_dmseen, ptr_alt, ndm)

    cpdef get_pyr_method(self, int n):
        """
            Get the pyramid centroiding method

        :parameters:
            n : (int) : pyr centroider number
        """

        cdef sutra_centroider_pyr * centro = NULL

        if(self.rtc.d_centro[n].is_type(b"pyrhr")):
            centro = dynamic_cast_centroider_pyr_ptr(self.rtc.d_centro[n])
            return centro.get_method_str()
        else:
            e = "Centroider should be pyrhr, got " + \
                self.rtc.d_centro[n].get_type()
            raise ValueError(e)

    cpdef set_pyr_method(self, int n, int method, list p_centroiders):
        """
            Set the pyramid centroiding method

        :parameters:
            n : (int) : pyr centroider number
            method : (int) : new centroiding method (0: nosinus global
                                                     1: sinus global
                                                     2: nosinus local
                                                     3: sinus local)
            p_centroiders : (list of Param_centroider) :
                list of centroider parameters
        """

        cdef sutra_centroider_pyr * centro = NULL

        if(self.rtc.d_centro[n].is_type(b"pyrhr")):
            if method >= Other:
                raise ValueError("method unknown")

            centro = dynamic_cast_centroider_pyr_ptr(self.rtc.d_centro[n])
            centro.set_method(method)
            p_centroiders[n].set_method(method)
        else:
            e = "Centroider should be pyrhr, got " + \
                self.rtc.d_centro[n].get_type()
            raise ValueError(e)

    cpdef set_pyr_thresh(self, int n, float threshold, list p_centroiders):
        """
            Set the pyramid threshold
        :parameters:
            n : (int) : pyr centroider number
            threshold : (float) : new threshold in photons
            p_centroiders : (list of Param_centroider) :
                list of centroider parameters
        """

        cdef sutra_centroider_pyr * centro = NULL
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.device, 1)

        if(self.rtc.d_centro[n].is_type(b"pyrhr")):
            centro = dynamic_cast_centroider_pyr_ptr(self.rtc.d_centro[n])
            centro.set_valid_thresh(threshold)

            p_centroiders[n].set_thresh(threshold)
        else:
            e = "Centroider should be pyrhr, got " + \
                self.rtc.d_centro[n].get_type()
            raise ValueError(e)

    cpdef set_pyr_ampl(self, int n, float ampli, list p_wfss, Param_tel p_tel):
        """
            Set the pyramid modulation amplitude
        :parameters:
            n : (int) : pyr centroider number
            ampli : (float) : new amplitude in units of lambda/D
            p_wfss : (list of Param_wfs) : list of wfs parameters
            p_tel : (Param_tel) : Telescope parameters
        """

        cdef np.ndarray[ndim= 1, dtype = np.float32_t] cx
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] cy
        cdef sutra_centroider * centro = NULL
        cdef sutra_wfs_pyr * pyr = NULL
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.device, 1)

        if(self.rtc.d_centro[n].is_type(b"pyrhr")):
            centro = self.rtc.d_centro.at(n)

            nwfs = centro.nwfs
            pwfs = p_wfss[nwfs]
            pwfs.set_pyr_ampl(ampli)
            pyr = dynamic_cast_wfs_pyr_ptr(centro.wfs)

            pixsize = (np.pi * pwfs._qpixsize) / (3600 * 180)
            scale_fact = 2 * np.pi / pwfs._Nfft * \
                (pwfs.Lambda * 1e-6 / p_tel.diam) / pixsize * ampli
            cx = scale_fact * \
                np.sin((np.arange(pwfs.pyr_npts, dtype=np.float32))
                       * 2. * np.pi / pwfs.pyr_npts)
            cy = scale_fact * \
                np.cos((np.arange(pwfs.pyr_npts, dtype=np.float32))
                       * 2. * np.pi / pwfs.pyr_npts)
            pwfs.set_pyr_cx(cx)
            pwfs.set_pyr_cy(cy)
            pyr.pyr_cx.fill_from(< float * >cx.data)
            pyr.pyr_cy.fill_from(< float * >cy.data)

            centro.scale = pwfs.Lambda * 1e-6 / p_tel.diam * ampli * 180. / np.pi * \
                3600.

        else:
            e = "Centroider should be pyrhr, got " + \
                self.rtc.d_centro[n].get_type()
            raise ValueError(e)

    def rmcontrol(self):
        """
            Remove a controller
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.device, 1)
        self.dms.rm_controller()

    def setthresh(self, int ncentro, float thresh):
        """
            Set threshold for the centroider #ncentro

        :parameters:
            ncentro: (int) : centroider's index
            thresh: (float) : threshold
        """

        cdef carma_context * context = &carma_context.instance()
        cdef sutra_centroider_tcog * tcog = NULL
        cdef sutra_centroider * centro = NULL
        context.set_activeDeviceForCpy(self.device, 1)
        if(self.rtc.d_centro[ncentro].is_type(b"tcog")):
            centro = self.rtc.d_centro.at(ncentro)
            tcog = dynamic_cast_centroider_tcog_ptr(centro)
            tcog.set_threshold(thresh)

    def sensors_initbcube(self, int ncentro):
        """
            Initialize npix in the sutra_centroider_corr object (useless ?)

        :parameters:
            ncentro: (int) : centroider's index
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.device, 1)
        cdef sutra_centroider * centro = NULL
        centro = self.rtc.d_centro.at(ncentro)
        cdef sutra_centroider_corr * corr
        if(self.rtc.d_centro[ncentro].is_type(b"corr")):
            corr = dynamic_cast_centroider_corr_ptr(centro)
            corr.init_bincube()
        return 1

    def sensors_initweights(self, int ncentro, np.ndarray[ndim=2, dtype=np.float32_t] w):
        """
            Load the weight array in sutra_centroider_wcog object

        :parameters:
            ncentro: (int) : centroider's index
            w: (np.ndarray[ndim=2, dtype=np.float32_t]) : weight
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.rtc.device, 1)
        cdef sutra_centroider_wcog * centro = \
            dynamic_cast_centroider_wcog_ptr(self.rtc.d_centro[ncentro])
        centro.init_weights()
        cdef np.ndarray w_F = w.flatten("F")
        centro.load_weights( < float * > w_F.data, int(w.ndim))

    def sensors_initcorr(self, int ncentro, np.ndarray[ndim=2, dtype=np.float32_t] w,
                         np.ndarray[ndim=2, dtype=np.float32_t] corr_norm,
                         int sizex, int sizey,
                         np.ndarray[ndim=2, dtype=np.float32_t] interpmat):
        """
            Initialize sutra_centroider_corr oblect

        :parameters:
            ncentro: (int) : centroider's index
            w: (np.ndarray[ndim=1,dtype=np.float32_t]) : weight
            corr_norm: (np.ndarray[ndim=2,dtype=np.float32_t]) :
            sizex: (int) :
            sizey: (int) :
            interpmat: ([ndim=2,dtype=np.float32_t]) :
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.rtc.device, 1)

        cdef sutra_centroider_corr * centro_corr = dynamic_cast_centroider_corr_ptr(
            self.rtc.d_centro[ncentro])
        cdef np.ndarray w_F = w.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] corr_norm_F = corr_norm.flatten("F")
        cdef np.ndarray[dtype = np.float32_t] interpmat_F = interpmat.flatten("F")
        centro_corr.init_corr(sizex, sizey, < float * > interpmat_F.data)
        centro_corr.load_corr(< float * > w_F.data, < float * > corr_norm_F.data, int(w.ndim))

    # TODO possible error -> check it
    cpdef getcentroids(self, int ncontrol, Sensors g_wfs=None, int nwfs=0):
        """
            Return the centroids computed by the sutra_rtc object
            If ncontrol <= d_control.size, return rtc.d_centroids
            Else, compute centroids from wfs[nwfs] with centroider[ncontrol]

        :parameters:
            ncontrol: (int) : controller's index
            g_wfs: (Sensors) : (optional) sutra_sensors object
            nwfs: (int) : (optional) number of wfs
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : centroids (arcsec)
        """

        cdef const long * dims
        cdef sutra_wfs * wfs
        cdef sutra_rtc * rtc = self.rtc
        cdef carma_obj[float] * d_tmp
        cdef carma_obj[float] * d_data
        cdef carma_context * context = &carma_context.instance()
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data

        if(ncontrol >= abs(rtc.d_control.size())):
            if(g_wfs is None):
                raise ValueError(
                    "Controller not initialized on the GPU, you have to specify the WFS")
            wfs = g_wfs.sensors.d_wfs[nwfs]
            dims = wfs.d_slopes.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            d_tmp = new carma_obj[float](context, wfs.d_subsum.getDims())
            d_data = new carma_obj[float](context, dims)
            rtc.d_centro[ncontrol].get_cog(
                d_tmp.getData(), d_data.getData(), True)
            d_data.device2host( < float * > data.data)

        else:
            dims = rtc.d_control[ncontrol].d_centroids.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            rtc.d_control[ncontrol].d_centroids.device2host( < float * > data.data)
        return data

    cpdef docentroids(self, int ncontrol=-1):
        """
            Compute the centroids with sutra_controller #ncontrol object

        :parameters:
            ncontrol: (optional) controller's index
        """

        if(ncontrol > -1):
            self.rtc.do_centroids(ncontrol)
        else:
            self.rtc.do_centroids()

    cpdef docentroids_geom(self, int ncontrol=-1):
        """
            Compute the geometric centroids with sutra_controller #ncontrol object

        :parameters:
            ncontrol: (optional) controller's index
        """

        if(ncontrol > -1):
            self.rtc.do_centroids_geom(ncontrol)
        else:
            self.rtc.do_centroids(0)

    cpdef init_proj(self, int ncontrol, Dms dms, np.ndarray[ndim=1, dtype=np.int32_t] indx_dm,
                    np.ndarray[ndim=1, dtype=np.float32_t] unitpervolt,
                    np.ndarray[ndim=1, dtype=np.int32_t] indx_pup, np.ndarray[ndim=1, dtype=np.int32_t] indx_mpup, int roket=0):
        """
            Initialize the projection matrix for sutra_controller_geo object.
            The projection matrix is (IFt.IF)**(-1) * IFt where IF is the DMs influence functions matrix

        :parameters:
            ncontrol: (int) : controller index
            dms: (Dms) : Dms object
            indx_dm: (np.ndarray[ndim=1,dtype=np.int32_t]) : indices of where(pup) on DM screen
            unitpervolt: (np.ndarray[ndim=1,dtype=np.float32_t]) : unitpervolt DM parameter
            indx_pup: (np.ndarray[ndim=1,dtype=np.int32_t]) : indices of where(pup) on ipupil screen
            roket : (int) : optimisation flag for ROKET
        """

        cdef carma_context * context = &carma_context.instance()
        cdef sutra_controller_geo * controller_geo = \
            dynamic_cast_controller_geo_ptr(self.rtc.d_control[ncontrol])

        context.set_activeDeviceForCpy(self.rtc.device, 1)
        controller_geo.init_proj_sparse(dms.dms, < int * > indx_dm.data,
                                        < float * > unitpervolt.data, < int * > indx_pup.data, < int * > indx_mpup.data, roket)

    cpdef init_modalOpti(self, int ncontrol, int nmodes, int nrec, np.ndarray[ndim=2, dtype=np.float32_t] M2V,
                         float gmin, float gmax, int ngain, float Fs):
        """
            Initialize the modal optimization controller : compute the slopes-to-modes matrix
            and the transfer functions

        :parameters:
            ncontrol: (int) : controller index
            nmodes: (int) : number of modes
            nrec: (int) : number of recorded open slopes measurements
            M2V: (np.ndarray[ndim=2,dtype=np.float32_t]) : modes-to-volt matrix
            gmin: (float) : minimum gain for modal optimization
            gmax: (float) : maximum gain for modal optimization
            ngain: (int) : Number of tested gains
            Fs: (float) : sampling frequency
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls

        cdef np.ndarray[dtype = np.float32_t] M2V_F = M2V.flatten("F")

        if(self.rtc.d_control[ncontrol].get_type() == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.init_modalOpti(nmodes, nrec, < float * > M2V_F.data, gmin, gmax, ngain, Fs)
        else:
            raise TypeError(
                "**** ERROR : Modal Optimization only for controller type ls ****")

    cpdef loadOpenLoop(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] ol_slopes):
        """
            Load an array of recoded open-loop measurements for modal optimization

        :parameters:
            ncontrol: (int) : controller index
            ol_slopes: (np.ndarray[ndim=2, dtype=np.float32_t]) : open-loop slopes
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls

        cdef np.ndarray[dtype = np.float32_t] slopes_F = ol_slopes.flatten("F")

        if(self.rtc.d_control[ncontrol].get_type() == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.loadOpenLoopSlp( < float * > slopes_F.data)
        else:
            raise TypeError("Controller type must be ls")

    cpdef modalControlOptimization(self, int ncontrol):
        """
            Compute the command matrix with modal control optimization

        :parameter:
            ncontrol: controller index
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls

        if(self.rtc.d_control[ncontrol].get_type() == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            if(controller_ls.is_modopti):
                controller_ls.modalControlOptimization()
            else:
                raise ValueError(
                    "**** ERROR : Modal Optimization not initialized ****")
        else:
            raise TypeError(
                "**** ERROR : Modal Optimization only for controller type ls ***")

    cpdef set_gain(self, int ncontrol, float gain):
        """
            Set the loop gain in sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            gain: (float) : loop gain
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_cured * controller_cured
        cdef sutra_controller_geo * controller_geo
        cdef sutra_controller_kalman * controller_kl

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.set_gain(gain)

        elif(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.set_gain(gain)

        elif(type_control == b"cured"):
            controller_cured = dynamic_cast_controller_cured_ptr(
                self.rtc.d_control[ncontrol])
            controller_cured.set_gain(gain)

        elif(type_control == b"geo"):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            controller_geo.set_gain(gain)
        elif(type_control == b"kalman_CPU" or type_control == b"kalman_GPU" or
             type_control == b"kalman_uninitialized"):
            controller_kl = dynamic_cast_controller_kl_ptr(
                self.rtc.d_control[ncontrol])
            controller_kl.set_gain(gain)
        else:
            raise TypeError(
                "Controller needs to be ls, mv, cured, geo, kalman_GPU or kalman_CPU. for generic (g=1.0) use set_mgain")

    cpdef set_mgain(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] mgain):
        """
            Set modal gains in sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            mgain: (np.ndarray[ndim=1,dtype=np.float32_t]) : modal gains
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.set_mgain( < float * > mgain.data)
        elif(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.set_mgain( < float * > mgain.data)
        elif(type_control == b"generic"):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.set_mgain( < float * > mgain.data)
        else:
            raise TypeError("Controller needs to be ls, mv or generic")

    def set_scalar_mgain(self, ncontrol, g):
        x = self.get_mgain(ncontrol)
        x[:] = g
        self.set_mgain(ncontrol, x.copy())

    cpdef setCom(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] comvec):
        """
            Set the command vector of a sutra_controller object to comvec

        :parameters:
            ncontrol: (int) : controller index
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        self.rtc.d_control[ncontrol].d_com.host2device( < float * > comvec.data)

    cpdef setCentroids(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] centro):
        """
            Set the centroids vector of a sutra_controller object to centro

        :parameters:
            ncontrol: (int) : controller index
            centro: (np.ndarray[ndim=1,dtype=np.float32_t]) : centroids vector
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        self.rtc.d_control[ncontrol].d_centroids.host2device( < float * > centro.data)

    cpdef get_gain(self, int ncontrol):
        """
            Return modal gains from sutra_controller

        :parameters:
            ncontrol: (int) : controller index
        :return:
            mgain : (np.ndarray[ndim=1,dtype=np.float32_t]) : modal gains
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_cured * controller_cured
        cdef sutra_controller_geo * controller_geo

        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
#        cdef float gain =  self.rtc.d_control[ncontrol].get_gain()
#        return gain
        cdef float gain
        raise TypeError("get_gain are not implemented")
    '''
        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(self.rtc.d_control[ncontrol])
            gain = np.float32(0)
            controller_ls.gain.device2host(gain)
            return gain

        elif(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(self.rtc.d_control[ncontrol])
            gain = np.float32(0)
            controller_mv.gain.device2host(gain)
            return gain

        elif(type_control == b"generic"):
            gain = np.float32(1.0)
            return gain

        elif(type_control == b"cured"):
            controller_cured = dynamic_cast_controller_cured_ptr(self.rtc.d_control[ncontrol])
            gain = np.float32(0)
            controller_cured.gain.device2host(gain)
            return gain

        elif(type_control == b"geo"):
            controller_geo = dynamic_cast_controller_geo_ptr(self.rtc.d_control[ncontrol])
            gain = np.float32(0)
            controller_geo.gain.device2host(gain)
            return gain

        else:
            raise TypeError("Controller needs to be ls, generic, mv, cured or geo")
    '''

    cpdef get_mgain(self, int ncontrol):
        """
            Return modal gains from sutra_controller

        :parameters:
            ncontrol: (int) : controller index
        :return:
            mgain : (np.ndarray[ndim=1,dtype=np.float32_t]) : modal gains
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef int size
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] mgain

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            size = controller_ls.d_gain.getNbElem()
            mgain = np.zeros(size, dtype=np.float32)
            controller_ls.d_gain.device2host( < float * > mgain.data)
            return mgain
        elif(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            size = controller_mv.d_gain.getNbElem()
            mgain = np.zeros(size, dtype=np.float32)
            controller_mv.d_gain.device2host( < float * > mgain.data)
            return mgain
        elif(type_control == b"generic"):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            size = controller_generic.d_gain.getNbElem()
            mgain = np.zeros(size, dtype=np.float32)
            controller_generic.d_gain.device2host( < float * > mgain.data)
            return mgain
        else:
            raise TypeError("Controller needs to be ls, generic or mv")

    cpdef set_imat(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            Set the interaction matrix on a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            data: (np.ndarray[ndim=2,dtype=np.float32_t]) : interaction matrix to use
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef np.ndarray[dtype = np.float32_t] data_F = data.flatten("F")

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.d_imat.host2device( < float * > data_F.data)
        elif(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.d_imat.host2device( < float * > data_F.data)
        else:
            raise TypeError("Controller needs to be ls or mv")

    cpdef get_imat(self, int ncontrol):
        """
            Return the interaction matrix of a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            imat : (np.ndarray[ndim=2,dtype=np.float32_t]) : interaction matrix
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_cured * controller_cured
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        cdef const long * dims = NULL
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] imat_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] imat

        if(type_control == b"generic" or type_control == b"geo"):
            raise TypeError("Generic controller doesn't have imat")

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.d_imat.getDims()
            imat_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_ls.d_imat.device2host( < float * > imat_F.data)

        elif(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_mv.d_imat.getDims()
            imat_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_mv.d_imat.device2host( < float * > imat_F.data)

        elif(type_control == b"cured"):
            controller_cured = dynamic_cast_controller_cured_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_cured.d_imat.getDims()
            imat_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_cured.d_imat.device2host( < float * > imat_F.data)

        imat = np.reshape(imat_F.flatten("F"), (dims[1], dims[2]))
        return imat

    cpdef set_cmat(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            Set the command matrix on a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            data: (np.ndarray[ndim=2,dtype=np.float32_t]) : command matrix to use
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef np.ndarray[dtype = np.float32_t] data_F = data.flatten("F")

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.d_cmat.host2device( < float * > data_F.data)
        elif(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.d_cmat.host2device( < float * > data_F.data)
        elif(type_control == b"generic"):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.d_cmat.host2device( < float * > data_F.data)
        else:
            raise TypeError("Controller needs to be ls, mv or generic")

    cpdef set_cmm(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            Set the Cmm matrix on a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index
            data: (np.ndarray[ndim=2,dtype=np.float32_t]) : Cmm matrix
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef np.ndarray[dtype = np.float32_t] data_F = data.flatten("F")

        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if type_control == b"mv":
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.d_Cmm.host2device( < float * > data_F.data)
        else:
            raise TypeError("Controller needs to be mv")

    cpdef get_cmm(self, int ncontrol):
        """
            Return the Cmm matrix from a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            cmm : (np.ndarray[ndim=2,dtype=np.float32_t]) : Cmm matrix
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F, cmm
        cdef const long * cdims

        if(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_mv.d_Cmm.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_mv.d_Cmm.device2host( < float * > data_F.data)
            cmm = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return cmm
        else:
            raise TypeError("Controller needs to be mv")

    cpdef get_cphim(self, int ncontrol):
        """
            Return the Cphim matrix from a sutra_controller_mv object

        :parameters;
            ncontrol: (int) : controller index

        :return:
            cphim : (np.ndarray[ndim=2,dtype=np.float32_t]) : Cphim matrix
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F, cphim
        cdef const long * cdims

        if(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_mv.d_Cphim.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_mv.d_Cphim.device2host( < float * > data_F.data)
            cphim = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return cphim
        else:
            raise TypeError("Controller needs to be mv")

    cpdef get_cmat(self, int ncontrol):
        """
            Return the command matrix from a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            cmat : (np.ndarray[ndim=2,dtype=np.float32_t]) : command matrix
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F, cmat
        cdef const long * cdims

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_ls.d_cmat.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_ls.d_cmat.device2host( < float * > data_F.data)
            cmat = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return cmat
        elif(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_mv.d_cmat.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_mv.d_cmat.device2host( < float * > data_F.data)
            cmat = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return cmat
        elif(type_control == b"generic"):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_generic.d_cmat.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_generic.d_cmat.device2host( < float * > data_F.data)
            cmat = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return cmat
        else:
            raise TypeError("Controller needs to be ls, mv or generic")

    cpdef set_decayFactor(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] decay):
        """
            Set the decay factor on a sutra_controller_generic object

        :parameters:
            ncontrol: (int) : controller index
            decay: (np.ndarray[ndim=1,dtype=np.float32_t]) : ask to Rico
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == b"generic"):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.set_decayFactor( < float * > decay.data)
        else:
            raise TypeError("Controller needs to be generic")

    cpdef set_matE(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] matE):
        """
            Set the matrix E on a sutra_controller_generic object

        :parameters:
            ncontrol: (int) : controller index
            matE: (np.ndarray[ndim=2,dtype=np.float32_t]) : ask to Rico
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        cdef np.ndarray[dtype = np.float32_t] matE_F = matE.flatten("F")

        if(type_control == b"generic"):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.set_matE( < float * > matE_F.data)
        else:
            raise TypeError("Controller needs to be generic")

    cpdef set_openloop(self, int ncontrol, int openloop):
        """
            Set the openloop state to a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            openloop: state of the controller
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller * controller = self.rtc.d_control[ncontrol]
        controller.set_openloop(openloop)

    cpdef set_commandlaw(self, int ncontrol, bytes law):
        """
            Set the law to a sutra_controller_generic object

        :parameters:
            ncontrol: (int) : controller index
            law: law of the controller
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == b"generic"):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.set_commandlaw(law)
        else:
            raise TypeError("Controller needs to be ls, mv or generic")

    cpdef doimat_geom(self, int ncontrol, Dms g_dms, int geom):
        """
            Compute the interaction matrix by using a geometric centroiding method

        :parameters:
            ncontrol: (int) : controller index
            g_dms: (Dms) : Dms object
            geom: (int) : type of geometric method (0 or 1)
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        self.rtc.do_imat_geom(ncontrol, g_dms.dms, geom)
        print("TODO call imat_geom")
        print("TODO set_imat")

    cpdef set_perturbcom(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] perturb):
        """
        :parameters:
            ncontrol: (int) : controller index
            perturb: (np.ndarray[ndim = 2, dtype = np.float32_t]) : perturb voltages
        """

        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        control.set_perturbcom( < float * > perturb.data, perturb.shape[0])

    cpdef set_centroids_ref(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] centroids_ref):
        """
        :parameters:
            ncontrol: (int) : controller index
        """

        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        control.set_centroids_ref( < float * > centroids_ref.data)

    cpdef get_centroids_ref(self, int ncontrol):
        """
        :parameters:
            ncontrol: (int) : controller index
        """

        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        cdef int nslope = control.nslope()
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] centroids_ref
        centroids_ref = np.zeros(nslope, dtype=np.float32)
        control.get_centroids_ref( < float * > centroids_ref.data)
        return centroids_ref

    cpdef do_centroids_ref(self, int ncontrol):
        '''
            This should be a docstring.
        '''

        cdef carma_obj[float] * phase
        cdef sutra_wfs * wfs
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] h_ref
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] h_rawslp
        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        cdef int nslope = control.nslope()
        cdef float tmp

        print("Doing refslp...")
        for idx_cntr in range(< int > self.rtc.d_centro.size()):
            wfs = self.rtc.d_centro[idx_cntr].wfs
            phase = wfs.d_gs.d_phase.d_screen
            phase.reset()
            tmp = wfs.noise
            wfs.noise = -1
            wfs.comp_image()
            wfs.noise = tmp
        self.rtc.do_centroids(ncontrol)
        h_ref = np.zeros(nslope, dtype=np.float32)
        self.rtc.get_centroids_ref(ncontrol, < float * > h_ref.data)
        h_rawslp = self.getCentroids(ncontrol) + h_ref
        self.rtc.set_centroids_ref(ncontrol, < float * > h_rawslp.data)

    cpdef doimat_kl(self, int ncontrol, Param_controller controller, Dms g_dms, p_dms, np.ndarray[ndim=2, dtype=np.float32_t] kl):
        '''
            Changer de strategie, on enregister l'imat sur python puis self.set_imat()
            puis appelle d'un fonction set Cmat modifdans le rtc_init à la place de la fonction set normal
        '''

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        cdef carma_obj[float] * d_imat = NULL
        cdef carma_obj[float] * d_imat_ret = NULL
        cdef long dims[3]
        cdef int nactu = control.nactu()
        cdef int nslope = control.nslope()

        if(control.get_type() == b"ls"):
            d_imat = dynamic_cast_controller_ls_ptr(control).d_imat
        elif(control.get_type() == b"mv"):
            d_imat = dynamic_cast_controller_mv_ptr(control).d_imat
        else:
            # Create a temporary imat to return

            dims[0] = 2
            dims[1] = nactu
            dims[2] = nslope
            d_imat_ret = new carma_obj[float](context, dims)
            d_imat = d_imat_ret
            print(
                "WARNING: the controller is NOT a LS or MV, the imat computed will be returned")

        cdef int inds1, j, idx_cntr, device
        cdef float tmp_noise
        cdef float tmp_nphot
        inds1 = 0
        cdef sutra_dm * dm
        cdef sutra_wfs * wfs
        cdef carma_obj[float] * screen
        cdef vector[sutra_dm * ].iterator it_dm

        cdef float * d_centroids
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] h_centroids

        cdef int rank
        IF USE_MPI:
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, & rank)
        ELSE:
            rank = 0

        cdef carma_obj[float] * phase

        klNorm = kl * 0.
        for k in range(kl.shape[1]):
            klmaxVal = np.abs(kl[:, k]).max()
            # Norm du max des modes a 1
            klNorm[:, k] = kl[:, k] / klmaxVal

        pp = 0
        ppz = np.zeros(np.size(p_dms), dtype=np.int32)  # boolean vector pzt dm
        ntt = 0
        ttz = 0  # tt dm position
        na = np.zeros(np.size(p_dms), dtype=np.int32)  # number actu for pzt_dm

        for i in range(np.size(p_dms)):
            if (p_dms[i].type_dm == b"pzt"):
                ppz[i] = 1
                pp += 1
                na[i] = p_dms[i]._ntotact
            elif ((p_dms[i].type_dm == b"tt") & (ntt == 0)):
                ttz = i
                ntt += 1
            else:
                raise type("need one 0 or 1 tt dm")
        cc = 0

        it_dm = control.d_dmseen.begin()  # ---> a changer
        print("Doing imat_kl...%d%%\r" % cc, end='')
        dm = deref(it_dm)

        pushkl = np.ones(kl.shape[0])
        a = 0
        for i in range(np.sum(ppz)):
            pushkl[a:a + na[i]] = p_dms[i].push4imat
            a = na[i]
        if ntt != 0:
            pushkl[-1] = p_dms[ttz].push4imat
            pushkl[-2] = p_dms[ttz].push4imat

        for j in range(kl.shape[1]):

            g_dms.set_full_comm(np.float32(klNorm[:, j].copy() * pushkl))

            for idx_cntr in range(< int > self.rtc.d_centro.size()):
                wfs = self.rtc.d_centro[idx_cntr].wfs
                screen = wfs.d_gs.d_phase.d_screen
                tmp_noise = wfs.noise
                tmp_nphot = wfs.nphot
                wfs.nphot = wfs.nphot4imat
                wfs.noise = -1
                wfs.kernconv = True
                wfs.sensor_trace(g_dms.dms, 1)
                IF USE_MPI:
                    Bcast(screen, 0)
                wfs.comp_image()
                wfs.noise = tmp_noise
                wfs.nphot = tmp_nphot
                wfs.kernconv = False

            self.rtc.do_centroids(ncontrol, True)

            h_centroids = self.getCentroids(ncontrol)
            control.d_centroids.host2device( < float * > h_centroids.data)

            device = control.d_centroids.getDevice()
            d_centroids = control.d_centroids.getData()
            if ((ntt != 0)and(j >= kl.shape[1] - 2)):
                convert_centro(d_centroids,
                               d_centroids, 0,
                               (1. / p_dms[ttz].push4imat),
                               control.d_centroids.getNbElem(),
                               context.get_device(device))
            else:
                convert_centro(d_centroids,
                               d_centroids, 0,
                               (1. / dm.push4imat),
                               control.d_centroids.getNbElem(),
                               context.get_device(device))

            control.d_centroids.copyInto(
                d_imat.getDataAt(inds1),
                control.nslope())
            for i in range(np.size(p_dms)):
                g_dms.resetdm(p_dms[i].type_dm, p_dms[i].alt)

            inds1 += control.nslope()
            cc = cc + 1
            print("Doing imat... #%d/%d \r" % (cc, kl.shape[1]), end=' ')

            #cont +=1
            # inc(it_dm)
        print("imat done")

        cdef np.ndarray[ndim = 2, dtype = np.float32_t] h_imat_ret
        if d_imat_ret != NULL:
            h_imat_ret = np.zeros((nactu, nslope), dtype=np.float32)
            d_imat_ret.device2host( < float * > h_imat_ret.data)
            del d_imat_ret
            return h_imat_ret

    cpdef doimat(self, int ncontrol, Dms g_dms):
        """
            Compute the interaction matrix

        :parameters:
            ncontrol: (int) : controller index
            g_dms: (Dms) : Dms object
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        cdef carma_obj[float] * d_imat = NULL
        cdef carma_obj[float] * d_imat_ret = NULL
        cdef long dims[3]
        cdef int nactu = control.nactu()
        cdef int nslope = control.nslope()

        if(control.get_type() == b"ls"):
            d_imat = dynamic_cast_controller_ls_ptr(control).d_imat
        elif(control.get_type() == b"mv"):
            d_imat = dynamic_cast_controller_mv_ptr(control).d_imat
        else:
            # Create a temporary imat to return
            dims[0] = 2
            dims[1] = nactu
            dims[2] = nslope
            d_imat_ret = new carma_obj[float](context, dims)
            d_imat = d_imat_ret
            print(
                "WARNING: the controller is NOT a LS or MV, the imat computed will be returned")

        cdef int inds1, j, idx_cntr, device
        cdef float tmp_noise
        cdef float tmp_nphot
        inds1 = 0
        cdef sutra_dm * dm
        cdef sutra_wfs * wfs
        cdef carma_obj[float] * screen
        cdef vector[sutra_dm * ].iterator it_dm

        cdef float * d_centroids
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] h_centroids

        cdef int rank
        IF USE_MPI:
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, & rank)
        ELSE:
            rank = 0

        cdef carma_obj[float] * phase

        cc = 0
        it_dm = control.d_dmseen.begin()
        print("Doing imat...%d%% \r" % cc, end='')
        while(it_dm != control.d_dmseen.end()):
            dm = deref(it_dm)
            for j in range(dm.ninflu):

                # push actu __________________________________________________

                dm.comp_oneactu(j, dm.push4imat)
                for idx_cntr in range(< int > self.rtc.d_centro.size()):
                    wfs = self.rtc.d_centro[idx_cntr].wfs
                    screen = wfs.d_gs.d_phase.d_screen
                    tmp_noise = wfs.noise
                    tmp_nphot = wfs.nphot
                    wfs.nphot = wfs.nphot4imat
                    wfs.noise = -1
                    wfs.kernconv = True
                    wfs.sensor_trace(g_dms.dms, 1)
                    IF USE_MPI:
                        Bcast(screen, 0)
                    wfs.comp_image()
                    wfs.noise = tmp_noise
                    wfs.nphot = tmp_nphot
                    wfs.kernconv = False

                self.rtc.do_centroids(ncontrol, True)

                h_centroids = self.getCentroids(ncontrol)
                control.d_centroids.host2device( < float * > h_centroids.data)

                device = control.d_centroids.getDevice()
                d_centroids = control.d_centroids.getData()

                convert_centro(d_centroids,
                               d_centroids, 0,
                               (0.5 / dm.push4imat),
                               control.d_centroids.getNbElem(),
                               context.get_device(device))

                control.d_centroids.copyInto(
                    d_imat.getDataAt(inds1),
                    control.nslope())
                dm.reset_shape()

                # pul actu __________________________________________________
                dm.comp_oneactu(j, -1. * dm.push4imat)

                for idx_cntr in range(< int > self.rtc.d_centro.size()):
                    wfs = self.rtc.d_centro[idx_cntr].wfs
                    tmp_noise = wfs.noise
                    tmp_nphot = wfs.nphot
                    wfs.nphot = wfs.nphot4imat
                    wfs.noise = -1
                    wfs.kernconv = True
                    wfs.sensor_trace(g_dms.dms, 1)
                    wfs.comp_image()
                    wfs.noise = tmp_noise
                    wfs.nphot = tmp_nphot
                    wfs.kernconv = False

                self.rtc.do_centroids(ncontrol, True)

                h_centroids = self.getCentroids(ncontrol)
                control.d_centroids.host2device( < float * > h_centroids.data)

                device = control.d_centroids.getDevice()
                d_centroids = control.d_centroids.getData()

                convert_centro(d_centroids,
                               d_centroids, 0,
                               (0.5 / dm.push4imat),
                               control.d_centroids.getNbElem(),
                               context.get_device(device))

                carma_axpy[float](context.get_cublasHandle(),
                                  control.d_centroids.getNbElem(), < float > -1.0,
                                  d_centroids, 1,
                                  d_imat.getDataAt(inds1), 1)

                dm.reset_shape()

                inds1 += control.nslope()
                cc = cc + 1
                print("Doing imat... #%d/%d \r" % (cc, nactu), end=' ')

            inc(it_dm)
        print("imat done\n")

        cdef np.ndarray[ndim = 2, dtype = np.float32_t] h_imat_ret
        if d_imat_ret != NULL:
            h_imat_ret = np.zeros((nactu, nslope), dtype=np.float32)
            d_imat_ret.device2host( < float * > h_imat_ret.data)
            del d_imat_ret
            return h_imat_ret

    cpdef sensors_compslopes(self, int ncentro, int nmax=-1, float thresh=-1):
        """
            Compute the slopes in a sutra_wfs object. This function is equivalent to
            docentroids() but the centroids are stored in the sutra_wfs object instead of
            the sutra_rtc object

        :parameters:
            ncentro: (int) : centroider index
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        self.rtc.d_centro[ncentro].get_cog()

    cpdef imat_svd(self, int ncontrol):
        """
            Compute the singular value decomposition of the interaction matrix

        :param ncontrol: controller index
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef sutra_controller_ls * controller_ls
        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])

            if(controller_ls.svdec_imat() == 1):
                raise RuntimeError(
                    "sutra controller has no SVD implementation")

        else:
            raise TypeError("Controller needs to be ls")

    cpdef setU(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] U):
        """
            Set the eigen modes matrix of the imat decomposition in a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index
            U: (np.ndarray[ndim=2,dtype=np.float32_t]) : eigen modes matrix
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        cdef np.ndarray[dtype = np.float32_t] data_F = U.flatten("F")
        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.d_U.host2device( < float * > data_F.data)

    cpdef setEigenvals(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] eigenvals):
        """
            Set the eigen values of the imat decomposition in a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index
            eigenvals: (np.ndarray[ndim=1,dtype=np.float32_t]) : eigen values
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.h_eigenvals.fill_from( < float * > eigenvals.data)
        else:
            raise TypeError("Controller needs to be ls")

    cpdef getU(self, int ncontrol):
        """
            Return the eigen modes matrix of the imat decomposition from a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            U : (np.ndarray[ndim=2,dtype=np.float32_t]) : eigen modes matrix
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.d_U.getDims()
            data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_ls.d_U.device2host( < float * > data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    cpdef getEigenvals(self, int ncontrol):
        """
            Return the eigen values of the imat decomposition in a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            eigenvals : (np.ndarray[ndim=1,dtype=np.float32_t]) : eigenvalues
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims
        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.h_eigenvals.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            controller_ls.h_eigenvals.fill_into( < float * > data.data)
        if(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_mv.h_Cmmeigenvals.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            controller_mv.h_eigenvals.fill_into( < float * > data.data)

        return data

    cpdef getCmmEigenvals(self, int ncontrol):
        """
            Return the eigen values of the Cmm decomposition in a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            eigenvals : (np.ndarray[ndim=1,dtype=np.float32_t]) : eigenvalues
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims
        if(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_mv.h_Cmmeigenvals.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            controller_mv.h_Cmmeigenvals.fill_into( < float * > data.data)

        return data

    cpdef getGeocov(self, int ncontrol):
        """
            Return the geocov matrix of the sutra_controller_geo object. In case of error_budget computation, this matrix is Btt basis
        :parameters:
            ncontrol: (int) : controller index
        :return:
            geocov : (np.ndarray[ndim=2,dtype=np.float32_t]) : geocov matrix
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_contro = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef const long * dims
        if(type_contro == b"geo"):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_geo.d_geocov.getDims()
            data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_geo.d_geocov.device2host( < float*>data_F.data)

            data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        else:
            raise 'controller must be geo'

        return data

    cpdef get_IFtt(self, int ncontrol):
        """
            Return the TT IF matrix of the sutra_controller_geo object (in case of error_budget computation)
        :parameters:
            ncontrol: (int) : controller index
        :return:
            geocov : (np.ndarray[ndim=2,dtype=np.float32_t]) : IF TT matrix
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_contro = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef const long * dims
        if(type_contro == b"geo"):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            if(controller_geo.Ntt):
                dims = controller_geo.d_TT.getDims()
                data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
                controller_geo.d_TT.device2host( < float*>data_F.data)

                data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))

                return data
            else:
                raise ValueError("TT not initialized : only with roket")

    cpdef getCenbuff(self, int ncontrol):
        """
            Return the centroids buffer from a sutra_controller_ls object.
            This buffer contains centroids from iteration i-delay to current iteration.

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : centroids buffer
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef const long * dims = NULL

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.d_cenbuff.getDims()
            data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_ls.d_cenbuff.device2host( < float * > data_F.data)

        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    cpdef getErr(self, int ncontrol):
        """
            Return the command increment (cmat*slopes) from a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : command increment
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims
        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.d_err.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            controller_ls.d_err.device2host( < float * > data.data)

        return data

    cpdef getCom(self, int ncontrol):
        """
            Return the command vector from a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : command vector
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims
        dims = self.rtc.d_control[ncontrol].d_com.getDims()
        data = np.zeros((dims[1]), dtype=np.float32)
        self.rtc.d_control[ncontrol].d_com.device2host( < float * > data.data)

        return data

    cpdef getolmeas(self, int ncontrol):
        """
            Return the reconstructed open-loop measurement from a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : reconstructed open-loop
        """

        cdef carma_context * context = &carma_context.instance()
        cdef sutra_controller_mv * controller_mv

        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims

        controller_mv = dynamic_cast_controller_mv_ptr(
            self.rtc.d_control[ncontrol])
        dims = controller_mv.d_olmeas.getDims()
        data = np.zeros((dims[1]), dtype=np.float32)
        controller_mv.d_olmeas.device2host( < float * > data.data)

        return data

    cpdef getVoltage(self, int ncontrol):
        """
            Return the voltage vector that will be effectively applied to the DMs

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : voltage vector
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims
        dims = self.rtc.d_control[ncontrol].d_voltage.getDims()
        data = np.zeros((dims[1]), dtype=np.float32)
        self.rtc.d_control[ncontrol].d_voltage.device2host( < float * > data.data)

        return data

    cpdef getCentroids(self, int ncontrol):
        """
            Return the centroids computed by the sutra_rtc object

        :parameters:
            ncontrol: (int) : controller's index

        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : centroids (arcsec)
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        cdef const long * dims
        dims = self.rtc.d_control[ncontrol].d_centroids.getDims()
        data = np.zeros((dims[1]), dtype=np.float32)
        self.rtc.d_control[ncontrol].d_centroids.device2host( < float * > data.data)

        IF USE_MPI:
            cdef np.ndarray[ndim = 1, dtype = np.float32_t] all_centroids
            cdef int comm_size, rank
            mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD, & comm_size)
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, & rank)

            cdef int * count = < int * > malloc(comm_size * sizeof(int))
            cdef int * disp = < int * > malloc((comm_size + 1) * sizeof(int))
            cdef int d = dims[1] / (comm_size * 2)
            cdef int i

            disp[0] = 0
            for i in range(comm_size):
                if(i < (dims[1] / 2) % comm_size):
                    count[i] = d + 1
                else:
                    count[i] = d

                disp[i + 1] = disp[i] + count[i]

            all_centroids = np.zeros(disp[comm_size] * 2, dtype=np.float32)

            cdef float * send = < float * > data.data
            cdef float * recv = < float * > all_centroids.data

            # gather centroids X axis
            mpi.MPI_Allgatherv(send, count[rank], mpi.MPI_FLOAT,
                               recv, count, disp, mpi.MPI_FLOAT,
                               mpi.MPI_COMM_WORLD)

            # gather centroids Y axis
            mpi.MPI_Allgatherv(& send[disp[comm_size]], count[rank], mpi.MPI_FLOAT,
                                & recv[disp[comm_size]], count, disp,
                                mpi.MPI_FLOAT, mpi.MPI_COMM_WORLD)

            free(count)
            free(disp)
            return all_centroids

        ELSE:
            return data

    def buildcmat_kl(self, ncontrol, KL2V, imat, gains):
        """
            Compute the command matrix in a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index

            nmode: (int) : number of modes to filter
        """

        type_control = self.rtc.d_control[ncontrol].get_type()

        print(" nmode:")
        print(KL2V.shape[1])
        print("imat size:")
        print(imat.shape)
        print("KL2V size :")
        print(KL2V.shape)

        if(type_control == b"ls"):

            KL2VN = KL2V * 0.
            for k in range(KL2V.shape[1]):
                # Norm du max des modes a 1
                KL2VN[:, k] = ((KL2V[:, k]) / (np.abs(KL2V[:, k]).max()))
            # filter imat
            D_filt = imat[:, :KL2V.shape[1]]

            # Direct inversion
            Dp_filt = np.linalg.inv(D_filt.T.dot(D_filt)).dot(D_filt.T)
            if (gains is not None):
                print("shape gains :")
                print(gains.shape[0])
                if (gains.shape[0] == KL2V.shape[1]):
                    for i in range(KL2V.shape[1]):
                        Dp_filt[:, i] *= gains[i]
                else:
                    print("Need size :")
                    print(KL2V.shape[1])
                    raise TypeError("incorect size for klgain vector")

            # Command matrix
            cmat_filt = KL2VN.dot(Dp_filt)

            self.set_cmat(0, cmat_filt.astype(np.float32).copy())
            print("cmat kl end processing")

        else:
            raise TypeError("Controller needs to be ls for kl imat")

    cpdef buildcmat(self, int ncontrol, int nfilt, int filt_tt=0):
        """
            Compute the command matrix in a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index

            nfilt: (int) : number of modes to filter

            filt_tt: (int) : (optional) flag to filter TT
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == b"ls"):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            if(filt_tt > 0):
                controller_ls.build_cmat(nfilt, True)
            else:
                controller_ls.build_cmat(nfilt)

    cpdef buildcmatmv(self, int ncontrol, float cond):
        """
            Compute the command matrix in a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index

            cond: (float) : conditioning factor for the Cmm inversion
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.build_cmat(cond)

    cpdef loadnoisemat(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] N):
        """
            Load the noise vector on a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index

            N: (np.ndarray[ndim=1,dtype=np.float32_t]) : noise vector
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == b"mv"):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.load_noisemat( < float * > N.data)

    cpdef doclipping(self, int ncontrol, float min, float max):
        """
            Clip the command to apply on the DMs on a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            min: (float) : minimum value for the command
            max: (float) : maximum value for the command
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.rtc.device, 1)
        self.rtc.do_clipping(ncontrol, min, max)

    cpdef docontrol(self, int ncontrol):
        """
        Compute the command to apply on the DMs on a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.rtc.device, 1)

        self.rtc.do_control(ncontrol)

    cpdef docontrol_geo(self, int ncontrol, Dms dms, Target target, int ntarget):
        """
            Compute the command to apply on the DMs on a sutra_controller_geo object
            for the target direction

        :parameters:
            ncontrol: (int) : controller index
            dms: (shesha_dm object) : shesha_dm
            target : (shesha_target) : shesha_target
            ntarget : (int) : target number
        """

        cdef carma_context * context = &carma_context.instance()
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        context.set_activeDevice(self.rtc.device, 1)

        if(type_control == b"geo"):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            controller_geo.comp_dphi(target.target.d_targets[ntarget], False)
            self.rtc.do_control(ncontrol)
        else:
            raise TypeError("Controller needs to be geo")

    cpdef docontrol_geo_onwfs(self, int ncontrol, Dms dms, Sensors sensors, int nwfs):
        """
            Compute the command to apply on the DMs on a sutra_controller_geo object
            for the wfs direction

        :parameters:
            ncontrol: (int) : controller index
            dms: (shesha_dm object) : shesha_dm
            sensors : (shesha_sensor) : shesha_sensor
            nwfs : (int) : wfs number
        """

        cdef carma_context * context = &carma_context.instance()
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        context.set_activeDevice(self.rtc.device, 1)

        if(type_control == b"geo"):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            controller_geo.comp_dphi(sensors.sensors.d_wfs[nwfs].d_gs, True)
            self.rtc.do_control(ncontrol)
        else:
            raise TypeError("Controller needs to be geo")

    cpdef load_Btt(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] Btt):
        """
            Load the Btt basis for sutra_controller_geo projection in case of error_budget

        :parameters:
            ncontrol: (int) : controller index
            Btt: (np.ndarray[ndim=2,dtype=np.float32_t]) : Btt basis
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_contro = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray TT = Btt[-2:, -2:].copy()
        cdef np.ndarray Btt_F = Btt[:-2, :-2].copy().flatten('F')
        cdef np.ndarray TT_F = TT.flatten('F')

        if(type_contro == b"geo"):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            controller_geo.load_Btt(< float*>Btt_F.data, < float*>TT_F.data)
        else:
            raise TypeError("Controller needs to be geo")

    cpdef applycontrol(self, int ncontrol, Dms dms):
        """
            Compute the DMs shapes from the commands computed in a sutra_controller_object.
            From the command vector, it computes the voltage command (adding perturbation voltages,
            taking delay into account) and then apply it to the dms

        :parameters:
            ncontrol: (int) : controller index
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.rtc.device, 1)
        self.rtc.apply_control(ncontrol, dms.dms)

    cpdef get_nfiltered(self, int ncontrol, Param_rtc p_rtc):
        """
            Get the number of filtered modes for cmat computation

        :parameters:
            ncontrol: (int) : controller index
            p_rtc: (Param_rtc) : rtc parameters
        """

        eigenv = self.getEigenvals(ncontrol)
        maxcond = p_rtc.controllers[ncontrol].maxcond
        if(eigenv[0] < eigenv[eigenv.shape[0] - 1]):
            mfilt = np.where(
                (eigenv / eigenv[eigenv.shape[0] - 3]) < 1. / maxcond)[0]
        else:
            mfilt = np.where((1. / (eigenv / eigenv[2])) > maxcond)[0]
        nfilt = mfilt.shape[0]

        return nfilt

    def get_IFsparse(self, int ncontrol):
        """
            Get the influence functions matrix computed by the geo controller
            Return a scipy.sparse object which shape is (nactus,Npts in the pupil)

        :parameters:
            ncontrol: (int) : controller index
        :return:
            IF : (scipy.sparse) : influence functions matrix
        """

        cdef sutra_controller_geo * controller_geo
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        sparse = naga_sparse_obj_Double()
        if(type_control == b"geo"):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            sparse.copy(controller_geo.d_IFsparse)
            return sparse.get_sparse()

        else:
            print(
                "This controller does not have a sparse_obj\n(or this function is not implemented for it yet")

    def __str__(self):
        print("RTC object:")

        cdef sutra_centroider * centro
        cdef sutra_controller * contro
        cdef int i

        info = "Contains " + str(self.rtc.d_centro.size()) + " Centroider(s)\n"
        info += "Centro # | Type  | nwfs | Nvalid\n"

        for i in range(< int > self.rtc.d_centro.size()):
            centro = self.rtc.d_centro[i]
            info += "%8d" % (i + 1) + " | " + "%5s" % centro.get_type() + " | " + "%4d" % (centro.nwfs + 1) + \
                    " | " + str(centro.nvalid) + "\n"

        info += "Contains " + \
            str(self.rtc.d_control.size()) + " Controller(s):\n"
        info += "Control # | Type  | Nslope | Nactu\n"

        for i in range(< int > self.rtc.d_control.size()):
            control = self.rtc.d_control[i]
            info += "%9d" % (i + 1) + " | " + "%5s" % control.get_type() + " | " + "%6d" % control.nslope() + \
                    " | " + str(control.nactu()) + "\n"

        info += "--------------------------------------------------------"
        return info


IF USE_BRAMA == 1:
    #################################################
    # P-Class Rtc_brama
    #################################################
    # child constructor must have the same prototype (same number of
    # non-optional arguments)
    cdef class Rtc_brama(Rtc):

        def __cinit__(self, Sensors sensor, Target target=None, device=-1):
            del self.rtc

            cdef carma_context * context = &carma_context.instance()
            if target is not None:
                self.rtc = new sutra_rtc_brama(context, sensor.sensors, target.target, "rtc_brama")
            else:
                self.rtc = new sutra_rtc_brama(context, sensor.sensors, NULL, "rtc_brama")

        def __dealloc__(self):
            pass  # del self.rtc # TODO

        cpdef publish(self):
            cdef sutra_rtc_brama * rtc = < sutra_rtc_brama * > (self.rtc)
            rtc.publish()
