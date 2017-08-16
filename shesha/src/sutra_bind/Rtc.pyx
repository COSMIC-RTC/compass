include "../../par.pxi"

import numpy as np
cimport numpy as np
# np.import_array()
import os
import shesha_constants as scons

from cython.operator cimport dereference as deref, preincrement as inc


cdef class Rtc:
    # child constructor must have the same prototype (same number of
    # non-optional arguments)
    def __cinit__(self, naga_context context, Sensors sensor=None, Target target=None, device=-1):
        self.context = context
        if(device == -1):
            device = context.get_activeDevice()
        context.set_activeDevice(device, 1)
        self.device = device

        self.rtc = new sutra_rtc(self.context.c)

    def __dealloc__(self):
        del self.rtc

    def __str__(self):
        print("RTC object:")

        cdef sutra_centroider * centro
        cdef sutra_controller * contro
        cdef int i

        info = "Contains " + str(self.rtc.d_centro.size()) + " Centroider(s)\n"
        info += "Centro # | Type  | nwfs | Nvalid\n"

        for i in range( < int > self.rtc.d_centro.size()):
            centro = self.rtc.d_centro[i]
            info += "%8d" % (i + 1) + " | " + "%5s" % centro.get_type() + " | " + "%4d" % (centro.nwfs + 1) + \
                    " | " + str(centro.nvalid) + "\n"

        info += "Contains " + \
            str(self.rtc.d_control.size()) + " Controller(s):\n"
        info += "Control # | Type  | Nslope | Nactu\n"

        for i in range( < int > self.rtc.d_control.size()):
            control = self.rtc.d_control[i]
            info += "%9d" % (i + 1) + " | " + "%5s" % control.get_type() + " | " + "%6d" % control.nslope() + \
                    " | " + str(control.nactu()) + "\n"

        info += "--------------------------------------------------------"
        return info

    def add_centroider(self, Sensors sensor, long nwfs, long nvalid,
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
        cdef int activeDevice = self.rtc.device
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        self.rtc.add_centroider(
            sensor.sensors,
            nwfs,
            nvalid,
            offset,
            scale,
            activeDevice,
            type_centro)

    def add_controller(self, int nactu, float delay, bytes type_control, Dms dms,
                       list dmseen, np.ndarray[ndim=1, dtype=np.float32_t] alt,
                       int ndm, long Nphi=-1, bool wfs_direction=False):
        """Add a controller in the sutra_controller vector of the RTC on the GPU

        :parameters:
            nactu: (int) : number of actuators

            delay: (float) : loop delay

            type_control: (str) : controller's type

            dms: (Dms) : sutra_dms object (GPU)

            type_dmseen: (list) : dms indices controled by the controller

            alt: (np.ndarray[ndim=1,dtype=np.float32_t]) : altitudes of the dms seen

            ndm: (int) : number of dms controled

            Nphi: (long) : number of pixels in the pupil (used in geo controler case only)
        """
        type_dmseen = < char ** > malloc(len(dmseen) * sizeof(char *))
        for j in range(len(dmseen)):
            type_dmseen[j] = dmseen[j]

        cdef float * ptr_alt = < float * > alt.data
        cdef char * ptr_dmseen = < char * > type_dmseen
        self.context.set_activeDeviceForCpy(self.device, 1)
        if(Nphi > -1):
            self.rtc.add_controller_geo(nactu, Nphi, delay, self.device, dms.dms, & ptr_dmseen, ptr_alt, ndm, wfs_direction)
        else:
            self.rtc.add_controller(nactu, delay, self.device, type_control, dms.dms, & ptr_dmseen, ptr_alt, ndm)

    def rm_controller(self):
        """Remove a controller"""
        self.context.set_activeDevice(self.device, 1)
        self.dms.rm_controller()

    def init_npix(self, int ncentro):
        """Initialize npix in the sutra_centroider_corr object (useless ?)

        :parameters:
            ncentro: (int) : centroider's index
        """

        self.context.set_activeDevice(self.device, 1)
        cdef sutra_centroider * centro = NULL
        centro = self.rtc.d_centro.at(ncentro)
        cdef sutra_centroider_corr * corr
        if(self.rtc.d_centro[ncentro].is_type(scons.CentroiderType.CORR)):
            corr = dynamic_cast_centroider_corr_ptr(centro)
            corr.init_bincube()

    def init_weights(self, int ncentro, np.ndarray[ndim=3, dtype=np.float32_t] w):
        """Load the weight array in sutra_centroider_wcog object

        :parameters:
            ncentro: (int) : centroider's index
            w: (np.ndarray[ndim=3, dtype=np.float32_t]) : weight
        """
        self.context.set_activeDevice(self.rtc.device, 1)
        cdef sutra_centroider_wcog * centro = \
            dynamic_cast_centroider_wcog_ptr(self.rtc.d_centro[ncentro])
        centro.init_weights()
        cdef np.ndarray w_F = w.T.copy()
        centro.load_weights(< float * > w_F.data, int(w.ndim))

    def init_corr(self, int ncentro, np.ndarray[ndim=2, dtype=np.float32_t] w,
                  np.ndarray[ndim=2, dtype=np.float32_t] corr_norm,
                  int sizex, int sizey,
                  np.ndarray[ndim=2, dtype=np.float32_t] interpmat):
        """Initialize sutra_centroider_corr oblect

        :parameters:
            ncentro: (int) : centroider's index
            w: (np.ndarray[ndim=1,dtype=np.float32_t]) : weight
            corr_norm: (np.ndarray[ndim=2,dtype=np.float32_t]) :
            sizex: (int) :
            sizey: (int) :
            interpmat: ([ndim=2,dtype=np.float32_t]) :
        """
        # TODO: verify the w dimensions needed...
        self.context.set_activeDevice(self.rtc.device, 1)

        cdef sutra_centroider_corr * centro_corr = dynamic_cast_centroider_corr_ptr(
            self.rtc.d_centro[ncentro])
        cdef np.ndarray w_F = w.T.copy()
        cdef np.ndarray[dtype= np.float32_t] corr_norm_F = corr_norm.T.copy()
        cdef np.ndarray[dtype= np.float32_t] interpmat_F = interpmat.T.copy()
        centro_corr.init_corr(sizex, sizey, < float * > interpmat_F.data)
        centro_corr.load_corr( < float * > w_F.data, < float * > corr_norm_F.data, int(w.ndim))

    def do_centroids(self, int ncontrol=-1):
        """Compute the centroids with sutra_controller #ncontrol object

        :parameters:
            ncontrol: (optional) controller's index
        """
        if(ncontrol > -1):
            self.rtc.do_centroids(ncontrol)
        else:
            self.rtc.do_centroids()

    def do_centroids_geom(self, int ncontrol=-1):
        """Compute the geometric centroids with sutra_controller #ncontrol object

        :parameters:
            ncontrol: (optional) controller's index
        """
        if(ncontrol > -1):
            self.rtc.do_centroids_geom(ncontrol)
        else:
            self.rtc.do_centroids(0)

    def init_proj(self, int ncontrol, Dms dms, np.ndarray[ndim=1, dtype=np.int32_t] indx_dm,
                  np.ndarray[ndim=1, dtype=np.float32_t] unitpervolt,
                  np.ndarray[ndim=1, dtype=np.int32_t] indx_pup, np.ndarray[ndim=1, dtype=np.int32_t] indx_mpup, int roket=0):
        """Initialize the projection matrix for sutra_controller_geo object.
        The projection matrix is (IFt.IF)**(-1) * IFt where IF is the DMs influence functions matrix

        :parameters:
            ncontrol: (int) : controller index

            dms: (Dms) : Dms object

            indx_dm: (np.ndarray[ndim=1,dtype=np.int32_t]) : indices of where(pup) on DM screen

            unitpervolt: (np.ndarray[ndim=1,dtype=np.float32_t]) : unitpervolt DM parameter

            indx_pup: (np.ndarray[ndim=1,dtype=np.int32_t]) : indices of where(pup) on ipupil screen

            roket : (int) : optimisation flag for ROKET
        """
        cdef sutra_controller_geo * controller_geo = \
            dynamic_cast_controller_geo_ptr(self.rtc.d_control[ncontrol])

        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        controller_geo.init_proj_sparse(dms.dms, < int * > indx_dm.data,
                                        < float * > unitpervolt.data, < int * > indx_pup.data, < int * > indx_mpup.data, roket)

    def init_modalOpti(self, int ncontrol, int nmodes, int nrec, np.ndarray[ndim=2, dtype=np.float32_t] M2V,
                       float gmin, float gmax, int ngain, float Fs):
        """Initialize the modal optimization controller : compute the slopes-to-modes matrix
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
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls

        cdef np.ndarray[dtype= np.float32_t] M2V_F = M2V.T.copy()
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.init_modalOpti(nmodes, nrec, < float * > M2V_F.data, gmin, gmax, ngain, Fs)
        else:
            raise TypeError(
                "**** ERROR : Modal Optimization only for controller type ls ****")

    def load_open_loop_slopes(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] ol_slopes):
        """Load an array of recoded open-loop measurements for modal optimization

        :parameters:
            ncontrol: (int) : controller index

            ol_slopes: (np.ndarray[ndim=2, dtype=np.float32_t]) : open-loop slopes
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls

        cdef np.ndarray[dtype= np.float32_t] slopes_F = ol_slopes.T.copy()
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.loadOpenLoopSlp(< float * > slopes_F.data)
        else:
            raise TypeError("Controller type must be ls")

    def modal_control_optimization(self, int ncontrol):
        """Compute the command matrix with modal control optimization

        :parameter:
            ncontrol: controller index
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        if(type_control == scons.ControllerType.LS):
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

    def init_cured(self, int ncontrol, int nxsub, np.ndarray[ndim=2, dtype=np.int32_t] isvalid,
                   int ndivs, bool tt_flag):
        """ Initialize the CURED controller

        :parameters:
            ncontrol: (int) : controller index
            nxsub: (int): number of subap. in the WFS diameter
            isvalid: (np.ndarray[ndim=2, dtype=np.int32]): array where 1 define a valid subap
            ndivs: (int): cured ndivs
            tt_flag: (bool): True if a TT DM exists
        """
        self.context.set_activeDevice(self.device, 1)
        cdef sutra_controller_cured * controller_cured

        controller_cured = dynamic_cast_controller_cured_ptr(
            self.rtc.d_control[0])
        controller_cured.init_cured(nxsub, < int * > isvalid.data, ndivs, tt_flag)

    def filter_cmat(self, int ncontrol, float cond):
        """ Filter the cmat from TT for MV controller

        :parameters:
            ncontrol : (int) : controller index
            cond: (float): TT conditioning
        """
        self.context.set_activeDevice(self.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.filter_cmat(cond)
        else:
            raise TypeError("Controller needs to be mv")

    def filter_cphim(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] F,
                     np.ndarray[ndim=2, dtype=np.float32_t] Nact):
        """ Filter the Cphim from piston and actuators coupling for MV controller

        :parameters:
            ncontrol : (int) : controller index
            F : (np.ndarray[ndim=2, dtype=np.float32]): piston filter matrix
            Nact : (np.ndarray[ndim=2, dtype=np.float32]): Coupling matrix
        """
        self.context.set_activeDevice(self.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.filter_cphim( < float * > F.data, < float * > Nact.data)
        else:
            raise TypeError("Controller needs to be mv")

    def compute_Cmm(self, int ncontrol, Atmos atmos, Sensors wfs,
                    np.ndarray[ndim=1, dtype=np.float64_t] L0,
                    np.ndarray[ndim=1, dtype=np.float64_t] frac,
                    np.ndarray[ndim=1, dtype=np.float64_t] alphaX,
                    np.ndarray[ndim=1, dtype=np.float64_t] alphaY,
                    float diam, float cobs):
        """
            Compute the Cmm matrix for the MV controller

        :parameters:
            ncontrol: (int): controller index
            atmos: (Atmos): Atmos object
            wfs: (Sensors): Wfs object
            L0: (np.ndarray[ndim=1, dtype=np.float64]): L0 of each layers
            frac : (np.ndarray[ndim=1, dtype=np.float64]): p_atmos.frac * (p_atmos.r0 ** (-5 / 3)))
            alphaX: (np.ndarray[ndim=1, dtype=np.float64]): X position of WFS GS in rad
            alphaY: (np.ndarray[ndim=1, dtype=np.float64]): Y position of WFS GS in rad
            diam :(float): telescope diameter
            cobs: (flaot): central obstruction
        """
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.compute_Cmm(atmos.s_a, wfs.sensors, < double * > L0.data,
                                      < double * > frac.data, < double * > alphaX.data,
                                      < double * > alphaY.data,
                                      < double > diam, < double > cobs)
        else:
            raise TypeError("Controller needs to be mv")

    def compute_Cphim(self, int ncontrol, Atmos atmos, Sensors wfs, Dms dms,
                      np.ndarray[ndim=1, dtype=np.float64_t] L0,
                      np.ndarray[ndim=1, dtype=np.float64_t] frac,
                      np.ndarray[ndim=1, dtype=np.float64_t] alphaX,
                      np.ndarray[ndim=1, dtype=np.float64_t] alphaY,
                      np.ndarray[ndim=1, dtype=np.float64_t] X,
                      np.ndarray[ndim=1, dtype=np.float64_t] Y,
                      np.ndarray[ndim=1, dtype=np.float64_t] Xactu,
                      np.ndarray[ndim=1, dtype=np.float64_t] Yactu,
                      float diam,
                      np.ndarray[ndim=1, dtype=np.float64_t] k2,
                      np.ndarray[ndim=1, dtype=np.int32_t] NlayersDM,
                      np.ndarray[ndim=1, dtype=np.int32_t] indlayersDM,
                      float FoV,
                      np.ndarray[ndim=1, dtype=np.float64_t] pitch,
                      np.ndarray[ndim=1, dtype=np.float64_t] alt_DM):
        """
            Compute the Cphim matrix for the MV controller

        :parameters:
            ncontrol: (int): controller index
            atmos: (Atmos): Atmos object
            wfs: (Sensors): Wfs object
            dms: (Dms): Dms object
            L0: (np.ndarray[ndim=1, dtype=np.float64]): L0 of each layers
            frac : (np.ndarray[ndim=1, dtype=np.float64]): p_atmos.frac * (p_atmos.r0 ** (-5 / 3)))
            alphaX: (np.ndarray[ndim=1, dtype=np.float64]): X position of WFS GS in rad
            alphaY: (np.ndarray[ndim=1, dtype=np.float64]): Y position of WFS GS in rad
            X: (np.ndarray[ndim=1, dtype=np.float64]): X position of WFS subap [m]
            Y: (np.ndarray[ndim=1, dtype=np.float64]): Y position of WFS subap [m]
            Xactu: (np.ndarray[ndim=1, dtype=np.float64]): X position of DMs actus [m]
            Yactu: (np.ndarray[ndim=1, dtype=np.float64]): Y position of DMs actus [m]
            diam :(float): telescope diameter
            k2: (np.ndarray[ndim=1, dtype=np.float64]): Normalizations constants
            NlayersDM: (np.ndarray[ndim=1, dtype=np.int32]): number of turbu. layers handled by each DMs
            indlayersDM: (np.ndarray[ndim=1, dtype=np.int32]): index of turbu. layers handled by each DMs
            FoV : (float): Field of view [rad]
            pitch: (np.ndarray[ndim=1, dtype=np.float64]): pitch of each DMs
            alt_DM: (np.ndarray[ndim=1, dtype=np.float64]): altitudes of each DMs
        """
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.compute_Cphim(atmos.s_a, wfs.sensors, dms.dms, < double * > L0.data,
                                        < double * > frac.data, < double * > alphaX.data,
                                        < double * > alphaY.data,
                                        < double * > X.data,
                                        < double * > Y.data,
                                        < double * > Xactu.data,
                                        < double * > Yactu.data,
                                        < double > diam,
                                        < double * > k2.data,
                                        < long * > NlayersDM.data,
                                        < long * > indlayersDM.data,
                                        < double > FoV,
                                        < double * > pitch.data,
                                        < double * > alt_DM.data)
        else:
            raise TypeError("Controller needs to be mv")

    def do_imat_geom(self, int ncontrol, Dms dms, int geom):
        """Compute the interaction matrix by using a geometric centroiding method

        :parameters:
            ncontrol: (int) : controller index

            dms: (Dms) : Dms object

            geom: (int) : type of geometric method (0 or 1)
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.rtc.device, 1)
        self.rtc.do_imat_geom(ncontrol, dms.dms, geom)

    def do_centroids_ref(self, int ncontrol):
        """
            TODO: docstring
        :parameters:
        ncontrol: (int) : controller index
        """

        cdef carma_obj[float] * phase
        cdef sutra_wfs * wfs
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] h_ref
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] h_rawslp
        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        cdef int nslope = control.nslope()
        cdef float tmp

        print("Doing refslp...")
        for idx_cntr in range( < int > self.rtc.d_centro.size()):
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

    def do_imat_kl(self, int ncontrol, Dms dms, p_dms, np.ndarray[ndim=2, dtype=np.float32_t] kl, bool ntt):
        """Compute the interaction matrix in the KL basis

        :parameters:
            ncontrol: (int) : controller index
            dms: (Dms) : Dms object
            ntt: (bool): True if a tt dm exists
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        cdef carma_obj[float] * d_imat = NULL
        cdef carma_obj[float] * d_imat_ret = NULL
        cdef long dims[3]
        cdef int nactu = control.nactu()
        cdef int nslope = control.nslope()
        cdef int inds1, j, idx_cntr, device
        cdef float tmp_noise
        cdef float tmp_nphot
        inds1 = 0
        cdef sutra_dm * dm
        cdef sutra_wfs * wfs
        cdef carma_obj[float] * screen
        cdef vector[sutra_dm *].iterator it_dm
        cdef float * d_centroids
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] h_centroids
        cdef carma_obj[float] * phase
        cdef bytes type_control = control.get_type()
        cdef bytes type_dm
        if(type_control == scons.ControllerType.LS):
            d_imat = dynamic_cast_controller_ls_ptr(control).d_imat
        elif(type_control == scons.ControllerType.MV):
            d_imat = dynamic_cast_controller_mv_ptr(control).d_imat
        else:
            # Create a temporary imat to return
            dims[0] = 2
            dims[1] = nactu
            dims[2] = nslope
            d_imat_ret = new carma_obj[float](self.context.c, dims)
            d_imat = d_imat_ret
            print(
                "WARNING: the controller is NOT a LS or MV, the imat computed will be returned")

        it_dm = control.d_dmseen.begin()
        dm = deref(it_dm)
        while(it_dm != control.d_dmseen.end()):
            dm = deref(it_dm)
            type_dm = dm.type
            if(type_dm == scons.DmType.TT):
                dmtt = deref(it_dm)
        cdef int cc = 0
        print("Doing imat_kl...%d%%\r" % cc, end='')

        for j in range(kl.shape[1]):
            # push actu __________________________________________________

            dms.set_full_comm(np.float32(kl[:, j].copy()))

            for idx_cntr in range( < int > self.rtc.d_centro.size()):
                wfs = self.rtc.d_centro[idx_cntr].wfs
                screen = wfs.d_gs.d_phase.d_screen
                tmp_noise = wfs.noise
                tmp_nphot = wfs.nphot
                wfs.nphot = wfs.nphot4imat
                wfs.noise = -1
                wfs.kernconv = True
                wfs.sensor_trace(dms.dms, 1)
                wfs.comp_image()
                wfs.noise = tmp_noise
                wfs.nphot = tmp_nphot
                wfs.kernconv = False

            self.rtc.do_centroids(ncontrol, True)

            h_centroids = self.getCentroids(ncontrol)
            control.d_centroids.host2device(< float * > h_centroids.data)

            device = control.d_centroids.getDevice()
            d_centroids = control.d_centroids.getData()
            if (ntt and (j >= kl.shape[1] - 2)):
                convert_centro(d_centroids,
                               d_centroids, 0,
                               (1. / dmtt.push4imat),
                               control.d_centroids.getNbElem(),
                               self.context.c.get_device(device))
            else:
                convert_centro(d_centroids,
                               d_centroids, 0,
                               (1. / dm.push4imat),
                               control.d_centroids.getNbElem(),
                               self.context.c.get_device(device))

            control.d_centroids.copyInto(
                d_imat.getDataAt(inds1),
                control.nslope())
            for i in range(np.size(p_dms)):
                dms.resetdm(p_dms[i].type_dm, p_dms[i].alt)

            inds1 += control.nslope()
            cc = cc + 1
            print("Doing imat... #%d/%d \r" % (cc, kl.shape[1]), end=' ')

        print("Done")

        cdef np.ndarray[ndim= 2, dtype = np.float32_t] h_imat_ret
        if d_imat_ret != NULL:
            h_imat_ret = np.zeros((nactu, nslope), dtype=np.float32)
            d_imat_ret.device2host(< float * > h_imat_ret.data)
            del d_imat_ret
            return h_imat_ret

    def do_imat(self, int ncontrol, Dms dms):
        """Compute the interaction matrix

        :parameters:
            ncontrol: (int) : controller index
            dms: (Dms) : Dms object
        """

        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        cdef carma_obj[float] * d_imat = NULL
        cdef carma_obj[float] * d_imat_ret = NULL
        cdef long dims[3]
        cdef int nactu = control.nactu()
        cdef int nslope = control.nslope()
        cdef bytes type_control = control.get_type()

        if(type_control == scons.ControllerType.LS):
            d_imat = dynamic_cast_controller_ls_ptr(control).d_imat
        elif(type_control == scons.ControllerType.MV):
            d_imat = dynamic_cast_controller_mv_ptr(control).d_imat
        else:
            # Create a temporary imat to return
            dims[0] = 2
            dims[1] = nactu
            dims[2] = nslope
            d_imat_ret = new carma_obj[float](self.context.c, dims)
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
        cdef vector[sutra_dm *].iterator it_dm
        cdef float * d_centroids
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] h_centroids
        cdef carma_obj[float] * phase

        cc = 0
        it_dm = control.d_dmseen.begin()
        print("Doing imat...%d%% \r" % cc, end='')
        while(it_dm != control.d_dmseen.end()):
            dm = deref(it_dm)
            for j in range(dm.ninflu):
                # push actu __________________________________________________
                dm.comp_oneactu(j, dm.push4imat)
                for idx_cntr in range( < int > self.rtc.d_centro.size()):
                    wfs = self.rtc.d_centro[idx_cntr].wfs
                    screen = wfs.d_gs.d_phase.d_screen
                    tmp_noise = wfs.noise
                    tmp_nphot = wfs.nphot
                    wfs.nphot = wfs.nphot4imat
                    wfs.noise = -1
                    wfs.kernconv = True
                    wfs.sensor_trace(dms.dms, 1)
                    wfs.comp_image()
                    wfs.noise = tmp_noise
                    wfs.nphot = tmp_nphot
                    wfs.kernconv = False

                self.rtc.do_centroids(ncontrol, True)

                h_centroids = self.getCentroids(ncontrol)
                control.d_centroids.host2device(< float * > h_centroids.data)

                device = control.d_centroids.getDevice()
                d_centroids = control.d_centroids.getData()

                convert_centro(d_centroids,
                               d_centroids, 0,
                               (0.5 / dm.push4imat),
                               control.d_centroids.getNbElem(),
                               self.context.c.get_device(device))

                control.d_centroids.copyInto(
                    d_imat.getDataAt(inds1),
                    control.nslope())
                dm.reset_shape()

                # pul actu __________________________________________________
                dm.comp_oneactu(j, -1. * dm.push4imat)

                for idx_cntr in range( < int > self.rtc.d_centro.size()):
                    wfs = self.rtc.d_centro[idx_cntr].wfs
                    tmp_noise = wfs.noise
                    tmp_nphot = wfs.nphot
                    wfs.nphot = wfs.nphot4imat
                    wfs.noise = -1
                    wfs.kernconv = True
                    wfs.sensor_trace(dms.dms, 1)
                    wfs.comp_image()
                    wfs.noise = tmp_noise
                    wfs.nphot = tmp_nphot
                    wfs.kernconv = False

                self.rtc.do_centroids(ncontrol, True)

                h_centroids = self.getCentroids(ncontrol)
                control.d_centroids.host2device(< float * > h_centroids.data)

                device = control.d_centroids.getDevice()
                d_centroids = control.d_centroids.getData()

                convert_centro(d_centroids,
                               d_centroids, 0,
                               (0.5 / dm.push4imat),
                               control.d_centroids.getNbElem(),
                               self.context.c.get_device(device))

                carma_axpy[float](self.context.c.get_cublasHandle(),
                                  control.d_centroids.getNbElem(), < float > -1.0,
                                  d_centroids, 1,
                                  d_imat.getDataAt(inds1), 1)

                dm.reset_shape()

                inds1 += control.nslope()
                cc = cc + 1
                print("Doing imat... #%d/%d \r" % (cc, nactu), end=' ')

            inc(it_dm)
        print("imat done\n")

        cdef np.ndarray[ndim= 2, dtype = np.float32_t] h_imat_ret
        if d_imat_ret != NULL:
            h_imat_ret = np.zeros((nactu, nslope), dtype=np.float32)
            d_imat_ret.device2host(< float * > h_imat_ret.data)
            del d_imat_ret
            return h_imat_ret

    def comp_slopes(self, int ncentro):
        """Compute the slopes in a sutra_wfs object. This function is equivalent to
        docentroids() but the centroids are stored in the sutra_wfs object instead of
        the sutra_rtc object

        :parameters:
            ncentro: (int) : centroider index
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)

        self.rtc.d_centro[ncentro].get_cog()

    def imat_svd(self, int ncontrol):
        """Compute the singular value decomposition of the interaction matrix

        :param ncontrol: controller index
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef sutra_controller_ls * controller_ls
        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            if(controller_ls.svdec_imat() == 1):
                raise RuntimeError("sutra controller has no SVD implementation")
        else:
            raise TypeError("Controller needs to be ls")

    def build_cmat(self, int ncontrol, int nfilt, int filt_tt=0):
        """Compute the command matrix in a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index

            nfilt: (int) : number of modes to filter

            filt_tt: (int) : (optional) flag to filter TT
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            if(filt_tt > 0):
                controller_ls.build_cmat(nfilt, True)
            else:
                controller_ls.build_cmat(nfilt)

    def build_cmat_mv(self, int ncontrol, float cond):
        """Compute the command matrix in a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index

            cond: (float) : conditioning factor for the Cmm inversion
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.build_cmat(cond)

    def load_Cn(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] N):
        """Load the noise vector on a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index

            N: (np.ndarray[ndim=1,dtype=np.float32_t]) : noise vector
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.load_noisemat(< float * > N.data)

    def do_clipping(self, int ncontrol, float min, float max):
        """Clip the command to apply on the DMs on a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            min: (float) : minimum value for the command
            max: (float) : maximum value for the command
        """
        self.context.set_activeDevice(self.rtc.device, 1)
        self.rtc.do_clipping(ncontrol, min, max)

    def do_control(self, int ncontrol):
        """Compute the command to apply on the DMs on a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        """
        self.context.set_activeDevice(self.rtc.device, 1)
        self.rtc.do_control(ncontrol)

    def do_control_geo(self, int ncontrol, Dms dms, Target target, int ntarget):
        """Compute the command to apply on the DMs on a sutra_controller_geo object
        for the target direction

        :parameters:
            ncontrol: (int) : controller index
            dms: (shesha_dm object) : shesha_dm
            target : (shesha_target) : shesha_target
            ntarget : (int) : target number
        """
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        self.context.set_activeDevice(self.rtc.device, 1)

        if(type_control == scons.ControllerType.GEO):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            controller_geo.comp_dphi(target.target.d_targets[ntarget], False)
            self.rtc.do_control(ncontrol)
        else:
            raise TypeError("Controller needs to be geo")

    def do_control_geo_on_wfs(self, int ncontrol, Dms dms, Sensors sensors, int nwfs):
        """Compute the command to apply on the DMs on a sutra_controller_geo object
        for the wfs direction

        :parameters:
            ncontrol: (int) : controller index
            dms: (shesha_dm object) : shesha_dm
            sensors : (shesha_sensor) : shesha_sensor
            nwfs : (int) : wfs number

        """
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        self.context.set_activeDevice(self.rtc.device, 1)

        if(type_control == scons.ControllerType.GEO):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            controller_geo.comp_dphi(sensors.sensors.d_wfs[nwfs].d_gs, True)
            self.rtc.do_control(ncontrol)
        else:
            raise TypeError("Controller needs to be geo")

    def load_Btt(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] Btt):
        """Load the Btt basis for sutra_controller_geo projection in case of error_budget

        :parameters:
            ncontrol: (int) : controller index
            Btt: (np.ndarray[ndim=2,dtype=np.float32_t]) : Btt basis
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_contro = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray TT = Btt[-2:, -2:].T.copy()
        cdef np.ndarray Btt_F = Btt[:-2, :-2].T.copy()
        cdef np.ndarray TT_F = TT.T.copy()

        if(type_contro == scons.ControllerType.GEO):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            controller_geo.load_Btt( < float*>Btt_F.data, < float*>TT_F.data)
        else:
            raise TypeError("Controller needs to be geo")

    def apply_control(self, int ncontrol, Dms dms):
        """Compute the DMs shapes from the commands computed in a sutra_controller_object.
        From the command vector, it computes the voltage command (adding perturbation voltages,
        taking delay into account) and then apply it to the dms

        :parameters:
            ncontrol: (int) : controller index
        """
        self.context.set_activeDevice(self.rtc.device, 1)
        self.rtc.apply_control(ncontrol, dms.dms)

    """
           _______. _______ .___________.                _______  _______ .___________.
          /       ||   ____||           |     ___       /  _____||   ____||           |
         |   (----`|  |__   `---|  |----`    ( _ )     |  |  __  |  |__   `---|  |----`
          \   \    |   __|      |  |         / _ \/\   |  | |_ | |   __|      |  |
      .----)   |   |  |____     |  |        | (_>  <   |  |__| | |  |____     |  |
      |_______/    |_______|    |__|         \___/\/    \______| |_______|    |__|

    """

    def get_pyr_method(self, int n):
        """Get the pyramid centroiding method
        :parameters:
        n : (int) : pyr centroider number
        """
        cdef sutra_centroider_pyr * centro = NULL

        if(self.rtc.d_centro[n].is_type(scons.WFSType.PYRHR)):
            centro = dynamic_cast_centroider_pyr_ptr(self.rtc.d_centro[n])
            return centro.get_method_str()
        else:
            e = "Centroider should be pyrhr, got " + \
                self.rtc.d_centro[n].get_type()
            raise ValueError(e)

    def set_pyr_method(self, int n, int method, list p_centroiders):
        """Set the pyramid centroiding method
        :parameters:
        n : (int) : pyr centroider number
        method : (int) : new centroiding method (0: nosinus global
                                                 1: sinus global
                                                 2: nosinus local
                                                 3: sinus local)
        p_centroiders : (list of Param_centroider) : list of centroider parameters
        """
        cdef sutra_centroider_pyr * centro = NULL

        if(self.rtc.d_centro[n].is_type(scons.WFSType.PYRHR)):
            if method >= Other:
                raise ValueError("method unknown")

            centro = dynamic_cast_centroider_pyr_ptr(self.rtc.d_centro[n])
            centro.set_method(method)
            p_centroiders[n].set_method(method)
        else:
            e = "Centroider should be pyrhr, got " + \
                self.rtc.d_centro[n].get_type()
            raise ValueError(e)

    def set_pyr_thresh(self, int n, float threshold, list p_centroiders):
        """Set the pyramid threshold
        :parameters:
        n : (int) : pyr centroider number
        threshold : (float) : new threshold in photons
        p_centroiders : (list of Param_centroider) : list of centroider parameters
        """
        cdef sutra_centroider_pyr * centro = NULL
        self.context.set_activeDevice(self.device, 1)

        if(self.rtc.d_centro[n].is_type(scons.WFSType.PYRHR)):
            centro = dynamic_cast_centroider_pyr_ptr(self.rtc.d_centro[n])
            centro.set_valid_thresh(threshold)

            p_centroiders[n].set_thresh(threshold)
        else:
            e = "Centroider should be pyrhr, got " + \
                self.rtc.d_centro[n].get_type()
            raise ValueError(e)

    def set_pyr_ampl(self, int n, np.ndarray[ndim=1, dtype=np.float32_t] cx,
                     np.ndarray[ndim=1, dtype=np.float32_t] cy, float scale):
        """Set the pyramid modulation amplitude
        :parameters:
        n : (int) : pyr centroider number
        cx : (np.ndarray[ndim=1, dtype=np.float32]) : X position of new modulation points
        cy : (np.ndarray[ndim=1, dtype=np.float32]) : Y position of new modulation points
        scale : (float) : scaling factor
        """
        self.context.set_activeDevice(self.device, 1)

        if(self.rtc.d_centro[n].is_type(scons.WFSType.PYRHR)):
            centro = self.rtc.d_centro.at(n)
            pyr = dynamic_cast_wfs_pyr_ptr(centro.wfs)
            pyr.pyr_cx.fill_from( < float * >cx.data)
            pyr.pyr_cy.fill_from( < float * >cy.data)
        else:
            e = "Centroider should be pyrhr, got " + \
                self.rtc.d_centro[n].get_type()
            raise ValueError(e)

    def set_thresh(self, int ncentro, float thresh):
        """set threshold for the centroider #ncentro

        :parameters:
            ncentro: (int) : centroider's index
            thresh: (float) : threshold
        """

        cdef sutra_centroider_tcog * tcog = NULL
        cdef sutra_centroider * centro = NULL
        self.context.set_activeDeviceForCpy(self.device, 1)
        if(self.rtc.d_centro[ncentro].is_type(scons.CentroiderType.TCOG)):
            centro = self.rtc.d_centro.at(ncentro)
            tcog = dynamic_cast_centroider_tcog_ptr(centro)
            tcog.set_threshold(thresh)

    def get_centroids(self, int ncontrol, Sensors g_wfs=None, int nwfs=0):
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

        cdef const long * dims
        cdef sutra_wfs * wfs
        cdef sutra_rtc * rtc = self.rtc
        cdef carma_obj[float] * d_tmp
        cdef carma_obj[float] * d_data
        cdef carma_context * context = &carma_context.instance()
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data

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
            d_data.device2host(< float * > data.data)

        else:
            dims = rtc.d_control[ncontrol].d_centroids.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            rtc.d_control[ncontrol].d_centroids.device2host(< float * > data.data)
        return data

    def set_gain(self, int ncontrol, float gain):
        """Set the loop gain in sutra_controller object

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

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.set_gain(gain)

        elif(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.set_gain(gain)

        elif(type_control == scons.ControllerType.CURED):
            controller_cured = dynamic_cast_controller_cured_ptr(
                self.rtc.d_control[ncontrol])
            controller_cured.set_gain(gain)

        elif(type_control == scons.ControllerType.GEO):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            controller_geo.set_gain(gain)
        else:
            raise TypeError(
                "Controller needs to be ls, mv, cured, geo, kalman_GPU or kalman_CPU. for generic (g=1.0) use set_mgain")

    def set_mgain(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] mgain):
        """Set modal gains in sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            mgain: (np.ndarray[ndim=1,dtype=np.float32_t]) : modal gains
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.set_mgain(< float * > mgain.data)
        elif(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.set_mgain(< float * > mgain.data)
        elif(type_control == scons.ControllerType.GENERIC):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.set_mgain(< float * > mgain.data)
        else:
            raise TypeError("Controller needs to be ls, mv or generic")

    def set_scalar_mgain(self, int ncontrol, float g):
        """
            Set modal gains all equal to the scalar given
        :parameters:
            ncontrol: (int) : controller index
            g: (float): gain
        """
        x = self.get_mgain(ncontrol)
        x[:] = g
        self.set_mgain(ncontrol, x.copy())

    def set_com(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] comvec):
        """Set the command vector of a sutra_controller object to comvec

        :parameters:
            ncontrol: (int) : controller index
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        self.rtc.d_control[ncontrol].d_com.host2device(< float * > comvec.data)

    def set_centroids(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] centro):
        """Set the centroids vector of a sutra_controller object to centro

        :parameters:
            ncontrol: (int) : controller index
            centro: (np.ndarray[ndim=1,dtype=np.float32_t]) : centroids vector
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        self.rtc.d_control[ncontrol].d_centroids.host2device(< float * > centro.data)

    def get_mgain(self, int ncontrol):
        """Return modal gains from sutra_controller

        :parameters:
            ncontrol: (int) : controller index
        :return:
            mgain : (np.ndarray[ndim=1,dtype=np.float32_t]) : modal gains
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef int size
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] mgain

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            size = controller_ls.d_gain.getNbElem()
            mgain = np.zeros(size, dtype=np.float32)
            controller_ls.d_gain.device2host(< float * > mgain.data)
            return mgain
        elif(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            size = controller_mv.d_gain.getNbElem()
            mgain = np.zeros(size, dtype=np.float32)
            controller_mv.d_gain.device2host(< float * > mgain.data)
            return mgain
        elif(type_control == scons.ControllerType.GENERIC):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            size = controller_generic.d_gain.getNbElem()
            mgain = np.zeros(size, dtype=np.float32)
            controller_generic.d_gain.device2host(< float * > mgain.data)
            return mgain
        else:
            raise TypeError("Controller needs to be ls, generic or mv")

    def set_imat(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """Set the interaction matrix on a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            data: (np.ndarray[ndim=2,dtype=np.float32_t]) : interaction matrix to use
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef np.ndarray[dtype= np.float32_t] data_F = data.T.copy()

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.d_imat.host2device(< float * > data_F.data)
        elif(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.d_imat.host2device(< float * > data_F.data)
        else:
            raise TypeError("Controller needs to be ls or mv")

    def get_imat(self, int ncontrol):
        """Return the interaction matrix of a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            imat : (np.ndarray[ndim=2,dtype=np.float32_t]) : interaction matrix

        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_cured * controller_cured
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        cdef const long * dims = NULL
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] imat_F

        if(type_control == scons.ControllerType.GENERIC or type_control == scons.ControllerType.GEO):
            raise TypeError("Generic controller doesn't have imat")

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.d_imat.getDims()
            imat_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_ls.d_imat.device2host(< float * > imat_F.data)

        elif(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_mv.d_imat.getDims()
            imat_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_mv.d_imat.device2host(< float * > imat_F.data)

        elif(type_control == scons.ControllerType.CURED):
            controller_cured = dynamic_cast_controller_cured_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_cured.d_imat.getDims()
            imat_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_cured.d_imat.device2host(< float * > imat_F.data)

        return imat_F.T.copy()

    def set_cmat(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """Set the command matrix on a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
            data: (np.ndarray[ndim=2,dtype=np.float32_t]) : command matrix to use
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef np.ndarray[dtype= np.float32_t] data_F = data.T.copy()

        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.d_cmat.host2device(< float * > data_F.data)
        elif(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.d_cmat.host2device(< float * > data_F.data)
        elif(type_control == scons.ControllerType.GENERIC):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.d_cmat.host2device(< float * > data_F.data)
        else:
            raise TypeError("Controller needs to be ls, mv or generic")

    def set_cmm(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """Set the Cmm matrix on a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index
            data: (np.ndarray[ndim=2,dtype=np.float32_t]) : Cmm matrix
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef np.ndarray[dtype= np.float32_t] data_F = data.T.copy()

        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            controller_mv.d_Cmm.host2device(< float * > data_F.data)
        else:
            raise TypeError("Controller needs to be mv")

    def get_cmm(self, int ncontrol):
        """Return the Cmm matrix from a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            cmm : (np.ndarray[ndim=2,dtype=np.float32_t]) : Cmm matrix
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef const long * cdims

        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_mv.d_Cmm.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_mv.d_Cmm.device2host(< float * > data_F.data)
            return data_F.T.copy()
        else:
            raise TypeError("Controller needs to be mv")

    def get_cphim(self, int ncontrol):
        """Return the Cphim matrix from a sutra_controller_mv object

        :parameters;
            ncontrol: (int) : controller index
        :return:
            cphim : (np.ndarray[ndim=2,dtype=np.float32_t]) : Cphim matrix
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef const long * cdims

        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_mv.d_Cphim.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_mv.d_Cphim.device2host(< float * > data_F.data)
            return data_F.T.copy()
        else:
            raise TypeError("Controller needs to be mv")

    def get_cmat(self, int ncontrol):
        """Return the command matrix from a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            cmat : (np.ndarray[ndim=2,dtype=np.float32_t]) : command matrix
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef const long * cdims
        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_ls.d_cmat.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_ls.d_cmat.device2host(< float * > data_F.data)
            return data_F.T.copy()
        elif(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_mv.d_cmat.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_mv.d_cmat.device2host(< float * > data_F.data)
            return data_F.T.copy()
        elif(type_control == scons.ControllerType.GENERIC):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            cdims = controller_generic.d_cmat.getDims()
            data_F = np.zeros((cdims[2], cdims[1]), dtype=np.float32)
            controller_generic.d_cmat.device2host(< float * > data_F.data)
            return data_F.T.copy()
        else:
            raise TypeError("Controller needs to be ls, mv or generic")

    def set_decayFactor(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] decay):
        """Set the decay factor on a sutra_controller_generic object

        :parameters:
            ncontrol: (int) : controller index
            decay: (np.ndarray[ndim=1,dtype=np.float32_t]) : ask to Rico
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.GENERIC):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.set_decayFactor(< float * > decay.data)
        else:
            raise TypeError("Controller needs to be generic")

    def set_matE(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] matE):
        """Set the matrix E on a sutra_controller_generic object

        :parameters:
            ncontrol: (int) : controller index
            matE: (np.ndarray[ndim=2,dtype=np.float32_t]) : ask to Rico
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        cdef np.ndarray[dtype= np.float32_t] matE_F = matE.T.copy()

        if(type_control == scons.ControllerType.GENERIC):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.set_matE(< float * > matE_F.data)
        else:
            raise TypeError("Controller needs to be generic")

    def set_openloop(self, int ncontrol, int openloop):
        """Set the openloop state to a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index

            openloop: state of the controller
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller * controller = self.rtc.d_control[ncontrol]
        controller.set_openloop(openloop)

    def set_commandlaw(self, int ncontrol, bytes law):
        """Set the law to a sutra_controller_generic object

        :parameters:
            ncontrol: (int) : controller index

            law: law of the controller
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)

        cdef sutra_controller_generic * controller_generic
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.GENERIC):
            controller_generic = dynamic_cast_controller_generic_ptr(
                self.rtc.d_control[ncontrol])
            controller_generic.set_commandlaw(law)
        else:
            raise TypeError("Controller needs to be ls, mv or generic")

    def set_perturbcom(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] perturb):
        """
        TODO: docstring

        :parameters:
            ncontrol: (int) : controller index
            perturb: (np.ndarray[ndim = 2, dtype = np.float32_t]) : perturb voltages

        """
        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        control.set_perturbcom(< float * > perturb.data, perturb.shape[0])

    def set_centroids_ref(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] centroids_ref):
        """
            TODO: docstring
        :parameters:
        ncontrol: (int) : controller index
        """
        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        control.set_centroids_ref(< float * > centroids_ref.data)

    def get_centroids_ref(self, int ncontrol):
        """
            TODO: docstring
        :parameters:
        ncontrol: (int) : controller index
        """
        cdef sutra_controller * control = self.rtc.d_control[ncontrol]
        cdef int nslope = control.nslope()
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] centroids_ref
        centroids_ref = np.zeros(nslope, dtype=np.float32)
        control.get_centroids_ref(< float * > centroids_ref.data)
        return centroids_ref

    def set_U(self, int ncontrol, np.ndarray[ndim=2, dtype=np.float32_t] U):
        """Set the eigen modes matrix of the imat decomposition in a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index
            U: (np.ndarray[ndim=2,dtype=np.float32_t]) : eigen modes matrix
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        cdef np.ndarray[dtype= np.float32_t] data_F = U.T.copy()
        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.d_U.host2device(< float * > data_F.data)

    def set_eigenvals(self, int ncontrol, np.ndarray[ndim=1, dtype=np.float32_t] eigenvals):
        """Set the eigen values of the imat decomposition in a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index
            eigenvals: (np.ndarray[ndim=1,dtype=np.float32_t]) : eigen values
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            controller_ls.h_eigenvals.fill_from(< float * > eigenvals.data)
        else:
            raise TypeError("Controller needs to be ls")

    def get_U(self, int ncontrol):
        """Return the eigen modes matrix of the imat decomposition from a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            U : (np.ndarray[ndim=2,dtype=np.float32_t]) : eigen modes matrix
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef const long * dims = NULL

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.d_U.getDims()
            data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_ls.d_U.device2host(< float * > data_F.data)

        return data_F.T.copy()

    def get_eigenvals(self, int ncontrol):
        """Return the eigen values of the imat decomposition in a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            eigenvals : (np.ndarray[ndim=1,dtype=np.float32_t]) : eigenvalues
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        cdef const long * dims
        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.h_eigenvals.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            controller_ls.h_eigenvals.fill_into(< float * > data.data)
        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_mv.h_Cmmeigenvals.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            controller_mv.h_eigenvals.fill_into(< float * > data.data)

        return data

    def get_cmm_eigenvals(self, int ncontrol):
        """Return the eigen values of the Cmm decomposition in a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            eigenvals : (np.ndarray[ndim=1,dtype=np.float32_t]) : eigenvalues
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_mv * controller_mv
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        cdef const long * dims
        if(type_control == scons.ControllerType.MV):
            controller_mv = dynamic_cast_controller_mv_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_mv.h_Cmmeigenvals.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            controller_mv.h_Cmmeigenvals.fill_into(< float * > data.data)

        return data

    def get_geocov(self, int ncontrol):
        """Return the geocov matrix of the sutra_controller_geo object. In case of error_budget computation, this matrix is Btt basis
        :parameters:
            ncontrol: (int) : controller index
        :return:
            geocov : (np.ndarray[ndim=2,dtype=np.float32_t]) : geocov matrix
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_contro = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef const long * dims
        if(type_contro == scons.ControllerType.GEO):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_geo.d_geocov.getDims()
            data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_geo.d_geocov.device2host(< float*>data_F.data)
        else:
            raise TypeError('controller must be geo')

        return data_F.T.copy()

    def get_IFtt(self, int ncontrol):
        """Return the TT IF matrix of the sutra_controller_geo object (in case of error_budget computation)
        :parameters:
            ncontrol: (int) : controller index
        :return:
            IFtt : (np.ndarray[ndim=2,dtype=np.float32_t]) : IF TT matrix
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_contro = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef const long * dims
        if(type_contro == scons.ControllerType.GEO):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            if(controller_geo.Ntt):
                dims = controller_geo.d_TT.getDims()
                data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
                controller_geo.d_TT.device2host(< float*>data_F.data)
                return data_F.T.copy()
            else:
                raise ValueError("TT not initialized : only with roket")

    def get_cenbuff(self, int ncontrol):
        """Return the centroids buffer from a sutra_controller_ls object.
        This buffer contains centroids from iteration i-delay to current iteration.

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=2,dtype=np.float32_t]) : centroids buffer
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef const long * dims = NULL

        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.d_cenbuff.getDims()
            data_F = np.zeros((dims[2], dims[1]), dtype=np.float32)
            controller_ls.d_cenbuff.device2host(< float * > data_F.data)

        return data_F.T.copy()

    def get_err(self, int ncontrol):
        """Return the command increment (cmat*slopes) from a sutra_controller_ls object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : command increment
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef sutra_controller_ls * controller_ls
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        cdef const long * dims
        if(type_control == scons.ControllerType.LS):
            controller_ls = dynamic_cast_controller_ls_ptr(
                self.rtc.d_control[ncontrol])
            dims = controller_ls.d_err.getDims()
            data = np.zeros((dims[1]), dtype=np.float32)
            controller_ls.d_err.device2host(< float * > data.data)

        return data

    def get_com(self, int ncontrol):
        """Return the command vector from a sutra_controller object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : command vector
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        cdef const long * dims
        dims = self.rtc.d_control[ncontrol].d_com.getDims()
        data = np.zeros((dims[1]), dtype=np.float32)
        self.rtc.d_control[ncontrol].d_com.device2host(< float * > data.data)

        return data

    def get_olmeas(self, int ncontrol):
        """Return the reconstructed open-loop measurement from a sutra_controller_mv object

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : reconstructed open-loop
        """
        cdef sutra_controller_mv * controller_mv
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        cdef const long * dims

        controller_mv = dynamic_cast_controller_mv_ptr(
            self.rtc.d_control[ncontrol])
        dims = controller_mv.d_olmeas.getDims()
        data = np.zeros((dims[1]), dtype=np.float32)
        controller_mv.d_olmeas.device2host(< float * > data.data)

        return data

    def get_voltage(self, int ncontrol):
        """Return the voltage vector that will be effectively applied to the DMs

        :parameters:
            ncontrol: (int) : controller index
        :return:
            data : (np.ndarray[ndim=1,dtype=np.float32_t]) : voltage vector
        """
        self.context.set_activeDeviceForCpy(self.rtc.device, 1)
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        cdef const long * dims
        dims = self.rtc.d_control[ncontrol].d_voltage.getDims()
        data = np.zeros((dims[1]), dtype=np.float32)
        self.rtc.d_control[ncontrol].d_voltage.device2host(< float * > data.data)

        return data

    def get_IFsparse(self, int ncontrol):
        """Get the influence functions matrix computed by the geo controller
        Return a scipy.sparse object which shape is (nactus,Npts in the pupil)

        :parameters:
            ncontrol: (int) : controller index
        :return:
            IF : (scipy.sparse) : influence functions matrix
        """
        cdef sutra_controller_geo * controller_geo
        cdef bytes type_control = self.rtc.d_control[ncontrol].get_type()
        sparse = naga_sparse_obj_Double()
        if(type_control == scons.ControllerType.GEO):
            controller_geo = dynamic_cast_controller_geo_ptr(
                self.rtc.d_control[ncontrol])
            sparse.copy(controller_geo.d_IFsparse)
            return sparse.get_sparse()
        else:
            print(
                "This controller does not have a sparse_obj\n(or this function is not implemented for it yet")


IF USE_BRAMA == 1:
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
            pass  # del self.rtc

        cpdef publish(self):
            cdef sutra_rtc_brama * rtc = < sutra_rtc_brama * > (self.rtc)
            rtc.publish()
