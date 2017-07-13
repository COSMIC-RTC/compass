include "../par.pxi"

import numpy as np
cimport numpy as np
# np.import_array()

import make_pupil as mkP

import os

IF USE_MPI == 1:
    # from mpi4py import MPI
    from mpi4py cimport MPI
    # C-level cdef, typed, Python objects
    # from mpi4py cimport mpi_c as mpi
    from mpi4py cimport libmpi as mpi

IF USE_MPI == 2:
    cimport mpi4py.MPI as MPI
    cimport mpi4py.libmpi as mpi


cdef class Sensors:
    """
        Constructor:
        Sensors(nsensors, type_data, npup, nxsub, nvalid, nphase,
                pdiam, npix, nrebin, nfft, nftota, nphot, lgs, odevice,
                comm_size, rank)

    :parameters:
        nsensors: (int) :
        type_data: list of strings) :
        npup: (np.ndarray[ndim=1,dtype=np.int64_t]) :
        nxsub: (np.ndarray[ndim=1,dtype=np.int64_t]) :
        nvalid: (np.ndarray[ndim=1,dtype=np.int64_t]) :
        nphase: (np.ndarray[ndim=1,dtype=np.int64_t]) :
        pdiam: (np.ndarray[ndim=1,dtype=np.float32_t) :
        npix: (np.ndarray[ndim=1,dtype=np.int64_t]) :
        nrebin: (np.ndarray[ndim=1,dtype=np.int64_t]) :
        nfft: (np.ndarray[ndim=1,dtype=np.int64_t]) :
        ntota: (np.ndarray[ndim=1,dtype=np.int64_t]) :
        nphot: (np.ndarray[ndim=1,dtype=np.float32_t]) :
        nphot4imat: (np.ndarray[ndim=1,dtype=np.float32_t]) :
        lgs: (np.ndarray[ndim=1,dtype=np.int32_t]) :
        odevice: (int) :
        comm_size: (int) : MPI communicator size
        rank: (int) : process rank
    """

    def __cinit__(self, int nsensors, Telescope tel, type_data,
                  np.ndarray[ndim=1, dtype=np.int64_t] npup,
                  np.ndarray[ndim=1, dtype=np.int64_t] nxsub,
                  np.ndarray[ndim=1, dtype=np.int64_t] nvalid,
                  np.ndarray[ndim=1, dtype=np.int64_t] nphase,
                  np.ndarray[ndim=1, dtype=np.float32_t] pdiam,
                  np.ndarray[ndim=1, dtype=np.int64_t] npix=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] nrebin=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] nfft=None,
                  np.ndarray[ndim=1, dtype=np.int64_t] ntota=None,
                  np.ndarray[ndim=1, dtype=np.float32_t] nphot=None,
                  np.ndarray[ndim=1, dtype=np.float32_t] nphot4imat=None,
                  np.ndarray[ndim=1, dtype=np.int32_t] lgs=None,
                  int odevice=-1,
                  int comm_size=1,
                  int rank=0,
                  bool error_budget=False
                  ):

        cdef char ** type_wfs = < char ** > malloc(len(type_data) * sizeof(char *))
        cdef int i
        for i in range(nsensors):
            type_wfs[i] = type_data[i]
        cdef carma_context * context = &carma_context.instance()
        if odevice < 0:
            odevice = context.get_activeDevice()
        if(type_data == b"geo"):
            self.sensors = new sutra_sensors(context, tel.telescope, nsensors,
                                             < long * > nxsub.data,
                                             < long * > nvalid.data,
                                             < long * > nphase.data,
                                             npup[0],
                                             < float * > pdiam.data,
                                             odevice)
        else:
            self.sensors = new sutra_sensors(context, tel.telescope, < char ** > type_wfs, nsensors,
                                             < long * > nxsub.data,
                                             < long * > nvalid.data,
                                             < long * > npix.data,
                                             < long * > nphase.data,
                                             < long * > nrebin.data,
                                             < long * > nfft.data,
                                             < long * > ntota.data,
                                             < long * > npup.data,
                                             < float * > pdiam.data,
                                             < float * > nphot.data,
                                             < float * > nphot4imat.data,
                                             < int * > lgs.data,
                                             odevice, error_budget)

        IF USE_MPI:
            self.sensors.define_mpi_rank(rank, comm_size)

        self.sensors.allocate_buffers()
        self.sensors.device = odevice
        free(type_wfs)

    def __dealloc__(self):
        del self.sensors

    cpdef sensors_initgs(self, np.ndarray[ndim=1, dtype=np.float32_t] xpos,
                         np.ndarray[ndim=1, dtype=np.float32_t] ypos,
                         np.ndarray[ndim=1, dtype=np.float32_t] Lambda,
                         np.ndarray[ndim=1, dtype=np.float32_t] mag,
                         float zerop,
                         np.ndarray[ndim=1, dtype=np.int64_t] size,
                         np.ndarray[ndim=1, dtype=np.float32_t] noise,
                         np.ndarray[ndim=1, dtype=np.int64_t] seed,
                         np.ndarray[ndim=1, dtype=np.float32_t] G,
                         np.ndarray[ndim=1, dtype=np.float32_t] thetaML,
                         np.ndarray[ndim=1, dtype=np.float32_t] dx,
                         np.ndarray[ndim=1, dtype=np.float32_t] dy):
        """
            Call the function sensors_initgs

        :parameters:
            xpos: (np.ndarray[ndim=1,dtype=np.float32_t]) :
            ypos: (np.ndarray[ndim=1,dtype=np.float32_t]) :
            Lambda: (np.ndarray[ndim=1,dtype=np.float32_t]) :
            mag: (np.ndarray[ndim=1,dtype=np.float32_t]) :
            zerop: (float) :
            size:  (np.ndarray[ndim=1,dtype=np.int64_t  ]) :
            noise: (np.ndarray[ndim=1,dtype=np.float32_t]) :
            seed:  (np.ndarray[ndim=1,dtype=np.int64_t  ]) :
            G:  (np.ndarray[ndim=1,dtype=np.float32_t  ]) :
            thetaML:  (np.ndarray[ndim=1,dtype=np.float32_t  ]) :
            dx:  (np.ndarray[ndim=1,dtype=np.float32_t  ]) :
            dy:  (np.ndarray[ndim=1,dtype=np.float32_t  ]) :
        """

        if(noise.size == 0):
            self.sensors.sensors_initgs(< float * > xpos.data, < float * > ypos.data,
                                         < float * > Lambda.data, < float * > mag.data,
                                         zerop, < long * > size.data,
                                         < float * > G.data, < float * > thetaML.data,
                                         < float * > dx.data, < float * > dy.data)

        elif(seed.size == 0):
            self.sensors.sensors_initgs(< float * > xpos.data, < float * > ypos.data,
                                         < float * > Lambda.data, < float * > mag.data,
                                         zerop, < long * > size.data,
                                         < float * > noise.data, < float * > G.data,
                                         < float * > thetaML.data,
                                         < float * > dx.data, < float * > dy.data)
        else:
            self.sensors.sensors_initgs(< float * > xpos.data, < float * > ypos.data,
                                         < float * > Lambda.data, < float * > mag.data,
                                         zerop, < long * > size.data,
                                         < float * > noise.data, < long * > seed.data,
                                         < float * > G.data, < float * > thetaML.data,
                                         < float * > dx.data, < float * > dy.data)

    cpdef sensors_initarr(self, int n, Param_wfs wfs):
        """
            Call the function wfs_initarrays from a sutra_wfs of the Sensors

        :parameters:
            n: (int) : index of the wfs
            wfs: (Param_wfs) :
        """

        cdef np.ndarray tmp_istart, tmp_jstart

        cdef sutra_wfs_geom * wfs_geom = NULL
        cdef sutra_wfs_sh * wfs_sh = NULL
        cdef sutra_wfs_pyr_roof * wfs_roof = NULL
        cdef sutra_wfs_pyr_pyr4 * wfs_pyr = NULL
        cdef sutra_wfs_pyr_pyrhr * wfs_pyrhr = NULL

        cdef np.ndarray[dtype= np.int32_t] phasemap_F = wfs._phasemap.flatten("F")
        cdef np.ndarray[dtype= np.int32_t] hrmap_F = wfs._hrmap.flatten("F")
        cdef np.ndarray[dtype= np.int32_t] binmap_F
        cdef int * binmap = NULL
        if(self.sensors.d_wfs[n].type == b"sh"):
            binmap_F = wfs._binmap.flatten("F")
            binmap = < int * > binmap_F.data

        cdef int * phasemap = < int * > phasemap_F.data
        cdef int * hrmap = < int * > hrmap_F.data
        cdef int * validx = < int * > wfs._validsubsx.data
        cdef int * validy = < int * > wfs._validsubsy.data
        tmp_istart = np.copy(wfs._istart + 1)
        cdef int * istart = < int * > tmp_istart.data
        tmp_jstart = np.copy(wfs._jstart + 1)
        cdef int * jstart = < int * > tmp_jstart.data
        cdef np.ndarray[dtype= np.float32_t] cx_F
        cdef np.ndarray[dtype= np.float32_t] cy_F
        cdef float * cx = NULL
        cdef float * cy = NULL

        if((self.sensors.d_wfs[n].type == b"pyrhr") or (self.sensors.d_wfs[n].type == b"pyr") or (self.sensors.d_wfs[n].type == b"roof")):
            cx_F = wfs._pyr_cx.astype(np.float32)
            cy_F = wfs._pyr_cy.astype(np.float32)
            cx = < float * > cx_F.data
            cy = < float * > cy_F.data

        # [dtype=np.float32_t] halfxy_F=wfs._halfxy.flatten("F")
        cdef np.ndarray halfxy_F
        cdef float * halfxy  # =<float*>halfxy_F.data
        # if(self.sensors.d_wfs[n].type == b"sh"):
        halfxy_F = wfs._halfxy.flatten("F")
        halfxy = < float * > halfxy_F.data

        cdef np.ndarray[dtype= np.float32_t] submask_F
        submask_F = wfs._submask.flatten("F").astype(np.float32)
        cdef np.ndarray[dtype= np.float32_t] sincar_F = wfs._sincar.flatten("F")

        cdef float * submask = < float * > submask_F.data
        cdef float * sincar = < float * > sincar_F.data

        cdef np.ndarray[dtype= np.complex64_t]ftkernel_F
        cdef float * ftkernel = NULL
        if(self.sensors.d_wfs[n].type == b"sh"):
            ftkernel_F = wfs._ftkernel.flatten("F")
            ftkernel = < float * > ftkernel_F.data

        cdef float * fluxPerSub = NULL
        cdef np.ndarray tmp
        cdef np.ndarray tmp_offset
        cdef cuFloatComplex * offset = NULL

        cdef np.ndarray tmp_halfxy
        cdef cuFloatComplex * pyr_halfxy = NULL

        # swap validx and validy due to transposition from yorick to python
        if(self.sensors.d_wfs[n].type == b"geo"):
            tmp = wfs._fluxPerSub.T[np.where(wfs._isvalid > 0)].copy()
            fluxPerSub = < float * > tmp.data
            wfs_geom = dynamic_cast_wfs_geom_ptr(self.sensors.d_wfs[n])
            wfs_geom.wfs_initarrays(phasemap, halfxy, fluxPerSub,
                                    validx, validy)

        elif(self.sensors.d_wfs[n].type == b"pyr"):
            tmp_offset = wfs._pyr_offsets.flatten("F").copy()
            offset = <cuFloatComplex * >tmp_offset.data
            wfs_pyr = dynamic_cast_wfs_pyr_pyr4_ptr(self.sensors.d_wfs[n])
            wfs_pyr.wfs_initarrays(< cuFloatComplex*>halfxy, offset,
                                    submask, cx, cy,
                                    sincar, phasemap, validx, validy)

        elif(self.sensors.d_wfs[n].type == b"pyrhr"):
            tmp = wfs._fluxPerSub.T[np.where(wfs._isvalid > 0)].copy()
            fluxPerSub = < float * > tmp.data
            tmp_halfxy = np.exp(
                1j * wfs._halfxy).astype(np.complex64).flatten("F").copy()
            pyr_halfxy = < cuFloatComplex * > tmp_halfxy.data
            wfs_pyrhr = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
            wfs_pyrhr.wfs_initarrays(
                pyr_halfxy,
                cx,
                cy,
                sincar,
                submask,
                validy,
                validx,
                phasemap,
                fluxPerSub)

        elif(self.sensors.d_wfs[n].type == b"roof"):
            tmp_offset = wfs.__pyr_offsets.flatten("F").copy()
            offset = < cuFloatComplex * > tmp_offset.data
            wfs_roof = dynamic_cast_wfs_pyr_roof_ptr(self.sensors.d_wfs[n])
            wfs_roof.wfs_initarrays(< cuFloatComplex * > halfxy, offset,
                                     submask, cx, cy,
                                     sincar, phasemap, validx, validy)

        elif(self.sensors.d_wfs[n].type == b"sh"):
            tmp = wfs._fluxPerSub.T[np.where(wfs._isvalid > 0)].copy()
            fluxPerSub = < float * > tmp.data
            wfs_sh = dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[n])
            wfs_sh.wfs_initarrays(phasemap, hrmap, binmap, halfxy,
                                  fluxPerSub, validx, validy,
                                  istart, jstart, < cuFloatComplex * > ftkernel)

    cpdef sensors_addlayer(self, int i, bytes type_dm, float alt,
                           float xoff, float yoff):
        """
            Call function add_layer from the sutra_source of a sutra_wfs of the Sensors

        :parameters:
            i: (int) :
            type_dm: (string) :
            alt: (float) :
            xoff: (float) :
            yoff: (float) :
        """

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.sensors.device, 1)
        self.sensors.d_wfs[i].d_gs.add_layer(< char * > type_dm, alt, xoff, yoff)

    def sensors_compimg(self, int n, bool noise=True):
        """
            Compute the wfs image

        :param n: (in) : index of the wfs
        """
        cdef carma_context * context = &carma_context.instance()
        cdef float tmp_noise
        context.set_activeDeviceForCpy(self.sensors.device, 1)
        if(noise):
            self.sensors.d_wfs[n].comp_image()
        else:
            tmp_noise = self.sensors.d_wfs[n].noise
            self.sensors.d_wfs[n].noise = -1.
            self.sensors.d_wfs[n].comp_image()
            self.sensors.d_wfs[n].noise = tmp_noise

    def get_offsets(self, int n):
        """
            Return the 'offset' array of a given wfs

        :param n: (int) : number of the wfs to get the 'offset' from
        """
        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        img = self.sensors.d_wfs[n].d_offsets
        cdims = img.getDims()
        data = np.empty((cdims[1]), dtype=np.float32)
        img.device2host(< float * > data.data)
        return data

    def get_subsum(self, int n):
        """
            Return the 'subsum' array of a given wfs

        :param n: (int) : number of the wfs to get the 'subsum' from
        """
        cdef carma_obj[float] * subsum
        cdef const long * cdims
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data
        subsum = self.sensors.d_wfs[n].d_subsum
        cdims = subsum.getDims()
        data = np.empty((cdims[1]), dtype=np.float32)
        subsum.device2host(< float * > data.data)
        return data

    def get_imgtele(self, int n, Telescope tel=None, Atmos atmos=None, Dms dms=None):
        """
            Return the 'image_telemetry' array of a given wfs

        :param n: (int) : number of the wfs to get the 'image_telemetry' from

        :options for raw image computation
            tel (Telescope) : shesha telescope
            atmos (Atmos) : shesha atmos
            dms (Dms) : shesha dms
        """
        cdef carma_host_obj[float] * img
        cdef sutra_wfs_sh * wfs = dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[n])
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data

        if (tel and atmos and dms):
            self.sensor_trace(n, "all", tel, atmos, dms)
            self.sensor_compimg(n)
        img = self.sensors.d_wfs[n].image_telemetry
        cdims = img.getDims()
        data = np.empty((cdims[1], cdims[2]), dtype=np.float32)

        wfs.fill_binimage(1)
        img.fill_into(< float * > data.data)
        data[np.where(data < 0)] = 0
        return data

    cpdef get_binimg(self, int n, Telescope tel=None, Atmos atmos=None, Dms dms=None):
        """
            Return the 'binimg' array of a given wfs

        :param
            n: (int) :number of the wfs to get the 'binimg' from
        :options for raw image computation
            tel (Telescope) : shesha telescope
            atmos (Atmos) : shesha atmos
            dms (Dms) : shesha dms
        """

        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        if (tel and atmos and dms):
            self.sensor_trace(n, "all", tel, atmos, dms)
            self.sensor_compimg(n)
        img = self.sensors.d_wfs[n].d_binimg
        cdims = img.getDims()
        data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
        data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
        dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[n]).fill_binimage(0)
        img.device2host(< float * > data_F.data)

        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        data[np.where(data < 0)] = 0
        return data

    cpdef get_binimg_notnoisy(self, int n):
        """
            Return the 'binimg_notnoisy' array of a given pyrhr wfs

        :param
            n: (int) :number of the wfs to get the 'binimg_notnoisy' from
        """
        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F

        if(self.sensors.error_budget):
            img = self.sensors.d_wfs[n].d_binimg_notnoisy
            cdims = img.getDims()
            data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
            data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
            img.device2host(< float * > data_F.data)
            data = np.reshape(
                data_F.flatten("F"), (cdims[1], cdims[2]))
            return data
        else:
            raise TypeError("the error budget analysis has to be enabled")

    def get_pyrimg(self, int n):
        """
            Return the image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from

        """

        return self._get_pyrimg(n)

    cdef _get_pyrimg(self, int n):
        """
            Return the image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from

        """
        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F

        cdef np.ndarray[ndim = 3, dtype = np.float32_t] bincube
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] pyrimg
        cdef bytes type_wfs = self.sensors.d_wfs[n].type
        cdef int npix

        if(type_wfs == b"pyrhr"):
            img = self.sensors.d_wfs[n].d_binimg
            cdims = img.getDims()
            data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
            data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
            img.device2host(< float * > data_F.data)
            data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return data

        if(type_wfs == b"pyr"):
            bincube = self.get_bincube(n)
            npix = bincube.shape[1]
            pyrimg = np.zeros((2 * npix + 3, 2 * npix + 3), dtype=np.float32)
            pyrimg[1:npix + 1, 1:npix + 1] = bincube[:, :, 0]
            pyrimg[npix + 2:2 * npix + 2, 1:npix + 1] = bincube[:, :, 1]
            pyrimg[1:npix + 1, npix + 2:2 * npix + 2] = bincube[:, :, 2]
            pyrimg[npix + 2:2 * npix + 2, npix +
                   2:2 * npix + 2] = bincube[:, :, 3]

            return pyrimg

        else:
            raise TypeError("wfs should be a pyr or pyrhr")

    def set_pyrimg(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            Set the image of a pyr wfsr

        :param n: (int) : number of the wfs to get the image from

        """
        return self._set_pyrimg(n, data)

    cdef _set_pyrimg(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            Set the image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from

        """
        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[dtype= np.float32_t] data_F = data.flatten("F")

        cdef bytes type_wfs = self.sensors.d_wfs[n].type

        if(type_wfs == b"pyrhr"):
            img = self.sensors.d_wfs[n].d_binimg
            img.host2device(< float * > data_F.data)
        else:
            raise TypeError("wfs should be a pyr")

    def set_submask(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            Set the field stop of a pyrhr wfs

        :parameters:
            n: (int) : WFS index
            data: (np.ndarray[ndim=2, dtype=np.float32_t]) : field stop
        """
        return self._set_submask(n, data)

    cdef _set_submask(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):

        cdef bytes type_wfs = self.sensors.d_wfs[n].type
        cdef np.ndarray[dtype= np.float32_t] data_F = data.flatten("F")
        cdef sutra_wfs_pyr_pyrhr * wfs

        if(type_wfs == b"pyrhr"):
            wfs = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
            wfs.set_submask(< float*> data_F.data)
        else:
            raise TypeError("WFS should be a pyrhr for using this function")

    def get_submask(self, int n):
        """
            Get the field stop of a pyrhr wfs

        :parameters:
            n: (int) : WFS index
        """
        return self._get_submask(n)

    cdef _get_submask(self, int n):
        cdef carma_obj[float] * submask
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        cdef sutra_wfs_pyr_pyrhr * wfs
        cdef bytes type_wfs = self.sensors.d_wfs[n].type

        if(type_wfs == b"pyrhr"):
            wfs = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
            submask = wfs.d_submask
            cdims = submask.getDims()
            data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
            data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
            submask.device2host(< float * > data_F.data)
            data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return data

        else:
            raise TypeError("WFS should be a pyrhr for using this function")

    def get_pyrimghr(self, int n):
        """
            Return the high res image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from
        """
        return self._get_pyrimghr(n)

    cdef _get_pyrimghr(self, int n):

        """
            Return the high res image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from
        """
        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F

        cdef bytes type_wfs = self.sensors.d_wfs[n].type
        cdef sutra_wfs_pyr_pyrhr * wfs

        if(type_wfs == b"pyrhr"):
            wfs = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
            img = wfs.d_hrimg
            cdims = img.getDims()
            data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
            data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
            img.device2host(< float * > data_F.data)
            data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return data

        else:
            raise TypeError("wfs should be a pyrhr")

    cpdef set_ncpa_phase(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data_F = data.flatten("F")
        self.sensors.d_wfs[n].set_ncpa_phase(< float*>data_F.data, data.size)

    cpdef get_ncpa_phase(self, int n):
        cdef carma_obj[float] * ph
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data

        ph = self.sensors.d_wfs[n].d_gs.d_phase.d_screen
        cdims = ph.getDims()
        data = np.empty((cdims[2], cdims[1]), dtype=np.float32)
        self.sensors.d_wfs[n].get_ncpa_phase(< float*>data.data, data.size)
        return np.reshape(data.flatten("F"), (cdims[1], cdims[2]))

    cpdef comp_modulation(self, int n, int cpt):
        """
            Return the high res image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from
        """

        cdef bytes type_wfs = self.sensors.d_wfs[n].type
        cdef sutra_wfs_pyr_pyrhr * wfs

        if(type_wfs == b"pyrhr"):
            wfs = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
            wfs.comp_modulation(cpt)

        else:
            raise TypeError("wfs should be a pyrhr")

    cdef _get_bincube(self, int n):
        """
            Return the 'bincube' array of a given wfs

        :param n: (int) : number of the wfs to get the 'bincube' from
        """
        cdef carma_obj[float] * cube
        cdef const long * cdims
        cdef np.ndarray[ndim= 3, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 3, dtype = np.float32_t] data_F
        cube = self.sensors.d_wfs[n].d_bincube
        cdims = cube.getDims()
        data = np.empty((cdims[1], cdims[2], cdims[3]), dtype=np.float32)
        data_F = np.empty((cdims[3], cdims[2], cdims[1]), dtype=np.float32)
        cube.device2host(< float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
        return data

    def get_bincube(self, int n):
        """
            Return the 'bincube' array of a given wfs

        :param n: (int) : number of the wfs to get the 'bincube' from
        """
        return self._get_bincube(n)

    cpdef set_bincube(self, int n, np.ndarray[ndim=3, dtype=np.float32_t] data):
        """
            Set the bincube of the WFS numner n

        :parameters:
            n: (int) : WFS number
            data: (np.ndarray[ndim=3,dtype=np.float32_t]) : bincube to use
        """
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDeviceForCpy(self.sensors.device, 1)
        cdef np.ndarray[dtype= np.float32_t] data_F = data.flatten("F")

        self.sensors.d_wfs[n].d_bincube.host2device(< float * > data_F.data)

    cpdef get_bincubeNotNoisy(self, int n):
        """
            Return the 'bincube_not_noisy' array of a given wfs. It's the bincube
            before noise has been added

        :param n: (int) : number of the wfs to get the 'bincube_not_noisy' from
        """
        cdef carma_obj[float] * cube
        cdef const long * cdims
        cdef np.ndarray[ndim= 3, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 3, dtype = np.float32_t] data_F
        if(self.sensors.error_budget):
            cube = self.sensors.d_wfs[n].d_bincube_notnoisy
            cdims = cube.getDims()
            data = np.empty((cdims[1], cdims[2], cdims[3]), dtype=np.float32)
            data_F = np.empty((cdims[3], cdims[2], cdims[1]), dtype=np.float32)
            cube.device2host(< float * > data_F.data)
            data = np.reshape(
                data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
            return data
        else:
            raise TypeError("the error budget analysis has to be enabled")

    def reset_phase(self, int n):
        """
            Reset the phase's array of a given wfs

        :param n: (int) : index of the given wfs
        """
        cdef carma_obj[float] * phase
        phase = self.sensors.d_wfs[n].d_gs.d_phase.d_screen
        phase.reset()

    def get_phase(self, int n):
        """
            Return the phase array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        """
        cdef carma_obj[float] * phase
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        phase = self.sensors.d_wfs[n].d_gs.d_phase.d_screen
        cdims = phase.getDims()
        data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
        data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
        phase.device2host(< float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        return data

    def get_lgskern(self, int n):
        """
            Return the lgskern array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        """
        cdef carma_obj[float] * lgs_kern
        cdef const long * cdims
        cdef np.ndarray[ndim= 3, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 3, dtype = np.float32_t] data_F
        if(self.sensors.d_wfs[n].lgs):
            lgs_kern = self.sensors.d_wfs[n].d_gs.d_lgs.d_lgskern
            cdims = lgs_kern.getDims()
            data = np.empty((cdims[1], cdims[2], cdims[3]), dtype=np.float32)
            data_F = np.empty((cdims[3], cdims[2], cdims[1]), dtype=np.float32)
            lgs_kern.device2host(< float * > data_F.data)
            data = np.reshape(
                data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
            return data
        else:
            raise TypeError("the WFS should be a LGS")

    def get_ftlgskern(self, int n):
        """
            Return the ftlgskern array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        """
        cdef carma_obj[cuFloatComplex] * ftlgs_kern
        cdef const long * cdims
        cdef np.ndarray[ndim= 3, dtype = np.complex64_t] data
        cdef np.ndarray[ndim= 3, dtype = np.complex64_t] data_F
        if(self.sensors.d_wfs[n].lgs):
            ftlgs_kern = self.sensors.d_wfs[n].d_gs.d_lgs.d_ftlgskern
            cdims = ftlgs_kern.getDims()
            data = np.empty((cdims[1], cdims[2], cdims[3]), dtype=np.complex64)
            data_F = np.empty(
                (cdims[3], cdims[2], cdims[1]), dtype=np.complex64)
            ftlgs_kern.device2host(< cuFloatComplex * > data_F.data)
            data = np.reshape(
                data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
            return data
        else:
            raise TypeError("the WFS should be a LGS")

    def set_phase(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """
            Set the phase array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        :param data: (np.ndarray) : the phase to set
        """
        # self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.sensors.d_wfs[n].d_gs

        cdef np.ndarray[dtype= np.float32_t] data_F = data.flatten("F")

        src.d_phase.d_screen.host2device(< float * > data_F.data)

    def comp_new_fstop(self, int n, Param_wfs wfs, float fssize, bytes fstop):
        """
            Compute a new field stop for pyrhr WFS

        :parameters:
            n : (int) : WFS index
            wfs : (Param_wfs) : WFS parameters
            fssize : (float) : field stop size [arcsec]
            fstop : (string) : "square" or "round" (field stop shape)
        """
        fsradius_pixels = long(fssize / wfs._qpixsize / 2.)
        if (fstop == b"round"):
            wfs.fstop = fstop
            focmask = mkP.dist(
                wfs._Nfft, xc=wfs._Nfft / 2. + 0.5, yc=wfs._Nfft / 2. + 0.5) < (fsradius_pixels)
            # fstop_area = np.pi * (wfs.fssize/2.)**2. #UNUSED
        elif (wfs.fstop == b"square"):
            wfs.fstop = fstop
            x, y = indices(wfs._Nfft)
            x -= (wfs._Nfft + 1.) / 2.
            y -= (wfs._Nfft + 1.) / 2.
            focmask = (np.abs(x) <= (fsradius_pixels)) * \
                (np.abs(y) <= (fsradius_pixels))
            # fstop_area = wfs.fssize**2. #UNUSED
        else:
            msg = "wfs " + str(n) + ". fstop must be round or square"
            raise ValueError(msg)

        # pyr_focmask = np.roll(focmask,focmask.shape[0]/2,axis=0)
        # pyr_focmask = np.roll(pyr_focmask,focmask.shape[1]/2,axis=1)
        pyr_focmask = focmask * 1.0  # np.fft.fftshift(focmask*1.0)
        wfs._submask = np.fft.fftshift(pyr_focmask).astype(np.float32)
        wfs_fssize = fssize
        self._set_submask(n, wfs._submask)

    def get_camplipup(self, int n):
        """
            Return the 'camplipup' array of a given wfs

        :param n: (int) : number of the wfs to get the 'camplipup' from
        """
        cdef carma_obj[cuFloatComplex] * amplipup
        cdef const long * cdims

        cdef np.ndarray[ndim= 3, dtype = np.complex64_t] data
        cdef np.ndarray[ndim= 3, dtype = np.complex64_t] data_F
        amplipup = self.sensors.d_wfs[n].d_camplipup
        cdims = amplipup.getDims()

        data = np.zeros((cdims[1], cdims[2], cdims[3]), dtype=np.complex64)
        data_F = np.zeros((cdims[3], cdims[2], cdims[1]), dtype=np.complex64)
        amplipup.device2host(< cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))

        return data

    def get_camplipup_pyr(self, int n):
        """
            Return the 'camplipup' array of a given wfs in the pyr case

        :param n: (int) : number of the wfs to get the 'camplipup' from

        """
        cdef carma_obj[cuFloatComplex] * amplipup
        cdef const long * cdims

        cdef np.ndarray[ndim= 2, dtype = np.complex64_t] data
        cdef np.ndarray[ndim= 2, dtype = np.complex64_t] data_F
        amplipup = self.sensors.d_wfs[n].d_camplipup
        cdims = amplipup.getDims()

        data = np.zeros((cdims[1], cdims[2]), dtype=np.complex64)
        data_F = np.zeros((cdims[2], cdims[1]), dtype=np.complex64)
        amplipup.device2host(< cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))

        return data

    def get_amplifoc(self, int n):
        """
            Return the 'amplifoc' array of a given wfs

        :param n: (int) : number of the wfs to get the 'amplifoc' from
        """
        cdef carma_obj[cuFloatComplex] * amplifoc
        cdef const long * cdims
        cdef np.ndarray[ndim= 3, dtype = np.complex64_t] data
        cdef np.ndarray[ndim= 3, dtype = np.complex64_t] data_F
        amplifoc = self.sensors.d_wfs[n].d_camplifoc
        cdims = amplifoc.getDims()
        data = np.zeros((cdims[1], cdims[2], cdims[3]), dtype=np.complex64)
        data_F = np.zeros((cdims[3], cdims[2], cdims[1]), dtype=np.complex64)
        amplifoc.device2host(< cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
        return data

    def get_amplifoc_pyr(self, int n):
        """
            Return the 'amplifoc' array of a given wfs

        :param n: (int) : number of the wfs to get the 'amplifoc' from
        """
        cdef carma_obj[cuFloatComplex] * amplifoc
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.complex64_t] data
        cdef np.ndarray[ndim= 2, dtype = np.complex64_t] data_F
        amplifoc = self.sensors.d_wfs[n].d_camplifoc
        cdims = amplifoc.getDims()
        data = np.zeros((cdims[1], cdims[2]), dtype=np.complex64)
        data_F = np.zeros((cdims[2], cdims[1]), dtype=np.complex64)
        amplifoc.device2host(< cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        return data

    def get_fttotim_pyr(self, int n):
        """
            Return the 'fttotim' array of a given wfs

        :param n: (int) : number of the wfs to get the 'amplifoc' from
        """
        cdef carma_obj[cuFloatComplex] * fttotim
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.complex64_t] data
        cdef np.ndarray[ndim= 2, dtype = np.complex64_t] data_F
        fttotim = self.sensors.d_wfs[n].d_fttotim
        cdims = fttotim.getDims()
        data = np.zeros((cdims[1], cdims[2]), dtype=np.complex64)
        data_F = np.zeros((cdims[2], cdims[1]), dtype=np.complex64)
        fttotim.device2host(< cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        return data

    def get_hrimg_pyr(self, int n):
        """
            Return the phase array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        """
        cdef carma_obj[float] * hrimg
        cdef const long * cdims
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim= 2, dtype = np.float32_t] data_F
        wfs = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
        hrimg = wfs.d_hrimg
        # hrimg=self.sensors.d_wfs[n].d_hrimg
        cdims = hrimg.getDims()
        data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
        data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
        hrimg.device2host(< float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        return data

    cpdef _get_slopesDims(self, int n):
        """
            Return the dimension of the slopes array of a given wfs

        :param n: (int) : number of the wfs to get the 'slopes' dimension from
        """
        cdef int comm_size, rank
        IF USE_MPI:
            mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD, & comm_size)
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, & rank)
        ELSE:
            comm_size = 1
            rank = 0

        cdef const long * cdims
        cdef long dim_tot
        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.sensors.device, 1)
        cdims = self.sensors.d_wfs[n].d_slopes.getDims()
        dim_tot = cdims[1]

        IF USE_MPI:
            mpi.MPI_Allreduce(mpi.MPI_IN_PLACE, & dim_tot, 1, mpi.MPI_LONG, mpi.MPI_SUM, mpi.MPI_COMM_WORLD)
        return dim_tot

    def get_slopes(self, int n):
        """
            Return the 'slopes' array of a given wfs

        :param n: (int) : number of the wfs to get the 'slopes' from
        """

        return self._get_slopes(n)

    cdef _get_slopes(self, int n):
        """
            Return the 'slopes' array of a given wfs

        :param n: (int) : number of the wfs to get the 'slopes' from
        """
        cdef carma_obj[float] * slopes
        cdef const long * cdims
        cdef np.ndarray[ndim= 1, dtype = np.float32_t] data

        cdef carma_context * context = &carma_context.instance()
        context.set_activeDevice(self.sensors.device, 1)

        slopes = self.sensors.d_wfs[n].d_slopes
        cdims = slopes.getDims()
        data = np.empty((cdims[1]), dtype=np.float32)
        slopes.device2host(< float * > data.data)

        IF USE_MPI:
            cdef int comm_size, rank
            mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD, & comm_size)
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, & rank)

            cdef int d = < int > (cdims[1] / 2)

            cdef int * count = < int * > malloc(comm_size * sizeof(int))
            mpi.MPI_Allgather(& d, 1, mpi.MPI_INT, count, 1, mpi.MPI_INT, mpi.MPI_COMM_WORLD)

            cdef int * disp = < int * > malloc((comm_size + 1) * sizeof(int))
            cdef int i, nvalid2
            disp[0] = 0
            for i in range(comm_size):
                disp[i + 1] = disp[i] + count[i]

            cdef np.ndarray[ndim= 1, dtype = np.float32_t] all_slopes = np.zeros(
                disp[comm_size] * 2, dtype=np.float32)

            cdef float * send = < float * > data.data
            cdef float * recv = < float * > all_slopes.data

            mpi.MPI_Allgatherv(send, count[rank], mpi.MPI_FLOAT,
                               recv, count, disp,
                               mpi.MPI_FLOAT, mpi.MPI_COMM_WORLD)

            mpi.MPI_Allgatherv(& send[count[rank]], count[rank], mpi.MPI_FLOAT,
                                & recv[disp[comm_size]], count, disp,
                                mpi.MPI_FLOAT, mpi.MPI_COMM_WORLD)

            return all_slopes

        ELSE:
            return data

    cpdef slopes_geom(self, int nsensor, int t):
        """
            Compute the geometric slopes in a sutra_wfs object

        :parameters:
            nsensor: (int) : wfs number

            :param t: (int) : method (0 or 1)
        """
        cdef sutra_wfs_sh * wfs_sh = NULL
        cdef sutra_wfs_pyr_pyrhr * wfs_pyrhr = NULL

        if(self.sensors.d_wfs[nsensor].type == b"sh"):
            #raise TypeError("wfs should be a SH")
            wfs_sh = dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[nsensor])
            wfs_sh.slopes_geom(t)
        else:
            if(self.sensors.d_wfs[nsensor].type == b"pyrhr"):
                wfs_pyrhr = dynamic_cast_wfs_pyr_pyrhr_ptr(
                    self.sensors.d_wfs[nsensor])
                wfs_pyrhr.slopes_geom(t)
            else:
                raise TypeError("wfs should be a SH or PYRHR")

    cpdef sensors_trace(self, int n, bytes type_trace, Telescope tel=None, Atmos atmos=None, Dms dms=None, int rst=0, int ncpa=0):
        """
            Does the raytracing for the wfs phase screen in sutra_wfs

        :parameters:
            n: (int) :
            type_trace: (str) : "all" : raytracing across atmos and dms seen
                                "dm"  : raytracing across dms seen only
                                "atmos" : raytracing across atmos only
                                "none" : raytracing across tel pupil only
            tel: (Telescope) :(optional) Telescope object
            atmos: (Atmos) :(optional) Atmos object
            dms: (Dms) : (optional) Dms object
            rst: (int) : (optional) reset before raytracing if rst = 1
            ncpa: (int) : (optional) NCPA raytracing if ncpa = 1
        """

        cdef carma_context * context = &carma_context.instance()
        cdef carma_obj[float] * d_screen = self.sensors.d_wfs[n].d_gs.d_phase.d_screen
        context.set_activeDeviceForce(self.sensors.device, 1)

        if(type_trace == b"all"):
            self.sensors.d_wfs[n].sensor_trace(atmos.s_a, dms.dms)
            rst = 0
        elif(type_trace == b"atmos"):
            self.sensors.d_wfs[n].sensor_trace(atmos.s_a)
            rst = 0
        elif(type_trace == b"dm"):
            self.sensors.d_wfs[n].sensor_trace(dms.dms, rst)
            rst = 0

        if ncpa:
            self.sensors.d_wfs[n].sensor_trace(rst)

        if tel is not None:
            d_screen.axpy(1.0, tel.telescope.d_phase_ab_M1_m, 1, 1)

    IF USE_MPI:
        cpdef Bcast_dscreen(self):
            """
                Broadcast the screen of every wfs on process 0 to all process
            """
            cdef carma_obj[float] * screen
            cdef float * ptr
            cdef int i, nsensors, size_i
            cdef long size

            nsensors = self.sensors.nsensors()
            for i in range(nsensors):
                screen = self.sensors.d_wfs[i].d_gs.d_phase.d_screen
                size = screen.getNbElem()
                size_i = size * 2
                ptr = < float * > malloc(size_i * sizeof(float))
                if(self._get_rank(0) == 0):
                    screen.device2host(ptr)
                mpi.MPI_Bcast(
                    ptr, size_i, mpi.MPI_FLOAT, 0, mpi.MPI_COMM_WORLD)
                screen.host2device(ptr)
                free(ptr)

        cpdef Bcast_dscreen_cuda_aware(self):
            """
                Broadcast the screen of every wfs on process 0 to all process
                using cuda_aware
            """
            cdef float * ptr
            cdef int i, nsensors, size_i
            cdef long size

            nsensors = self.sensors.nsensors()
            for i in range(nsensors):
                size = self.sensors.d_wfs[i].d_gs.d_phase.d_screen.getNbElem()
                size_i = size  # convert from long to int
                ptr = self.sensors.d_wfs[i].d_gs.d_phase.d_screen.getData()
                mpi.MPI_Bcast(
                    ptr, size_i, mpi.MPI_FLOAT, 0, mpi.MPI_COMM_WORLD)

        cpdef gather_bincube(self, int n):
            """
                Gather the carma object 'bincube' of a wfs on the process 0

            :parameters:
                comm: (MPI.Intracomm) :communicator mpi
                n: (int) : number of the wfs where the gather will occured
            """

            cdef int * count_bincube = self.sensors.d_wfs[n].count_bincube
            cdef int * displ_bincube = self.sensors.d_wfs[n].displ_bincube
            cdef const long * cdims = self.sensors.d_wfs[n].d_bincube.getDims()

            cdef float * ptr
            ptr = < float * > malloc(cdims[1] * cdims[2] * cdims[3] * sizeof(float))
            self.sensors.d_wfs[n].d_bincube.device2host(ptr)

            cdef int i, j

            if(self._get_rank(n) == 0):
                mpi.MPI_Gatherv(mpi.MPI_IN_PLACE, count_bincube[0], mpi.MPI_FLOAT,
                                ptr, count_bincube, displ_bincube, mpi.MPI_FLOAT,
                                0, mpi.MPI_COMM_WORLD)

            else:
                mpi.MPI_Gatherv(ptr, count_bincube[self.get_rank(n)], mpi.MPI_FLOAT,
                                ptr, count_bincube, displ_bincube, mpi.MPI_FLOAT,
                                0, mpi.MPI_COMM_WORLD)

            if(self._get_rank(n) == 0):
                self.sensors.d_wfs[n].d_bincube.host2device(ptr)

            free(ptr)

        cpdef gather_bincube_cuda_aware(self, int n):
            """
                Gather the carma object 'bincube' of a wfs on the process 0
                using mpi cuda_aware

            :param n: (int) : number of the wfs where the gather will occured
            """
            cdef float * recv_bin = self.sensors.d_wfs[n].d_bincube.getData()
            cdef float * send_bin = self.sensors.d_wfs[n].d_bincube.getData()
            cdef int nx = self.sensors.d_wfs[n].npix
            cdef int nz = self.sensors.d_wfs[n].nvalid
            cdef int * count_bincube = self.sensors.d_wfs[n].count_bincube
            cdef int * displ_bincube = self.sensors.d_wfs[n].displ_bincube

            mpi.MPI_Gatherv(send_bin, nx * nx * nz, mpi.MPI_FLOAT,
                            recv_bin, count_bincube, displ_bincube, mpi.MPI_FLOAT,
                            0, mpi.MPI_COMM_WORLD)

    cdef _get_rank(self, int n):
        """
            Return the rank of one of the sensors wfs

        :param n: (int): index of the wfs to get the rank for
        """
        IF USE_MPI:
            return self.sensors.d_wfs[n].rank
        ELSE:
            return 0

    def get_rank(self, int n):
        """Return the rank of one of the sensors wfs

        :param n: (int) : index of the wfs to get the rank for
        """
        return self.sensors.d_wfs[n].rank

    def __str__(self):
        info = "Sensors object:\n"
        info += "Contains " + str(self.sensors.nsensors()) + " WFS(s):\n"
        info += "WFS # |  Nsubaps  | Nvalid | Npix | Nphase | Nfft | Nrebin | Ntot | Npup\n"
        cdef int i
        cdef sutra_wfs * wfs
        for i in range(< int > self.sensors.nsensors()):
            wfs = self.sensors.d_wfs.at(i)
            info += "%5d" % (i + 1) + " | " + "%3d" % wfs.nxsub + " x " + "%-3d" % wfs.nxsub + " | "\
                "%6d" % wfs.nvalid + " | " + "%4d" % wfs.npix + " | " + "%6d" % wfs.nphase + " | " + \
                "%4d" % wfs.nfft + " | " + "%6d" % wfs.nrebin + " | " + \
                    "%4d" % wfs.ntot + " | " + "%4d" % wfs.npup + "\n"

        info += "--------------------------------------------------------"
        return info

    '''
    #for profiling purpose
    @cython.profile(True)
    cdef gather_bincube_prof(self,int n):
        cdef int nx=self.sensors.d_wfs[n].npix
        cdef int ny=nx
        cdef int nz=self.sensors.d_wfs[n].nvalid
        cdef int nz_t=self.sensors.d_wfs[n].nvalid_tot
        cdef int size=nx*ny*nz
        cdef int *count_bincube=self.sensors.d_wfs[n].count_bincube
        cdef int *displ_bincube=self.sensors.d_wfs[n].displ_bincube
        cdef const long *cdims=self.sensors.d_wfs[n].d_bincube.getDims()

        cdef float *ptr
        ptr=<float*>malloc(cdims[1]*cdims[2]*cdims[3]*sizeof(float))

        self.wait1_prof()
        self.d2h_prof(ptr,n)
        cdef int i,j

        self.wait2_prof()
        self.gather_prof( ptr, size, count_bincube, displ_bincube)

        if(self._get_rank(n)==0):
            self.h2d_prof(ptr,n)
        free(ptr)

    @cython.profile(True)
    cdef wait1_prof(self):
        mpi.MPI_Barrier(mpi.MPI_COMM_WORLD)

    @cython.profile(True)
    cdef wait2_prof(self):
        mpi.MPI_Barrier(mpi.MPI_COMM_WORLD)

    @cython.profile(True)
    cdef d2h_prof(self,float* ptr,n):
        self.sensors.d_wfs[n].d_bincube.device2host(ptr)

    @cython.profile(True)
    cdef h2d_prof(self,float* ptr,n):
        self.sensors.d_wfs[n].d_bincube.host2device(ptr)

    @cython.profile(True)
    cdef gather_prof(self,float *ptr, int size, int *count, int  *displ):
        mpi.MPI_Gatherv(ptr,size,mpi.MPI_FLOAT,
                    ptr, count , displ ,
                    mpi.MPI_FLOAT,0,mpi.MPI_COMM_WORLD)

    '''

    cdef _get_hrmap(self, int n):
        """
            Return the 'bincube' array of a given wfs

        :param n: (int) : number of the wfs to get the 'bincube' from
        """
        cdef carma_obj[int] * cube
        cdef const long * cdims
        cdef np.ndarray[ndim= 1, dtype = np.int32_t] data
        cube = self.sensors.d_wfs[n].d_hrmap
        cdims = cube.getDims()
        data = np.empty((cdims[1]), dtype=np.float32)
        cube.device2host(< int * > data.data)
        return data

    cpdef set_noise(self, int nwfs, np.float32_t noise, np.int64_t seed=-1):
        """
            Set the noise of a given wfs

        :param nwfs: (int) : number of the wfs wanted
        :param noise: (np.float32_t) : desired noise : < 0 = no noise
                                                         0 = photon only
                                                       > 0 = photon + ron in ?
        :param seed: (long) : seed of the new noise
        TODO: fill the noise unit
        """
        if seed < 0:
            seed = 1234 * nwfs
        self.sensors.set_noise(nwfs, noise, seed)

"""
    cdef getDims(self):
        cdef const long *dims_ampli
        cdef const long *dims_fttot

        dims_ampli=self.sensors.d_camplifoc.getDims()
        dims_fttot=self.sensors.d_fttotim.getDims()

        d_amp=np.array([dims_ampli[0],dims_ampli[1],dims_ampli[2],dims_ampli[3]])
        d_tot=np.array([dims_fttot[0],dims_fttot[1],dims_fttot[2],dims_fttot[3]])

        return d_amp,d_tot
"""
