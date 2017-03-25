
import numpy as np
cimport numpy as np
# np.import_array()

import h5py

import hdf5_utils as h5u

import pandas

import os
import iterkolmo as itK

from cython.operator cimport dereference as deref, preincrement as inc
from subprocess import check_output


def atmos_init(naga_context c, Param_atmos atm, Param_tel tel, Param_geom geom,
               Param_loop loop, wfss=None, Param_target target=None,
               int rank=0, int clean=1, dict load={}):
    """atmos_init(naga_context c, Param_atmos atm, Param_tel tel,  Param_geom geom,
                  Param_loop loop, wfss=None, Param_target target=None,
                  int rank=0, int clean=1, dict load={})

    Create and initialise an atmos object

    :parameters:
        c: (naga_context) : context

        tel: (Param_tel) : telescope settings

        geom: (Param_geom) : geometry settings

        loop: (Param_loop) : loop settings

        wfss: (list of Param_wfs) : (optional) wfs settings

        sensors: (Sensors) : (optional) Sensors object on GPU

        target: (Param_target) : (optional) target_settings

        overwrite: (int) : (optional) overwrite data files if overwite=1 (default 1)

        rank: (int) : (optional) rank of the process (default=0)

        clean: (clean) : (optional) ? (default=1)

        load: (dict) : (optional) ? (default={})
    """

    cdef double ittime = loop.ittime
    if(atm.r0 is None):
        atm.r0 = 0

    cdef int i
    cdef float dtor = np.pi / 180

    # ajust layers alt using zenith angle
    atm.alt = (atm.alt / np.cos(geom.zenithangle * dtor))
    # pixel size in meter
    atm.pupixsize = tel.diam / geom.pupdiam

    cdef long max_size
    if (wfss is not None):
        norms = [
            np.linalg.norm([wfss[i].xpos, wfss[i].ypos], axis=0) for i in range(len(wfss))]
        indmax = np.where(norms == np.max(norms))[0][0]
        wfs = wfss[indmax]
    else:
        wfs = None

    # compute total fov using targets and wfs gs
    if ((wfs is not None) and (target is not None)):
        max_size = np.max([np.max(np.linalg.norm([target.xpos, target.ypos], axis=0)),
                           np.max(np.linalg.norm([wfs.xpos, wfs.ypos], axis=0))])
    elif (target is not None):
        max_size = np.max(np.linalg.norm([target.xpos, target.ypos], axis=0))
    elif (wfs is not None):
        max_size = np.max(np.linalg.norm([wfs.xpos, wfs.ypos], axis=0))
    else:
        max_size = 0

    # compute corresponding meta-pupil diameters at each alt
    cdef np.ndarray patch_diam
    patch_diam = geom._n + 2 * \
        (max_size * 4.84814e-6 * atm.alt) / atm.pupixsize + 4
    atm.dim_screens = (patch_diam + patch_diam % 2).astype(np.int64)

    # compute phase screens speed in pixels / iteration
    atm.deltax = geom.pupdiam / tel.diam * atm.windspeed * \
        np.cos(dtor * geom.zenithangle) * ittime
    atm.deltay = atm.deltax * np.sin(dtor * (atm.winddir) + np.pi)
    atm.deltax = atm.deltax * np.cos(dtor * (atm.winddir) + np.pi)

    if(atm.frac.size == 1):
        atm.frac = np.ones((1), dtype=np.float32)
    else:
        atm.frac /= sum(atm.frac)

    if (atm.L0 is None):
        atm.L0 = np.ones(atm.nscreens, dtype=np.float32) * \
            (1.e5)  # infinite L0
        L0_pix = np.copy(atm.L0)
    else:
        if (atm.L0.shape[0] == 1):
            atm.L0 = np.ones(atm.nscreens, dtype=np.float32) * atm.L0[0]
        L0_pix = np.copy(atm.L0)
        for i in range(atm.nscreens):
            frac_l0 = tel.diam / atm.L0[i]
            L0_pix[i] = geom.pupdiam / frac_l0

    if(atm.seeds is None):
        seeds = (np.arange(atm.nscreens, dtype=np.int64) + 1) * 1234
    else:
        seeds = atm.seeds

    return atmos_create(c, atm.nscreens, atm.r0, L0_pix, atm.pupixsize,
                        atm.dim_screens, atm.frac, atm.alt, atm.windspeed,
                        atm.winddir, atm.deltax, atm.deltay, seeds, rank, clean, load)


#################################################
# P-Class atmos
#################################################
cdef class Atmos:
    def __cinit__(self):
        self.context = None

    def __dealloc__(self):
        if(self.s_a != NULL):
            del self.s_a

    cdef realinit(self, naga_context ctxt, int nscreens,
                  np.ndarray[ndim=1, dtype=np.float32_t] r0,
                  np.ndarray[ndim=1, dtype=np.int64_t] size,
                  np.ndarray[ndim=1, dtype=np.float32_t] altitude,
                  np.ndarray[ndim=1, dtype=np.float32_t] windspeed,
                  np.ndarray[ndim=1, dtype=np.float32_t] winddir,
                  np.ndarray[ndim=1, dtype=np.float32_t] deltax,
                  np.ndarray[ndim=1, dtype=np.float32_t] deltay,
                  int device):
        """Create a sutra_atmos object

        :parameters:
            c: (naga_context): context

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

        cdef np.ndarray[ndim = 1, dtype = np.int64_t]size2
        size2 = compute_size2(size)

        self.s_a = new sutra_atmos(ctxt.c, nscreens,
                                   < np.float32_t * > r0.data,
                                   < long * > size.data,
                                   < long * > size2.data,
                                   < np.float32_t * > altitude.data,
                                   < np.float32_t * > windspeed.data,
                                   < np.float32_t * > winddir.data,
                                   < np.float32_t * > deltax.data,
                                   < np.float32_t * > deltay.data,
                                   device)
        self.context = ctxt

    def get_screen(self, float alt):
        """Return a numpy array containing the turbulence at a given altitude

        :param alt: (float) :altitude of the screen to get
        """
        cdef carma_obj[float] * screen = self.s_a.d_screens[alt].d_tscreen.d_screen
        self.context.set_activeDevice(screen.getDevice(), 1)
        cdef const long * dims
        dims = screen.getDims()
        cdef np.ndarray data_F = np.empty((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray data = np.empty((dims[1], dims[2]), dtype=np.float32)
        screen.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (dims[1], dims[2]))
        return data

    def disp(self, float alt):
        """Display the screen phase at a given altitude

        :param alt: (float) : altitude of the screen to display
        """
        cdef carma_obj[float] * c_phase = self.s_a.d_screens[alt].d_tscreen.d_screen
        cdef const long * dims = c_phase.getDims()
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data = np.zeros((dims[2], dims[1]), dtype=np.float32)
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] phase = np.ndarray((dims[1], dims[2]),
                                                                           dtype=np.float32)

        c_phase.device2host( < float * > data.data)
        phase = np.reshape(data.flatten("F"), (dims[1], dims[2]))

        pl.ion()
        pl.clf()
        pl.imshow(phase, cmap="Blues")
        pl.show()

    def add_screen(self, long size, float amplitude, float altitude,
                   float windspeed, float winddir, float deltax, float deltay, int device):
        """Add a screen to the atmos object.

        :parameters:
            size: (float) : dimension of the screen (size x size)

            amplitude: (float) : frac

            altitude: (float) : altitude of the screen in meters

            windspeed: (float) : windspeed of the screen [m/s]

            winddir: (float) : wind direction (deg)

            deltax: (float) : extrude deltax pixels in the x-direction at each iteration

            deltay: (float) : extrude deltay pixels in the y-direction at each iteration

            device: (int) : device number
        """
        cdef long size2 = compute_size2(np.array([size], dtype=np.int64))[0]

        if(self.s_a.d_screens.find(altitude) != self.s_a.d_screens.end()):
            print "There is already a screen at this altitude"
            print "No screen created"
            return

        cdef sutra_tscreen * screen = new sutra_tscreen(self.s_a.current_context, size, size2, amplitude, altitude, windspeed, winddir, deltax, deltay, device)

        cdef pair[float, sutra_tscreen * ] p
        p.first, p.second = altitude, screen
        self.s_a.d_screens.insert(p)
        self.s_a.nscreens += 1

    def del_screen(self, float alt):
        """Delete a screen from the atmos object

        :param alt: (float) : altitude of the screen to delete
        """
        if(self.s_a.d_screens.find(alt) == self.s_a.d_screens.end()):
            print "No screen at this altitude"
            print "No screen deleted"
            return
        self.s_a.nscreens -= 1
        self.s_a.d_screens.erase(alt)

    def list_alt(self):
        """Display the list of the screens altitude"""

        cdef map[float, sutra_tscreen * ].iterator it
        cdef int i = 0
        cdef np.ndarray alt = np.zeros(self.s_a.nscreens, dtype=np.float32)
        it = self.s_a.d_screens.begin()

        while it != self.s_a.d_screens.end():
            alt[i] = deref(it).first
            inc(it)
            i += 1
        print alt

    def move_atmos(self):
        """Move the turbulence in the atmos screen following previous loaded
        paramters such as windspeed and wind direction
        """

        self.s_a.move_atmos()

    def __str__(self):
        cdef map[float, sutra_tscreen * ].iterator it = self.s_a.d_screens.begin()
        cdef sutra_tscreen * screen
        cdef int i = 1
        info = "Atmos obect:\n"
        info += "Contains " + \
            str(self.s_a.nscreens) + " turbulent screen(s):\n"
        info += "Screen # | alt.(m) | speed (m/s) | dir.(deg) | r0 (pix) | deltax | deltay\n"
        while it != self.s_a.d_screens.end():
            screen = deref(it).second
            info += "%8d" % i + " | " + "%7.4f" % screen.altitude + " | " + "%11.4f" % screen.windspeed + \
                " | " + "%9.4f" % screen.winddir + " | " + "%8.4f" % screen.amplitude ** -(6. / 5.) + \
                " | " + "%6.4f" % screen.deltax + \
                " | " + "%6.4f" % screen.deltay + "\n"
            i = i + 1
            inc(it)

        info += "--------------------------------------------------------"

        return info


cdef atmos_create(naga_context c, int nscreens,
                  float r0,
                  np.ndarray[dtype=np.float32_t] L0,
                  float pupixsize,
                  np.ndarray[ndim=1, dtype=np.int64_t] dim_screens,
                  np.ndarray[ndim=1, dtype=np.float32_t] frac,
                  np.ndarray[ndim=1, dtype=np.float32_t] alt,
                  np.ndarray[ndim=1, dtype=np.float32_t] windspeed,
                  np.ndarray[ndim=1, dtype=np.float32_t] winddir,
                  np.ndarray[ndim=1, dtype=np.float32_t] deltax,
                  np.ndarray[ndim=1, dtype=np.float32_t] deltay,
                  np.ndarray[ndim=1, dtype=np.int64_t] seeds,
                  int verbose, int clean, dict load):
    """atmos_create(naga_context c, int nscreens,
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
                    np.ndarray[ndim=1,dtype=np.int64_t] seeds,
                    int verbose, int clean, dict load)

    Create and initialise an atmos object.

    :parameters:
        c: (naga_context) : context

        nscreens: (float) : number of turbulent layers

        r0: (float) : global r0

        L0: (np.ndarray[ndim=1, dtype=np.float32_t]) : L0

        pupixsize: (float) : pixel size [m]

        dim_screens: (np.ndarray[ndim=1,dtype=np.int64_t]) : screens dimensions

        frac: (np.ndarray[ndim=1,dtype=np.float32_t]) : fraction of r0

        alt: (np.ndarray[ndim=1,dtype=np.float32_t]) : altitudes [m]

        windspeed: (np.ndarray[ndim=1,dtype=np.float32_t]) : wind speed [m/s]

        winddir: (np.ndarray[ndim=1,dtype=np.float32_t]) : wind direction [deg]

        deltax: (np.ndarray[ndim=1,dtype=np.float32_t]) : extrude deltax pixels in the x-direction at each iteration

        deltay: (np.ndarray[ndim=1,dtype=np.float32_t]) : extrude deltay pixels in the y-direction at each iteration

        seeds: (np.ndarray[ndim=1,dtype=np.float32_t]) : seed for each screen

        verbose: (int) : 0 or 1
    """

    # get fraction of r0 for corresponding layer
    r0_layers = r0 / (frac ** (3. / 5.) * pupixsize)
    # create atmos object on gpu
    atmos_obj = Atmos()

    cdef carma_context * context = &carma_context.instance()
    cdef int device = context.get_activeDevice()
    atmos_obj.realinit(naga_context(), nscreens, r0_layers, dim_screens, alt,
                       windspeed, winddir, deltax, deltay, device)

    cdef int i, j
    cdef np.ndarray[ndim = 2, dtype = np.float32_t] A, B, A_F, B_F
    cdef np.ndarray[ndim = 1, dtype = np.uint32_t] istx, isty
    cdef sutra_tscreen * tscreen

    cdef str file_A, file_B, file_istx, file_isty

    for i in range(nscreens):
        if(load.has_key("A")):
            print "loading", load["A"]
            f = h5py.File(load["A"])
            A = f["A"][:]
            f.close()
            print "loading", load["B"]
            f = h5py.File(load["B"])
            B = f["B"][:]
            f.close()
            print "loading", load["istx"]
            f = h5py.File(load["istx"])
            istx = f["istx"][:]
            f.close()
            print "loading", load["isty"]
            f = h5py.File(load["isty"])
            isty = f["isty"][:]
            f.close()

        else:
            A, B, istx, isty = itK.AB(dim_screens[i], L0[i], deltax[i], deltay[i], verbose)
            if not(clean):
                version = check_output(["git", "rev-parse", "HEAD"]).replace("\n", "")
                print "writing files and updating database"
                df = pandas.read_hdf(
                    shesha_savepath + "/matricesDataBase.h5", "A")
                ind = len(df.index) - 1
                savename = shesha_savepath + "/turbu/A_r" + \
                    version + "_" + str(ind) + ".h5"
                h5u.save_hdf5(savename, "A", A)

                savename = shesha_savepath + "/turbu/B_r" + \
                    version + "_" + str(ind) + ".h5"
                h5u.save_hdf5(savename, "B", B)

                savename = shesha_savepath + "/turbu/istx_r" + \
                    version + "_" + str(ind) + ".h5"
                h5u.save_hdf5(savename, "istx", istx)

                savename = shesha_savepath + "/turbu/isty_r" + \
                    version + "_" + str(ind) + ".h5"
                h5u.save_hdf5(savename, "isty", isty)

        A_F = np.reshape(A.flatten("F"), (A.shape[0], A.shape[1]))
        B_F = np.reshape(B.flatten("F"), (B.shape[0], B.shape[1]))

        tscreen = atmos_obj.s_a.d_screens[alt[i]]
        tscreen.init_screen(< float * > (A_F.data), < float * > (B_F.data),
                             < unsigned int * > istx.data, < unsigned int * > isty.data, seeds[i])
        for j in range(2 * tscreen.screen_size):
            tscreen.extrude(1 * np.sign(deltax[i]))
    return atmos_obj


cdef compute_size2(np.ndarray[ndim=1, dtype=np.int64_t] size):
    """compute_size2(np.ndarray[ndim=1, dtype=np.int64_t] size)

    Compute the size of a stencil, given the screen size

    :parameters:
        size: (np.ndarray[ndim=1,dtype=np.int64_t]) :screen size
    """

    cdef n = size.shape[0]
    cdef np.ndarray[ndim = 1, dtype = np.int64_t] size2 = np.zeros(n, dtype=np.int64)

    cdef int i
    for i in range(n):
        size2[i] = itK.stencil_size(size[i])
    return size2
