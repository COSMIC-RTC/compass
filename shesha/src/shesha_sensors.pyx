include "../par.pxi"

import numpy as np
cimport numpy as np
np.import_array()

import make_pupil as mkP

import os

def rebin(a, shape):
    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    return a.reshape(sh).mean(-1).mean(1)

cpdef prep_lgs_prof(Param_wfs p_wfs, int nsensors, Param_tel p_tel,
                    np.ndarray[dtype=np.float32_t] prof,
                    np.ndarray[dtype=np.float32_t] h,
                    float beam, Sensors sensors,
                    bytes center= < bytes > "", int imat=0):
    """The function returns an image array(double,n,n) of a laser beacon elongated by perpective
effect. It is obtaind by convolution of a gaussian of width "lgsWidth" arcseconds, with the
line of the sodium profile "prof". The altitude of the profile is the array "h".



    :parameters:
        nsensors: (int) : wfs index

        p_tel: (Param_tel) : telescope settings

        prof: (np.ndarray[dtype=np.float32]) : Na profile intensity, in arbitrary units

        h: (np.ndarray[dtype=np.float32]) : altitude, in meters. h MUST be an array with EQUALLY spaced elements.

        beam: (float) : size in arcsec of the laser beam

        center: (string) : either "image" or "fourier" depending on where the centre should be.

Computation of LGS spot from the sodium profile:
Everything is done here in 1D, because the Na profile is the result of the convolution of a function
P(x,y) = profile(x) . dirac(y)
by a gaussian function, for which variables x and y can be split :
exp(-(x^2+y^2)/2.s^2)  =  exp(-x^2/2.s^2) * exp(-y^2/2.s^2)
The convolution is (symbol $ denotes integral)
C(X,Y) = $$ exp(-x^2/2.s^2) * exp(-y^2/2.s^2) * profile(x-X) * dirac(y-Y)  dx  dy
First one performs the integration along y
C(X,Y) = exp(-Y^2/2.s^2)  $ exp(-x^2/2.s^2) * profile(x-X)  dx
which shows that the profile can be computed by
- convolving the 1-D profile
- multiplying it in the 2nd dimension by a gaussian function

If one has to undersample the inital profile, then some structures may be "lost". In this case,
it's better to try to "save" those structures by re-sampling the integral of the profile, and
then derivating it afterwards.
Now, if the initial profile is a coarse one, and that one has to oversample it, then a
simple re-sampling of the profile is adequate.
"""

    p_wfs._prof1d = prof
    p_wfs._profcum = np.zeros(prof.size + 1, dtype=np.float32)
    p_wfs._profcum[1:] = prof.cumsum()
    cdef float subapdiam = p_tel.diam / p_wfs.nxsub  # diam of subap
    cdef np.ndarray[dtype = np.float32_t] xsubs, ysubs, dOffAxis, g
    if(p_wfs.nxsub > 1):
        xsubs = np.linspace((subapdiam - p_tel.diam) / 2, (p_tel.diam - subapdiam) / 2,
                            p_wfs.nxsub).astype(np.float32)
    else:
        xsubs = np.zeros(1, dtype=np.float32)
    ysubs = xsubs.copy().astype(np.float32)

    # cdef int nP=prof.shape[0]  #number of points of the profile #UNUSED
    # center of gravity of the profile
    cdef float hG = np.sum(h * prof) / np.sum(prof)
    cdef np.ndarray[dtype = np.float32_t] x = np.arange(p_wfs._Ntot).astype(np.float32)\
        - p_wfs._Ntot / 2
    # x expressed in pixels. (0,0) is in the fourier-center

    x = x * p_wfs._qpixsize  # x expressed in arcseconds
    # cdef float dx=x[1]-x[0] #UNUSED
    cdef float dh = h[1] - h[0]

    if(p_wfs.nxsub > 1):
        dOffAxis = np.sqrt((xsubs[p_wfs._validsubsx] - p_wfs.lltx) ** 2 +
                           (ysubs[p_wfs._validsubsy] - p_wfs.llty) ** 2)
    else:
        dOffAxis = np.sqrt(
            (xsubs - p_wfs.lltx) ** 2 + (ysubs - p_wfs.llty) ** 2)

    if(imat > 0):
        dOffAxis *= 0.

    cdef w = beam / 2.35482005
    if(w == 0):
        # TODO what is n
        n = 1
        g = np.zeros(n, dtype=np.float32)
        if(center == "image"):
            g[n / 2 - 1] = 0.5
            g[n / 2] = 0.5
        else:
            g[n / 2] = 1

    else:
        if(center == "image"):
            g = np.exp(-(x + p_wfs._qpixsize / 2) ** 2 / (2 * w ** 2.))
        else:
            g = np.exp(-x ** 2 / (2 * w ** 2.))

    p_wfs._ftbeam = np.fft.fft(g, axis=0).astype(np.complex64)
    p_wfs._beam = g
    # convolved profile in 1D.

    if(xsubs.size > 1):
        azimuth = np.arctan2(ysubs[p_wfs._validsubsy] - p_wfs.llty,
                             xsubs[p_wfs._validsubsx] - p_wfs.lltx)
    else:
        azimuth = np.arctan2(ysubs - p_wfs.llty,
                             xsubs - p_wfs.lltx)

    p_wfs._azimuth = azimuth

    cdef carma_context * context = carma_context.instance()
    context.set_activeDevice(sensors.sensors.device, 1)

    cdef sutra_lgs * lgs = sensors.sensors.d_wfs[nsensors].d_gs.d_lgs

    lgs.lgs_init(p_wfs._prof1d.size, hG, h[0], dh, p_wfs._qpixsize,
                 < float * > dOffAxis.data, < float * > p_wfs._prof1d.data,
                 < float * > p_wfs._profcum.data, < float * > p_wfs._beam.data,
                 < cuFloatComplex * > p_wfs._ftbeam.data, < float * > p_wfs._azimuth.data)

    lgs.lgs_update(context.get_device(sensors.sensors.device))
    lgs.lgs_makespot(context.get_device(sensors.sensors.device), 0)


cpdef make_lgs_prof1d(p_wfs, Param_tel p_tel,
                      np.ndarray[dtype=np.float32_t] prof, np.ndarray[dtype=np.float32_t] h,
                      float beam, bytes center= < bytes > ""):
    """same as prep_lgs_prof but cpu only. original routine from rico

    :parameters:
        p_tel: (Param_tel) : telescope settings

        prof: (np.ndarray[dtype=np.float32]) : Na profile intensity, in arbitrary units

        h: (np.ndarray[dtype=np.float32]) : altitude, in meters. h MUST be an array with EQUALLY spaced elements.

        beam: (float) : size in arcsec of the laser beam

        center: (string) : either "image" or "fourier" depending on where the centre should be.
    """

    p_wfs._prof1d = prof
    p_wfs._profcum = np.zeros(prof.shape[1] + 1, dtype=np.float32)
    p_wfs._profcum[1:] = prof.cumsum()

    cdef float subapdiam = p_tel.diam / p_wfs.nxsub  # diam of subap
    cdef np.ndarray[dtype = np.float32_t] xsubs, ysubs, dOffAxis
    cdef np.ndarray[dtype = np.float32_t] g
    if(p_wfs.nxsub > 1):
        xsubs = np.linspace((subapdiam - p_tel.diam) / 2, (p_tel.diam - subapdiam) / 2,
                            p_wfs.nxsub).astype(np.float32)
    else:
        xsubs = np.zeros(1, dtype=np.float32)
    ysubs = xsubs.copy().astype(np.float32)

    # cdef int nP=prof.shape[0] #UNUSED
    cdef float hG = np.sum(h * prof) / np.sum(prof)
    cdef np.ndarray[dtype = np.float32_t] x = np.arange(p_wfs._Ntot).astype(np.float32) - p_wfs._Ntot / 2
    # x expressed in pixels. (0,0) is in the fourier-center.
    x = x * p_wfs._qpixsize  # x expressed in arcseconds
    # cdef float dx=x[1]-x[0] #UNUSED
    # cdef float dh=h[1]-h[0] #UNUSED

    if(p_wfs.nxsub > 1):
        dOffAxis = np.sqrt((xsubs[p_wfs._validsubsy] - p_wfs.lltx) ** 2 +
                           (ysubs[p_wfs._validsubsx] - p_wfs.llty) ** 2)
    else:
        dOffAxis = np.sqrt(
            (xsubs - p_wfs.lltx) ** 2 + (ysubs - p_wfs.llty) ** 2)

    cdef np.ndarray[ndim = 2, dtype = np.float32_t] profi,
    cdef np.ndarray[ndim = 1, dtype = np.float32_t] zhc, avg_zhc, avg_x
    profi = np.zeros((p_wfs._Ntot, p_wfs._nvalid), dtype=np.float32)

    # cdef np.ndarray[ndim=1,dtype=np.int32_t] subsdone, dif2do #UNUSED
    cdef np.ndarray[ndim = 1, dtype = np.int64_t] inds
    subsdone = np.ones(p_wfs._nvalid, dtype=np.int32)
    # dif2do =np.zeros(p_wfs._nvalid,dtype=np.int32) #UNUSED
    cdef float tmp, dzhc

    cdef int i

    # BUG put a correct value to tmp
    raise "tmp not initialized... in yorick, it was: ...\
           tmp = array(float,[3,y_wfs(numwfs)._Ntot,y_wfs(numwfs)._Ntot,y_wfs(numwfs)._nvalid])"
    zhc = (h.astype(np.float32) - hG) * (206265. * tmp / hG ** 2)
    dzhc = zhc[1] - zhc[0]

    tmp = dOffAxis[np.where(subsdone)][0]
    zhc = (h - hG) * (206265. * tmp / hG ** 2)
    avg_zhc = np.zeros(zhc.size + 1, dtype=np.float32)
    avg_zhc[0] = zhc[0]
    avg_zhc[avg_zhc.size - 1] = zhc[zhc.size - 1]
    avg_zhc[1:-1] = 0.5 * (zhc[1:] + zhc[:-1])
    avg_x = np.zeros(x.size + 1, dtype=np.float32)
    avg_x[0] = x[0]
    avg_x[avg_x.size - 1] = x[x.size - 1]
    avg_x[1:-1] = 0.5 * (x[1:] + x[:-1])

    while(np.any(subsdone)):
        tmp = dOffAxis[np.where(subsdone)][0]
        inds = np.where(dOffAxis == tmp)[0]
        # height, translated in arcsec due to perspective effect
        zhc = (h - hG) * (206265. * tmp / hG ** 2)
        dzhc = zhc[1] - zhc[0]

        if(p_wfs._qpixsize > dzhc):
            avg_zhc = np.zeros(zhc.size + 1, dtype=np.float32)
            avg_zhc[0] = zhc[0]
            avg_zhc[avg_zhc.size - 1] = zhc[zhc.size - 1]
            avg_zhc[1:-1] = 0.5 * (zhc[1:] + zhc[:-1])
            avg_x = np.zeros(x.size + 1, dtype=np.float32)
            avg_x[0] = x[0]
            avg_x[avg_x.size - 1] = x[x.size - 1]
            avg_x[1:-1] = 0.5 * (x[1:] + x[:-1])

            for i in range(inds.size):
                profi[:, inds[i]] = np.diff(
                    np.interp(avg_x, avg_zhc, p_wfs._profcum)).astype(np.float32)

        else:
            for i in range(inds.size):
                profi[:, inds[i]] = np.interp(x, zhc, prof)
        subsdone[inds] = 0

    cdef w = beam / 2.35482005
    if(w == 0):
        # TODO what is n
        n = 1
        g = np.zeros(n, dtype=np.float32)
        if(center == "image"):
            g[n / 2 - 1] = 0.5
            g[n / 2] = 0.5
        else:
            g[n / 2] = 1

    else:
        if(center == "image"):
            if((p_wfs.npix * p_wfs._nrebin) % 2 != p_wfs._Nfft % 2):
                g = np.exp(-(x + p_wfs._qpixsize) ** 2 / (2 * w ** 2.))
            else:
                g = np.exp(-(x + p_wfs._qpixsize / 2) ** 2 / (2 * w ** 2.))

        else:
            g = np.exp(-x ** 2 / (2 * w ** 2.))

    p_wfs._ftbeam = np.fft.fft(g).astype(np.complex64)
    p_wfs._beam = g.astype(np.float32)
    # convolved profile in 1D.

    cdef np.ndarray[ndim = 2, dtype = np.float32_t] p1d
    cdef np.ndarray[ndim = 2, dtype = np.float32_t] g_extended = np.tile(g, (p_wfs._nvalid, 1)).T

    p1d = np.fft.ifft(np.fft.fft(profi, axis=0) *
                      np.fft.fft(g_extended, axis=0),
                      axis=0).real.astype(np.float32)
    p1d = p1d * p1d.shape[0]
    p1d = np.roll(p1d, int(p_wfs._Ntot / 2. + 0.5), axis=0)
    p1d = np.abs(p1d)

    cdef np.ndarray[ndim = 3, dtype = np.float32_t] im
    im = np.zeros((p1d.shape[1], p1d.shape[0], p1d.shape[0]), dtype=np.float32)

    cdef int l, c
    for i in range(p1d.shape[1]):
        for l in range(p1d.shape[0]):
            for c in range(p1d.shape[0]):
                im[i, l, c] = g[l] * p1d[c, i]

    if(ysubs.size > 1):
        azimuth = np.arctan2(ysubs[p_wfs._validsubsy] - p_wfs.llty,
                             xsubs[p_wfs._validsubsx] - p_wfs.lltx)
    else:
        azimuth = np.arctan2(ysubs - p_wfs.llty,
                             xsubs - p_wfs.lltx)

    p_wfs._azimuth = azimuth

    cdef float xcent, ycent

    if(center == "image"):
        xcent = p_wfs._Ntot / 2. + 0.5
        ycent = xcent
    else:
        xcent = p_wfs._Ntot / 2. + 1
        ycent = xcent

    cdef np.ndarray[ndim = 1, dtype = np.float32_t] max_im
    if(ysubs.size > 0):
        # TODO rotate
        im = rotate3d(im, azimuth * 180 / np.pi, xcent, ycent)
        max_im = np.max(im, axis=(1, 2))
        im = (im.T / max_im).T
    else:
        im = rotate(im, azimuth * 180 / np.pi, xcent, ycent)
        im = im / np.max(im)

    p_wfs._lgskern = im.T


cpdef type_present(liste, int pyr, int roof, int sh, int geo):
    """Check the present types in a list

    :parameters:
        liste: (list of str) : list of types

        pyr: (int) : set to 1 if the list contains "pyr" (0 else)

        roof: (int): set to 1 if the list contains "roof" (0 else)

        sh: (int) : set to 1 if the list contains "sh" (0 else)

        geo: (int) : set to 1 if the list contains "geo" (0 else)

    return 1 if the wfs type is present (0 else)
    """

    cdef int i, l

    cdef bytes wfs_type

    l = len(liste)
    for i in range(l):
        wfs_type = liste[i]
        if(wfs_type == "pyr"):
            pyr = 1
        elif (wfs_type == "roof"):
            roof = 1
        elif(wfs_type == "sh"):
            sh = 1
        elif(wfs_type == "geo"):
            geo = 1


cpdef wheremax(liste):
    """return the index of the maximum value of the list

    :param liste: (list of values) : values to get the index of the maximum from
    """
    cdef int i, j
    cdef int l = len(liste)
    cdef float m = 0
    cdef float val
    j = -1
    for i in range(l):
        val = liste[i]
        if(m < val):
            j = i
            m = val
    return j


def wfs_init(wfs, Param_atmos p_atmos, Param_tel p_tel, Param_geom p_geom,
             Param_target p_target, Param_loop p_loop, dm=None, int comm_size=1, int rank=0):
    """
    Create and initialise  a Sensors object

    :parameters:
        wfs: (list of Param_wfs) : wfs settings

        p_atmos: (Param_atmos) : atmos settings

        p_tel: (Param_tel) : telescope settings

        p_geom: (Param_geom) : geom settings

        p_target: (Param_target) : target settings

        p_loop: (Param_loop) : loop settings

        dm: (list of Param_dm) : (optional) dms settings [=None]

        comm_size: (int) : (optional) communicator size [=1]

        rank: (int) : (optional) process rank [=0]
    """
    cdef int nsensors = len(wfs)
    cdef int i
    cdef Param_wfs o
    cdef int any_pyr = 0
    cdef int any_roof = 0
    cdef int any_sh = 0
    cdef int any_geo = 0

    liste = [o.type_wfs for o in wfs]
    type_present(liste, any_pyr, any_roof, any_sh, any_geo)

    # dm = None
    if(wfs[0].dms_seen is None and dm is not None):
        for i in range(nsensors):
            if(not wfs[i].openloop):
                wfs[i].set_dms_seen(np.arange(len(dm), dtype=np.int32))

    cdef int indmax
    # first get the wfs with max # of subaps
    # we'll derive the geometry from the requirements in terms of sampling
    if(any_sh):
        indmax = wheremax([o.nxsub for o in wfs if o.type_wfs == "sh"])
    else:
        indmax = wheremax([o.nxsub for o in wfs])

    # init geometry
    if(any_geo == 0):
        init_wfs_geom(wfs[indmax], wfs[0], indmax, p_atmos, p_tel, p_geom, p_target,
                      p_loop, init=1, verbose=rank)
    else:
        init_wfs_geom(wfs[indmax], wfs[0], indmax, p_atmos, p_tel, p_geom, p_target,
                      p_loop, init=0, verbose=rank)

    # #do the same for other wfs
    for i in range(nsensors):
        if(i != indmax):
            init_wfs_geom(wfs[i], wfs[i], i, p_atmos, p_tel, p_geom, p_target,
                          p_loop, init=0, verbose=rank)
    # create sensor object on gpu
    # and init sensor gs object on gpu

    # arrays needed to call Sensors constructor
    t_wfs = [o.type_wfs for o in wfs]
    # cdef np.ndarray t_wfs  = np.array([o.type_wfs  for o in
    # wfs],dtype=np.str)
    cdef np.ndarray nxsub = np.array([o.nxsub for o in wfs], dtype=np.int64)
    cdef np.ndarray nvalid = np.array([o._nvalid for o in wfs], dtype=np.int64)
    cdef np.ndarray nphase = np.array([o._pdiam for o in wfs], dtype=np.int64)
    cdef np.ndarray pdiam = np.array([o._subapd for o in wfs], dtype=np.float32)
    cdef np.ndarray npix = np.array([o.npix for o in wfs], dtype=np.int64)
    cdef np.ndarray nrebin = np.array([o._nrebin for o in wfs], dtype=np.int64)
    cdef np.ndarray nfft = np.array([o._Nfft for o in wfs], dtype=np.int64)
    cdef np.ndarray ntota = np.array([o._Ntot for o in wfs], dtype=np.int64)
    cdef np.ndarray nphot = np.array([o._nphotons for o in wfs], dtype=np.float32)
    cdef np.ndarray nphot4imat = np.array([o.nphotons4imat for o in wfs], dtype=np.float32)
    cdef np.ndarray lgs = np.array([o.gsalt > 0 for o in wfs], dtype=np.int32)

    # arrays needed to call sensors_initgs
    cdef np.ndarray xpos = np.array([o.xpos for o in wfs], dtype=np.float32)
    cdef np.ndarray ypos = np.array([o.ypos for o in wfs], dtype=np.float32)
    cdef np.ndarray Lambda = np.array([o.Lambda for o in wfs], dtype=np.float32)
    cdef np.ndarray mag
    cdef float zerop = wfs[0].zerop
    cdef np.ndarray size = np.zeros(nsensors, dtype=np.int64) + p_geom._n
    cdef np.ndarray noise
    cdef np.ndarray seed = np.array([], dtype=np.int64)
    cdef np.ndarray npup = (np.zeros((nsensors)) + p_geom._n).astype(np.int64)

    cdef np.ndarray tmp

    error_budget_flag = [w.error_budget for w in wfs]
    if(True in error_budget_flag):
        error_budget_flag = True
    else:
        error_budget_flag = False

    telescope = Telescope(p_geom._spupil.shape[0], np.where(p_geom._spupil > 0)[0].size,
                          p_geom._spupil *
                          p_geom._apodizer, p_geom._phase_ab_M1,
                          p_geom._mpupil.shape[0], p_geom._mpupil, p_geom._phase_ab_M1_m)

    if(wfs[0].type_wfs == "sh"):
        g_wfs = Sensors(nsensors, telescope, t_wfs, npup, nxsub, nvalid, nphase, pdiam, npix, nrebin,
                        nfft, ntota, nphot, nphot4imat, lgs, comm_size=comm_size, rank=rank, error_budget=error_budget_flag)

        mag = np.array([o.gsmag for o in wfs], dtype=np.float32)
        noise = np.array([o.noise for o in wfs], dtype=np.float32)
        g_wfs.sensors_initgs(xpos, ypos, Lambda, mag, zerop, size, noise, seed)

    elif(wfs[0].type_wfs=="pyr" or wfs[0].type_wfs == "roof"):
        npup = np.array([wfs[0].pyr_npts])
        g_wfs = Sensors(nsensors, telescope, t_wfs, npup, nxsub, nvalid, nphase, pdiam, npix, nrebin,
                        nfft, ntota, nphot, nphot4imat, lgs, comm_size=comm_size, rank=rank, error_budget=error_budget_flag)

        mag = np.array([o.gsmag for o in wfs], dtype=np.float32)
        noise = np.array([o.noise for o in wfs], dtype=np.float32)
        g_wfs.sensors_initgs(xpos, ypos, Lambda, mag, zerop, size, noise, seed)

    elif(wfs[0].type_wfs == "pyrhr"):
        npup = np.array([wfs[0].pyr_npts])
        g_wfs = Sensors(nsensors, telescope, t_wfs, npup, nxsub, nvalid, nphase, pdiam, npix, nrebin,
                        nfft, ntota, nphot, nphot4imat, lgs, comm_size=comm_size, rank=rank, error_budget=error_budget_flag)

        mag = np.array([o.gsmag for o in wfs], dtype=np.float32)
        noise = np.array([o.noise for o in wfs], dtype=np.float32)
        g_wfs.sensors_initgs(xpos, ypos, Lambda, mag, zerop, size, noise, seed)

    elif(wfs[0].type_wfs == "geo"):
        npup = np.array([wfs[0].p_geom._n])
        g_wfs = Sensors(nsensors, telescope, wfs[0].type_wfs, npup, nxsub, nvalid, nphase, pdiam,
                        comm_size=comm_size, rank=rank)

        mag = np.zeros(nsensors - 1, dtype=np.float32)
        noise = np.zeros(nsensors - 1, dtype=np.float32) - 1
        g_wfs.sensors_initgs(xpos, ypos, Lambda, mag, zerop, size, noise, seed)

    # fill sensor object with data
    for i in range(nsensors):
        g_wfs.sensors_initarr(i, wfs[i])

    # lgs case
    for i in range(nsensors):
        if(wfs[i].gsalt > 0):
            # lgs mode requested
            if(wfs[i].proftype is None or wfs[i].proftype == ""):
                wfs[i].set_proftype("Gauss1")

            if(wfs[i].proftype == "Gauss1"):
                profilename = "allProfileNa_withAltitude_1Gaussian.npy"
            elif(wfs[i].proftype == "Gauss2"):
                profilename = "allProfileNa_withAltitude_2Gaussian.npy"
            elif(wfs[i].proftype == "Gauss3"):
                profilename = "allProfileNa_withAltitude_3Gaussian.npy"
            elif(wfs[i].proftype == "Exp"):
                profilename = "allProfileNa_withAltitude.npy"
            else:
                error = "Param_wfs proftype unknown: got '" + \
                    wfs[i].proftype + \
                    "', expect one of: \n''\n'Gauss1'\n'Gauss2'\n'Gauss3'\n'Exp'"
                raise ValueError(error)
            print "reading Na profile from", shesha_savepath + profilename
            prof = np.load(shesha_savepath + profilename)
            wfs[i].set_altna(prof[0, :].astype(np.float32))
            wfs[i].set_profna(np.mean(prof[1:, :], axis=0).astype(np.float32))
            # init sensor gs object with necessary data
            prep_lgs_prof(wfs[i], i, p_tel, wfs[i]._profna, wfs[i]._altna,
                          wfs[i].beamsize, g_wfs)

    return g_wfs, telescope

cpdef init_wfs_geom(Param_wfs wfs, Param_wfs wfs0, int n, Param_atmos atmos,
                    Param_tel tel, Param_geom geom, Param_target p_target,
                    Param_loop loop, int init=0, int verbose=0):
    """Compute the geometry of WFSs: valid subaps, positions of the subaps,
    flux per subap, etc...

    :parameters:
        wfs: (Param_wfs) : wfs settings

        wfs0: (Param_wfs) : reference wfs settings

        n: (int) : index of the wfs (diplay information purpose only)

        atmos: (Param_atmos) : atmos settings

        tel: (Param_tel) : telescope settings

        geom: (Param_geom) : geom settings

        target: (Param_target) : target settings

        loop: (Param_loop) : loop settings

        init: (int) : (optional)

        verbose: (int) : (optional) display informations if 0


    """

    if(verbose == 0):
        print "*-----------------------"
    if(verbose == 0):
        print "Doing inits on WFS", n

    cdef long pdiam = 0

    if(init == 0):
        if(wfs.type_wfs == "sh"):
            pdiam = geom.pupdiam / wfs.nxsub
            if(geom.pupdiam % wfs.nxsub > 0):
                pdiam += 1
        if(wfs.type_wfs == "geo"):
            if(geom.pupdiam == 0):
                pdiam = 20
            else:
                pdiam = geom.pupdiam / wfs.nxsub
                if(geom.pupdiam % wfs.nxsub > 0):
                    pdiam += 1
    else:
        pdiam = -1

    cdef int Nfft = 0
    cdef int Ntot = 0
    cdef int nrebin = 0
    cdef float pixsize = 0
    cdef float qpixsize = 0
    # , validsubs
    cdef np.ndarray pyr_focmask, pupvalid, pupreb, istart, jstart
    cdef np.ndarray phase_shift, pshift, cx, cy, phasemap, tmp, sincar, fluxPerSub, halfxy
    cdef np.ndarray validx, validy, isvalid
    cdef np.ndarray binmap, binindices

    cdef np.ndarray x = np.empty((1, 1), dtype=np.float32)
    cdef np.ndarray y = np.empty((1, 1), dtype=np.float32)

    cdef int i, j, indi  # , indj
    cdef long n1, n2, npup
    cdef float coef1, coef2

    # TODO define psize properly
    # defined in yoga_rtc.i /and yoga_dm.i as
    #   (1): y_geom.pupdiam
    #   (2): y_tel.diam/y_geom.pupdiam
    cdef int psize = 0  # int(tel.diam/geom.pupdiam)
    init_wfs_size(wfs, n, atmos, tel, psize, & pdiam, & Nfft, & Ntot, & nrebin, & pixsize, & qpixsize, verbose)

    if(wfs.type_wfs != "geo"):
        wfs.pixsize = pixsize
        wfs._Nfft = Nfft
        wfs._Ntot = Ntot
        wfs._nrebin = nrebin
        wfs._qpixsize = qpixsize

    wfs._subapd = tel.diam / wfs.nxsub

    wfs._pdiam = pdiam

    if(wfs.type_wfs == "pyr" or wfs.type_wfs == "roof"):
        wfs.npix = pdiam

    if(init == 1 or (wfs.type_wfs == "geo" and n == 1)):
        # this is the wfs with largest # of subaps
        # the overall geometry is deduced from it
        geom.geom_init(tel, pdiam * wfs.nxsub, p_target.apod)

    if(wfs.type_wfs=="pyr" or wfs.type_wfs == "roof"):
        padding = 2
        npup = wfs._Ntot
        n1 = geom.ssize / 2 - geom.pupdiam / 2 - padding * wfs.npix
        n2 = n1 + geom.pupdiam + 2 * padding * wfs.npix

        geom._mpupil = geom._ipupil[n1:n2, n1:n2]
        geom._n1 = n1
        geom._n2 = n2
        geom._n = npup
        geom._phase_ab_M1_m = mkP.pad_array(
            geom._phase_ab_M1, geom._n).astype(np.float32)

        # pup   = pup(ii,ii);
        # phase = phase(ii,ii);
        mod_ampl_pixels = wfs.pyr_ampl / wfs._qpixsize  # modulation in pixels
        fsradius_pixels = long(wfs.fssize / wfs._qpixsize / 2.)

        if (wfs.fstop == "round"):
            focmask = mkP.dist(
                npup, xc=npup / 2. + 0.5, yc=npup / 2. + 0.5) < (fsradius_pixels)
            # fstop_area = np.pi * (wfs.fssize/2.)**2. #UNUSED
        elif (wfs.fstop == "square"):
            x, y = indices(npup)
            x -= (npup + 1.) / 2.
            y -= (npup + 1.) / 2.
            focmask = (np.abs(x) <= (fsradius_pixels)) * \
                (np.abs(y) <= (fsradius_pixels))
            # fstop_area = wfs.fssize**2. #UNUSED
        else:
            msg = "wfs " + str(n) + ". fstop must be round or square"
            raise ValueError(msg)

        pyr_focmask = np.roll(focmask, focmask.shape[0] / 2, axis=0)
        pyr_focmask = np.roll(pyr_focmask, focmask.shape[1] / 2, axis=1)
        wfs._submask = pyr_focmask

        pup = geom._mpupil

        pupreb = bin2d(pup * 1., wfs.npix) / wfs.npix ** 2.
        wsubok = np.where(pupreb >= wfs.fracsub)
        pupvalid = pupreb * 0.
        pupvalid[wsubok] = 1
        wfs._nvalid = wsubok[0].size
        validx = np.where(pupvalid)[1].astype(np.int32)
        validy = np.where(pupvalid)[0].astype(np.int32)
        wfs._validsubsx = validx
        wfs._validsubsy = validy

        istart = (np.linspace(0.5, geom.pupdiam + 2 * padding *
                              wfs.npix + 0.5, wfs.nxsub + 2 * padding) + 1).astype(np.int32)[:-1]
        jstart = np.copy(istart)
        wfs._istart = istart.astype(np.int32)
        wfs._jstart = jstart.astype(np.int32)

        x, y = indices(npup)
        x -= (npup + 1) / 2.
        y -= (npup + 1) / 2.
        phase_shift = np.roll(
            np.exp(1j * 2 * np.pi * (0.5 * (x + y)) / npup), x.shape[0] / 2, axis=0)
        phase_shift = np.roll(phase_shift, x.shape[1] / 2, axis=1)
        wfs._halfxy = phase_shift

        x, y = indices(Nfft)
        x -= (Nfft + 1) / 2.
        y -= (Nfft + 1) / 2.
        if(wfs.nxsub * nrebin % 2 == 1):
            coef1 = 0.
        else:
            coef1 = -0.5
        if(wfs.nxsub % 2 == 1 and wfs.npix * nrebin % 2 == 1):
            coef2 = 1
        else:
            coef2 = 0.5

        pshift = np.exp(
            1j * 2 * np.pi * (coef1 / Nfft + coef2 * nrebin / wfs.npix / Nfft) * (x + y))
        wfs._pyr_offsets = pshift

        if(wfs.pyrtype == "Pyramid"):
            if(wfs.pyr_pos == None):
                cx = np.round(
                    mod_ampl_pixels * np.sin((np.arange(wfs.pyr_npts) + 1) * 2. * np.pi / wfs.pyr_npts))
                cy = np.round(
                    mod_ampl_pixels * np.cos((np.arange(wfs.pyr_npts) + 1) * 2. * np.pi / wfs.pyr_npts))
                # mod_npts = wfs.pyr_npts #UNUSED
            else:
                if(verbose == 0):
                    print "Using user-defined positions for the pyramid modulation"
                cx = np.round(wfs.pyr_pos[:, 0] / qpixsize)
                cy = np.round(wfs.pyr_pos[:, 1] / qpixsize)
                # mod_npts=cx.shape[0] #UNUSED
        elif(wfs.pyrtype == "RoofPrism"):
            cx = np.round(
                2. * mod_ampl_pixels * ((np.arange(wfs.pyr_npts) + 1) - (wfs.pyr_npts + 1) / 2.) / wfs.pyr_npts)
            cy = cx
            # mod_npts = wfs.pyr_npts #UNUSED
        else:
            if(wfs.pyr_pos==None):
                cx = np.round(mod_ampl_pixels*np.sin((np.arange(wfs.pyr_npts)+1)*2.*np.pi/wfs.pyr_npts))
                cy = np.round(mod_ampl_pixels*np.cos((np.arange(wfs.pyr_npts)+1)*2.*np.pi/wfs.pyr_npts))
                #mod_npts = wfs.pyr_npts #UNUSED
            else:
                if(verbose==0):print "Using user-defined positions for the pyramid modulation"
                cx=np.round(wfs.pyr_pos[:,0]/qpixsize)
                cy=np.round(wfs.pyr_pos[:,1]/qpixsize)
                #mod_npts=cx.shape[0] #UNUSED

        wfs._pyr_cx = cx.astype(np.int32)
        wfs._pyr_cy = cy.astype(np.int32)

        wfs._nphotons = wfs.zerop * \
            2.51189 ** (-wfs.gsmag) * loop.ittime * wfs.optthroughput

        # spatial filtering by the pixel extent:
        # *2/2 intended. min should be 0.40 = sinc(0.5)^2.
        x = x / (Nfft - 1) * 2 / 2
        y = y / (Nfft - 1) * 2 / 2
        sincar = np.roll(np.sinc(x) * np.sinc(y), x.shape[0] / 2, axis=0)
        sincar = np.roll(sincar, y.shape[0] / 2, axis=1)

        # sincar = np.roll(np.pi*x*np.pi*y,x.shape[1],axis=1)
        wfs._sincar = sincar.astype(np.float32)
        # must be set even if unused
        wfs._hrmap = np.array([0], dtype=np.int32)

        # this defines how we cut the phase into subaps
        phasemap = np.zeros((pdiam, pdiam, wfs._nvalid), dtype=np.int32)
        x, y = indices(geom._n)  # we need c-like indice
        x -= 1
        y -= 1
        tmp = x + y * geom._n
        for i in range(wfs._nvalid):
            indi = istart[wfs._validsubsx[i]] + 1  # +2-1 (yorick->python
            indj = jstart[wfs._validsubsy[i]] + 1
            phasemap[:, :, i] = tmp[indi:indi + pdiam, indj:indj + pdiam]

        wfs._phasemap = phasemap

    if(wfs.type_wfs == "pyrhr"):
        # nrebin[0]  = pdiam[0] / wfs.nxsub
        # Ntot[0]    = Nfft[0] /  pdiam[0] * wfs.nxsub
        # pixsize[0] = qpixsize[0] * nrebin[0]
        # pdiam[0]   = pdiam[0] / wfs.nxsub
        wfs.npix = pdiam

        # Initialize the intermediate pupil support size
        padding = nrebin
        n1 = geom.ssize / 2 - geom.pupdiam / 2 - padding
        n2 = n1 + geom.pupdiam + 2 * padding
        npup = geom.pupdiam + 2 * padding

        geom._mpupil = geom._ipupil[n1:n2, n1:n2]
        geom._n1 = n1
        geom._n2 = n2
        geom._n = npup
        geom._phase_ab_M1_m = mkP.pad_array(
            geom._phase_ab_M1, geom._n).astype(np.float32)
        wfs._Ntot = npup

        # mod_ampl_pixels = wfs.pyr_ampl / wfs._qpixsize # modulation in pixels

        # Creating field stop mask
        fsradius_pixels = long(wfs.fssize / wfs._qpixsize / 2.)
        if (wfs.fstop == "round"):
            focmask = mkP.dist(
                wfs._Nfft, xc=wfs._Nfft / 2. + 0.5, yc=wfs._Nfft / 2. + 0.5) < (fsradius_pixels)
            # fstop_area = np.pi * (wfs.fssize/2.)**2. #UNUSED
        elif (wfs.fstop == "square"):
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
        wfs._submask = pyr_focmask

        # Creating pyramid mask
        pyrsize = wfs._Nfft
        cobs = tel.cobs
        rpup = geom.pupdiam / 2.0
        dpup = geom.pupdiam
        nrebin = wfs._nrebin
        fracsub = wfs.fracsub
        pup_sep = int(wfs.nxsub)
        # number of pix betw two centers two pupil images

        # pyrsize = 512
        # cobs = 0.12
        # rpup = 64.0
        # dpup = 128.0
        # nrebin = 4
        # fracsub = 0.7
        # pup_sep = 16

        y = np.tile(np.arange(pyrsize) - pyrsize / 2, (pyrsize, 1))
        x = y.T

        # Pangle = int(3*dpup/4)+1 #Pyramid angle
        Pangle = pup_sep * nrebin  # Pyramid angle in HR pixels
        K = 2 * np.pi * Pangle / pyrsize
        pyr = K * (np.abs(x) + np.abs(y))  # Pyramide
        pyr = np.fft.fftshift(pyr)

        wfs._halfxy = pyr

        # Valid pixels identification
        mypup = mkP.dist(
            pyrsize, xc=pyrsize / 2. + 0.5, yc=pyrsize / 2. + 0.5) < (rpup)
        # mypup = geom._ipupil
        # mypyr = np.abs(np.fft.ifft2(np.fft.fft2(mypup)*np.exp(1j*pyr)))**2

        xcc = pyrsize / 2 - Pangle + 0.5
        ycc = pyrsize / 2 + Pangle + 0.5

        mymsk = mkP.dist(pyrsize, xc=xcc, yc=xcc)
        # mypup = (mymsk < rpup)
        mypup = (mymsk < rpup) & (mymsk > cobs * rpup)
        mskreb = rebin(mypup.copy(), [pyrsize / nrebin, pyrsize / nrebin])
        validx = np.where(mskreb >= fracsub)[1].astype(np.int32)
        validy = np.where(mskreb >= fracsub)[0].astype(np.int32)
        nvalid = validx.size
        mskrebtot = mskreb

        mymsk = mkP.dist(pyrsize, xc=ycc, yc=ycc)
        # mypup = (mymsk < rpup)
        mypup = (mymsk < rpup) & (mymsk > cobs * rpup)
        mskreb = rebin(mypup * 1., [pyrsize / nrebin, pyrsize / nrebin])
        tmpx = np.where(mskreb >= fracsub)[1].astype(np.int32)
        tmpy = np.where(mskreb >= fracsub)[0].astype(np.int32)
        validx = np.concatenate((validx, tmpx))
        validy = np.concatenate((validy, tmpy))
        mskrebtot += mskreb

        mymsk = mkP.dist(pyrsize, xc=xcc, yc=ycc)
        # mypup = (mymsk < rpup)
        mypup = (mymsk < rpup) & (mymsk > cobs * rpup)
        mskreb = rebin(mypup * 1., [pyrsize / nrebin, pyrsize / nrebin])
        tmpx = np.where(mskreb >= fracsub)[1].astype(np.int32)
        tmpy = np.where(mskreb >= fracsub)[0].astype(np.int32)
        validx = np.concatenate((validx, tmpx))
        validy = np.concatenate((validy, tmpy))
        mskrebtot += mskreb

        mymsk = mkP.dist(pyrsize, xc=ycc, yc=xcc)
        # mypup = (mymsk < rpup)
        mypup = (mymsk < rpup) & (mymsk > cobs * rpup)
        mskreb = rebin(mypup * 1., [pyrsize / nrebin, pyrsize / nrebin])
        tmpx = np.where(mskreb >= fracsub)[1].astype(np.int32)
        tmpy = np.where(mskreb >= fracsub)[0].astype(np.int32)
        validx = np.concatenate((validx, tmpx))
        validy = np.concatenate((validy, tmpy))
        mskrebtot += mskreb

        wfs._nvalid = nvalid
        wfs._validsubsx = validx
        wfs._validsubsy = validy
        wfs._hrmap = mskrebtot.astype(np.int32)

        if(wfs.pyr_pos == None):
            pixsize = (np.pi * qpixsize) / (3600 * 180)
            scale_fact = 2 * np.pi / npup * \
                (wfs.Lambda / tel.diam / 4.848) / pixsize * wfs.pyr_ampl
#             Proposition de Flo
#             scale_fact = 2 * np.pi / npup * \
#                 (wfs.Lambda * 1e-6 / tel.diam) / pixsize * wfs.pyr_ampl
            cx = scale_fact * \
                np.sin((np.arange(wfs.pyr_npts)) * 2. * np.pi / wfs.pyr_npts)
            cy = scale_fact * \
                np.cos((np.arange(wfs.pyr_npts)) * 2. * np.pi / wfs.pyr_npts)
            # mod_npts = wfs.pyr_npts #UNUSED
        else:
            if(verbose == 0):
                print "Using user-defined positions for the pyramid modulation"
            cx = wfs.pyr_pos[:, 0] / qpixsize
            cy = wfs.pyr_pos[:, 1] / qpixsize
            # mod_npts=cx.shape[0] #UNUSED

        wfs._pyr_cx = cx.astype(np.int32)
        wfs._pyr_cy = cy.astype(np.int32)

        wfs._nphotons = wfs.zerop * \
            2.51189 ** (-wfs.gsmag) * loop.ittime * wfs.optthroughput

        # spatial filtering by the pixel extent:
        # *2/2 intended. min should be 0.40 = sinc(0.5)^2.
        y = np.tile(np.arange(pyrsize) - pyrsize / 2, (pyrsize, 1))
        x = y.T
        x = x * 1.0 / (pyrsize - 1) * 2 / 2
        y = y * 1.0 / (pyrsize - 1) * 2 / 2
        sincar = np.roll(np.sinc(x) * np.sinc(y), x.shape[0] / 2, axis=0)
        sincar = np.roll(sincar, y.shape[0] / 2, axis=1)

        # sincar = np.roll(np.pi*x*np.pi*y,x.shape[1],axis=1)
        wfs._sincar = sincar.astype(np.float32)

        pup = geom._mpupil
        pupreb = rebin(
            pup * 1., [geom._n / nrebin, geom._n / nrebin])
        wsubok = np.where(pupreb >= wfs.fracsub)
        pupvalid = pupreb * 0.
        pupvalid[wsubok] = 1
        wfs._isvalid = pupvalid.astype(np.int32)

        validsubsx = np.where(pupvalid)[1].astype(np.int32)
        validsubsy = np.where(pupvalid)[0].astype(np.int32)

        istart = (
            (np.linspace(0.5, geom._n + 0.5, wfs.nxsub + 2))[:-1]).astype(np.int64)

        jstart = np.copy(istart)
        wfs._istart = istart.astype(np.int32)
        wfs._jstart = jstart.astype(np.int32)

        # sorting out valid subaps
        fluxPerSub = np.zeros((wfs.nxsub+2, wfs.nxsub+2), dtype=np.float32)

        for i in range(wfs.nxsub):
            indi = istart[i] + 1  # +2-1 (yorick->python)
            for j in range(wfs.nxsub):
                indj = jstart[j] + 1  # +2-1 (yorick->python)
                fluxPerSub[i, j] = np.sum(
                    geom._mpupil[indi:indi +  wfs.npix, indj:indj +  wfs.npix])
                # fluxPerSub[i,j] = np.where(geom._mpupil[indi:indi+pdiam,indj:indj+pdiam] > 0)[0].size

        fluxPerSub = fluxPerSub /  wfs.nxsub** 2.

        wfs._fluxPerSub = fluxPerSub
 
        phasemap = np.zeros((wfs.npix, wfs.npix, wfs._nvalid), dtype=np.int32)
        x,y=indices(geom._n) #we need c-like indice
        x-=1
        y-=1
        tmp=x+y*geom._n
        for i in range(wfs._nvalid):
            indi=istart[validsubsx[i]]+1 #+2-1 (yorick->python
            indj=jstart[validsubsy[i]]+1
            phasemap[:,:,i]=tmp[indi:indi+wfs.npix, indj:indj+wfs.npix]

        wfs._phasemap = phasemap

        # x,y= indices(Nfft)
        # x-=(Nfft+1)/2.
        # y-=(Nfft+1)/2.
        # if(wfs.nxsub*nrebin%2==1):
        #    coef1=0.
        # else:
        #    coef1=-0.5
        # if(wfs.nxsub%2==1 and wfs.npix*nrebin%2==1):
        #    coef2=1
        # else:
        #    coef2=0.5

        # pshift = np.exp(1j*2*np.pi*(coef1/Nfft+coef2*nrebin/wfs.npix/Nfft)*(x+y))
        wfs._pyr_offsets = np.array([0], dtype=np.int32)  # pshift
        # this defines a phase for half a pixel shift in x and y
        # defined on mpupil
        # x,y = indices(geom._n)
        # x-=(geom._n+1)/2.
        # y-=(geom._n+1)/2.
        # phase_shift = np.roll( np.exp(1j*2*np.pi*(0.5*(x+y))/geom._n), x.shape[0]/2, axis=0 )
        # phase_shift = np.roll( phase_shift, x.shape[1]/2, axis=1 )

    if(wfs.type_wfs == "sh" or wfs.type_wfs == "geo"):
        # this is the i,j index of lower left pixel of subap
        istart = (
            (np.linspace(0.5, geom.pupdiam + 0.5, wfs.nxsub + 1) + 1)[:-1]).astype(np.int64)

        jstart = np.copy(istart)
        wfs._istart = istart.astype(np.int32)
        wfs._jstart = jstart.astype(np.int32)

        # sorting out valid subaps
        fluxPerSub = np.zeros((wfs.nxsub, wfs.nxsub), dtype=np.float32)

        for i in range(wfs.nxsub):
            indi = istart[i] + 1  # +2-1 (yorick->python)
            for j in range(wfs.nxsub):
                indj = jstart[j] + 1  # +2-1 (yorick->python)
                fluxPerSub[i, j] = np.sum(
                    geom._mpupil[indi:indi + pdiam, indj:indj + pdiam])
                # fluxPerSub[i,j] = np.where(geom._mpupil[indi:indi+pdiam,indj:indj+pdiam] > 0)[0].size

        fluxPerSub = fluxPerSub / pdiam ** 2.

        pupvalid = (fluxPerSub >= wfs.fracsub) * 1
        pupvalid = pupvalid.T
        wfs._isvalid = pupvalid.astype(np.int32)
        wfs._nvalid = np.sum(pupvalid)
        wfs._fluxPerSub = fluxPerSub
        validx = np.where(pupvalid)[1].astype(np.int32)
        validy = np.where(pupvalid)[0].astype(np.int32)
        wfs._validsubsx = validx
        wfs._validsubsy = validy

        # this defines how we cut the phase into subaps
        phasemap = np.zeros((pdiam * pdiam, wfs._nvalid), dtype=np.int32)

        x, y = indices(geom._n)
        x -= 1
        y -= 1
        tmp = x + y * geom._n

        n = wfs._nvalid
        # for i in range(wfs._nvalid):
        for i in range(n):
            indi = istart[wfs._validsubsy[i]] + 1  # +2-1 (yorick->python)
            indj = jstart[wfs._validsubsx[i]] + 1
            phasemap[:, i] = tmp[
                indi:indi + pdiam, indj:indj + pdiam].flatten("C")
        wfs._phasemap = phasemap

        # this is a phase shift of 1/2 pix in x and y
        if(wfs.type_wfs == "sh"):
            halfxy = np.linspace(
                0, 2 * np.pi, wfs._Nfft + 1)[0:wfs._pdiam] / 2.
            halfxy = np.tile(halfxy, (wfs._pdiam, 1))
            halfxy += halfxy.T
            # wfs._halfxy = <float*>(halfxy*0.).data #dont work: half*0 is temp
            # python obj

            if(wfs.npix % 2 == 1 and wfs._nrebin % 2 == 1):
                # wfs._halfxy = <float*>(halfxy*0.)
                halfxy = np.zeros((wfs._pdiam, wfs._pdiam), dtype=np.float32)
                wfs._halfxy = halfxy.astype(np.float32)
            else:
                wfs._halfxy = halfxy.astype(np.float32)

        else:
            halfxy = np.linspace(
                0, 2 * np.pi, wfs._pdiam + 1)[1:wfs._pdiam] / 2
            halfxy = np.tile(halfxy, (wfs._pdiam, 1))
            halfxy = halfxy * 0.
            wfs._halfxy = halfxy.astype(np.float32)

        # must be set even if unused
        wfs._submask = np.array([0], dtype=np.float32)

    if (wfs.type_wfs == "sh"):
        # this defines how we create a larger fov if required
        if(wfs._Ntot != wfs._Nfft):
            indi = long((wfs._Ntot - wfs._Nfft) / 2.)  # +1 -1 (yorick>python)
            indj = long(indi + wfs._Nfft)
            x, y = indices(wfs._Nfft)
            # hrpix
            tmp = np.zeros((wfs._Ntot, wfs._Ntot))
            tmp[indi:indj, indi:indj] = np.roll(
                x + (y - 1) * wfs._Nfft, wfs._Nfft / 2, axis=0)
            tmp[indi:indj, indi:indj] = np.roll(
                tmp[indi:indj, indi:indj], wfs._Nfft / 2, axis=1)
            # hrmap=roll(hrpix)
            tmp = np.roll(tmp, wfs._Ntot / 2, axis=0)
            tmp = np.roll(tmp, wfs._Ntot / 2, axis=1)

            tmp = np.where(tmp.flatten())[0]

            wfs._hrmap = np.copy(tmp.astype(np.int32))
            # must be set even if unused
            wfs._sincar = np.array([0], dtype=np.float32)

        else:
            tmp = np.zeros((1))
            wfs._hrmap = np.copy(tmp.astype(np.int32))
            # must be set even if unused
            wfs._sincar = np.array([0], dtype=np.float32)

        if(wfs._nrebin * wfs.npix % 2 != wfs._Ntot % 2):
            # +2-1 yorick>python
            indi = long((wfs._Ntot - wfs._nrebin * wfs.npix) / 2.) + 1
        else:
            indi = long((wfs._Ntot - wfs._nrebin * wfs.npix) / 2.) + 0
        indj = long(indi + wfs._nrebin * wfs.npix - 1)

        x, y = indices(wfs._nrebin * wfs.npix)
        x = ((x - 1) / wfs._nrebin).astype(np.int64)
        y = ((y - 1) / wfs._nrebin).astype(np.int64)

        # binindices
        binindices = np.zeros((wfs._Ntot, wfs._Ntot))
        binindices[indi:indj + 1, indi:indj + 1] = x + y * wfs.npix + 1

        binmap = np.zeros((wfs._nrebin * wfs._nrebin, wfs.npix * wfs.npix))

        x, y = indices(wfs._Ntot)
        x -= 1
        y -= 1
        tmp = x + y * wfs._Ntot

        if(wfs.gsalt <= 0):
            binindices = np.roll(binindices, binindices.shape[0] / 2, axis=0)
            binindices = np.roll(binindices, binindices.shape[1] / 2, axis=1)

        for i in range(wfs.npix * wfs.npix):
            binmap[:, i] = tmp[np.where(binindices == i + 1)]
        # binmap=np.reshape(binmap.flatten("F"),(binmap.shape[0],binmap.shape[1]))
        wfs._binmap = np.copy(binmap.astype(np.int32))

        dr0 = tel.diam / atmos.r0 * \
            (0.5 / wfs.Lambda) ** 1.2 / np.cos(geom.zenithangle * dtor) ** 0.6
        fwhmseeing = wfs.Lambda / \
            (tel.diam / np.sqrt(wfs.nxsub ** 2. + (dr0 / 1.5) ** 2.)) / 4.848
        kernelfwhm = np.sqrt(fwhmseeing ** 2. + wfs.kernel ** 2.)

        tmp = makegaussian(wfs._Ntot, kernelfwhm / wfs._qpixsize, wfs._Ntot / 2 + 1,
                           wfs._Ntot / 2 + 1).astype(np.float32)

        tmp = np.roll(tmp, tmp.shape[0] / 2, axis=0)
        tmp = np.roll(tmp, tmp.shape[1] / 2, axis=1)

        tmp[0,
            0] = 1.  # this insures that even with fwhm=0, the kernel is a dirac
        tmp = tmp / np.sum(tmp)
        tmp = np.fft.fft2(tmp).astype(np.complex64) / (wfs._Ntot * wfs._Ntot)
        wfs._ftkernel = np.copy(tmp).astype(np.complex64)

        # dealing with photometry
        telSurf = np.pi / 4. * (1 - tel.cobs ** 2.) * tel.diam ** 2.

        # from the guide star
        if(wfs.gsalt == 0):
            if(wfs.zerop == 0):
                wfs.zerop = 1e11
            wfs._nphotons = wfs.zerop * 10 ** (-0.4 * wfs.gsmag) * \
                wfs.optthroughput * \
                (tel.diam / wfs.nxsub) ** 2. / telSurf * \
                loop.ittime
# include throughput to WFS
# for unobstructed subaperture
# per iteration

        else:  # we are dealing with a LGS
            wfs._nphotons = wfs.lgsreturnperwatt * \
                wfs.laserpower * \
                wfs.optthroughput * \
                (tel.diam / wfs.nxsub) ** 2. * 1e4 * \
                loop.ittime
# detected by WFS
# ... for given power
# include throughput to WFS
# for unobstructed subaperture
# per iteration

        if(verbose == 0):
            print "nphotons : ", wfs._nphotons


cdef init_wfs_size(Param_wfs wfs, int n, Param_atmos atmos,
                   Param_tel tel, int psize, long * pdiam, int * Nfft, int * Ntot, int * nrebin,
                   float * pixsize, float * qpixsize, int verbose=0):
    """Compute all the parameters usefull for further WFS image computation (array sizes)

    :parameters:
        wfs: (Param_wfs) : wfs settings

        n: (int) : WFS number

        atmos: (Param_atmos) : atmos settings

        tel: (Param_tel) : telescope settings

        psize: (int) : unused TODO: remove it

        pdiam: (long) : pupil diam for each subap (pixels) (output)

        Nfft: (int*) : fft size for a subap (pixels) (output)

        Ntot: (int) : hr image size for a subap (pixels) (output)

        nrebin: (int*) : rebin factor for a subap (output)

        pixsize: (float*) : simulated pixel size for a subap (arcsec) (output)

        qpixsize: (float*) : quantum pixel size for a subap (arcsec) (output)

        verbose: (int) : (optional) display informations if 0

    Scheme to determine arrays sizes
    sh :
    k = 6
    p = k * d/r0
    n = long(2*d*v/lambda/RASC)+1
    N = fft_goodsize(k*n/v*lambda/r0*RASC)
    u = k * lambda / r0 * RASC / N
    n = v/u - long(v/u) > 0.5 ? long(v/u)+1 : long(v/u)
    v = n * u
    Nt = v * Npix

    pyr :
    Fs = field stop radius in arcsec
    N size of big array for FFT in pixels
    P pupil diameter in pixels
    D diameter of telescope in m
    Nssp : number of pyr measurement points in the pupil

    Rf = Fs . N . D / lambda / P
    ideally we choose : Fs = lambda / D . Nssp / 2

    if we want good sampling of r0 (avoid aliasing of speckles)
    P > D / r0 . m
    with m = 2 or 3

    to get reasonable space between pupil images : N > P.(2 + 3S) 
    with S close to 1
    N must be a power of 2

    to ease computation lets assume that Nssp is a divider of P
    scaling factor between high res pupil images in pyramid model
    and actual pupil size on camera images would be P / Nssp

    """

    cdef float r0 = atmos.r0 * (wfs.Lambda * 2) ** (6. / 5)
    cdef int k, padding
    cdef long nphase, npix,
    cdef float subapdiam

    cdef np.ndarray w, npix_ok

    if(r0 != 0):
        if(verbose == 0):
            print "r0 for WFS :", "%3.2f" % r0, " m"
        # seeing = RASC * (wfs.lambda * 1.e-6) / r0
        if(verbose == 0):
            print "seeing for WFS : ", "%3.2f" % (RASC * (wfs.Lambda * 1.e-6) / r0), "\""

    if(pdiam[0] <= 0):
        # this case is usualy for the wfs with max # of subaps
        # we look for the best compromise between pixsize and fov
        subapdiam = tel.diam / float(wfs.nxsub)  # diam of subap
        k = 6
        pdiam[0] = long(k * subapdiam / r0)  # number of phase points per subap
        if (pdiam[0] < 16):
            pdiam[0] = 16

        # Must be even to keep ssp and actuators grids aligned in the pupil
        if((pdiam[0] * wfs.nxsub) % 2):
            pdiam[0] += 1

        if(wfs.type_wfs == "sh"):
            nrebin[0] = long(
                2 * subapdiam * wfs.pixsize / (wfs.Lambda * 1.e-6) / RASC) + 1
            nrebin[0] = max(2, nrebin[0])
            # first atempt on a rebin factor

            # since we clipped pdiam we have to be carreful in nfft computation
            Nfft[0] = fft_goodsize(long(
                pdiam[0] / subapdiam * nrebin[0] / wfs.pixsize * RASC * (wfs.Lambda * 1.e-6)))
            # size of the support in fourier domain

            # qpixsize = k * (wfs.Lambda*1.e-6) / r0 * RASC / Nfft
            qpixsize[0] = (
                pdiam[0] * (wfs.Lambda * 1.e-6) / subapdiam * RASC) / Nfft[0]

        if(wfs.type_wfs=="pyr" or wfs.type_wfs == "roof"):
            # while (pdiam % wfs.npix != 0) pdiam+=1;
            padding = 2
            nphase = pdiam[0] * wfs.nxsub + 2 * padding * pdiam[0]
            qpixsize[0] = (
                pdiam[0] * wfs.nxsub / float(nphase)) * wfs.Lambda / tel.diam / 4.848

            fssize_pixels = long(wfs.fssize / qpixsize[0] / 2.)
            # nrebin  = pdiam / wfs.npix

            npix = (wfs.nxsub + 2 * padding)
            npix_ok = npix * (np.arange(100) + 1)
            npix_ok = npix_ok[np.where(npix_ok.flatten() <= (nphase / 2.))[0]]

            # now we want the quadrant image dim to be > fssize_pixels:
            w = np.where(npix_ok.flatten() >= fssize_pixels)[0]
            if (w.size == 0):
                # maxfs = npix_ok[0]*2*psize
                msg = "wfs ", n, ". ffsize too large "  # (max=",maxfs,")!"
                raise ValueError(msg)

            # npix_ok = npix_ok[w[1]]

            nrebin[0] = npix_ok[w[0]] / npix
            Nfft[0] = npix_ok[w[0]]
            Ntot[0] = nphase
            pixsize[0] = qpixsize[0] * nrebin[0]

        if(wfs.type_wfs == "pyrhr"):
            # while (pdiam % wfs.npix != 0) pdiam+=1;
            k = 4
            pdiam[0] = long(tel.diam / r0 * k)
            while (pdiam[0] % wfs.nxsub != 0):
                pdiam[0] += 1  # we choose to have a multiple of wfs.nxsub

            m = 3
            # fft_goodsize( m * pdiam[0])
            Nfft[0] = long(2 ** np.ceil(np.log2(m * pdiam[0])))

            nrebin[0] = pdiam[0] / wfs.nxsub
            while(Nfft[0] % nrebin[0] != 0):
                nrebin[0] += 1  # we choose to have a divisor of Nfft
                pdiam[0] = nrebin[0] * wfs.nxsub
                nphase = pdiam[0] + 2 * padding
                Nfft[0] = long(2 ** np.ceil(np.log2(m * pdiam[0])))

            qpixsize[0] = (
                pdiam[0] * (wfs.Lambda * 1.e-6) / tel.diam * RASC) / Nfft[0]

            padding = 2
            nphase = pdiam[0] + 2 * padding

            fssize_pixels = long(wfs.fssize / qpixsize[0] / 2.)

            Ntot[0] = Nfft[0] / pdiam[0] * wfs.nxsub
            pixsize[0] = qpixsize[0] * nrebin[0]
            pdiam[0] = pdiam[0] / wfs.nxsub

        # quantum pixel size
    else:
        # this case is for a wfs with fixed # of phase points
        subapdiam = tel.diam / float(wfs.nxsub)  # diam of subap
        if (wfs.type_wfs != "geo"):
            Nfft[0] = fft_goodsize(2 * pdiam[0])
            # size of the support in fourier domain

            qpixsize[0] = pdiam[0] * \
                (wfs.Lambda * 1.e-6) / subapdiam * RASC / Nfft[0]
            # quantum pixel size

    if (wfs.type_wfs == "sh"):
        # actual rebin factor
        if(wfs.pixsize / qpixsize[0] - long(wfs.pixsize / qpixsize[0]) > 0.5):
            nrebin[0] = long(wfs.pixsize / qpixsize[0]) + 1
        else:
            nrebin[0] = long(wfs.pixsize / qpixsize[0])

        # actual pixel size
        pixsize[0] = nrebin[0] * qpixsize[0]

        if (pixsize[0] * wfs.npix > qpixsize[0] * Nfft[0]):
            Ntot[0] = fft_goodsize(
                long(pixsize[0] * wfs.npix / qpixsize[0]) + 1)
        else:
            Ntot[0] = Nfft[0]

        if (Ntot[0] % 2 != Nfft[0] % 2):
            Ntot[0] += 1

    if (wfs.type_wfs == "geo"):
        if(verbose == 0):
            print "quantum pixsize : ", "%5.4f" % qpixsize[0], "\""
        if(verbose == 0):
            print "simulated FoV : ", "%3.2f" % (Ntot[0] * qpixsize[0]), "\" x ", "%3.2f" % (Ntot[0] * qpixsize[0]), "\""
        if(verbose == 0):
            print "actual pixsize : ", "%5.4f" % pixsize[0]
        if(verbose == 0):
            print "actual FoV : ", "%3.2f" % (pixsize[0] * wfs.npix), "\" x ", "%3.2f" % (pixsize[0] * wfs.npix), "\""
        if(verbose == 0):
            print "number of phase points : ", pdiam[0]
        if(verbose == 0):
            print "size of fft support : ", Nfft[0]
        if (Ntot[0] > Nfft[0]):
            if(verbose == 0):
                print "size of HR spot support : ", Ntot[0]

    if (wfs.type_wfs == "pyrhr"):
        if(verbose == 0):
            print "quantum pixsize in pyr image : ", "%5.4f" % qpixsize[0], "\""
        if(verbose == 0):
            print "simulated FoV : ", "%3.2f" % (Nfft[0] * qpixsize[0]), "\" x ", "%3.2f" % (Nfft[0] * qpixsize[0]), "\""
        # if(verbose==0):print "actual pixsize : ", "%5.4f"%pixsize[0]
        if(verbose == 0):
            print "number of phase points : ", pdiam[0] * wfs.nxsub
        if(verbose == 0):
            print "size of fft support : ", Nfft[0]
        # if(verbose==0):print "size of pyramid images : ",Ntot[0]

cpdef noise_cov(int nw, Param_wfs p_wfs, Param_atmos p_atmos, Param_tel p_tel):
    """Compute the diagonal of the noise covariance matrix for a SH WFS (arcsec^2)
    Photon noise: (pi^2/2)*(1/Nphotons)*(d/r0)^2 / (2*pi*d/lambda)^2
    Electronic noise: (pi^2/3)*(wfs.noise^2/N^2photons)*wfs.npix^2*(wfs.npix*wfs.pixsize*d/lambda)^2 / (2*pi*d/lambda)^2

    :parameters:
        nw: wfs number
        p_wfs: (Param_wfs) : wfs settings
        p_atmos: (Param_atmos) : atmos settings
        p_tel: (Param_tel) : telescope settings
    :return:
        cov : (np.ndarray(ndim=1,dtype=np.float64)) : noise covariance diagonal

    """
    cov = np.zeros(2 * p_wfs._nvalid)
    if(p_wfs.noise >= 0):
        m = p_wfs._validsubsy
        n = p_wfs._validsubsx
        ind = m * p_wfs.nxsub + n
        flux = np.copy(p_wfs._fluxPerSub)
        flux = flux.reshape(flux.size, order='F')
        flux = flux[ind]
        Nph = flux * p_wfs._nphotons

        r0 = (p_wfs.Lambda / 0.5) ** (6.0 / 5.0) * p_atmos.r0

        sig = (np.pi ** 2 / 2) * (1 / Nph) * \
            (1. / r0) ** 2  # Photon noise in m^-2
        # Noise variance in rad^2
        sig /= (2 * np.pi / (p_wfs.Lambda * 1e-6)) ** 2
        sig *= RASC ** 2

        Ns = p_wfs.npix  # Number of pixel
        Nd = (p_wfs.Lambda * 1e-6) * RASC / p_wfs.pixsize
        sigphi = (np.pi ** 2 / 3.0) * (1 / Nph ** 2) * (p_wfs.noise) ** 2 * \
            Ns ** 2 * (Ns / Nd) ** 2  # Phase variance in m^-2
        # Noise variance in rad^2
        sigsh = sigphi / (2 * np.pi / (p_wfs.Lambda * 1e-6)) ** 2
        sigsh *= RASC ** 2  # Electronic noise variance in arcsec^2

        cov[:len(sig)] = sig + sigsh
        cov[len(sig):] = sig + sigsh

    return cov


cpdef fft_goodsize(long s):
    """find best size for a fft from size s

    :parameters:
         s: (long) size
    """
    return 2 ** (long(np.log2(s)) + 1)


cpdef bin2d(np.ndarray data_in, int binfact):
    """
    Returns the input 2D array "array", binned with the binning factor "binfact".
    The input array X and/or Y dimensions needs not to be a multiple of
    "binfact"; The final/edge pixels are in effect replicated if needed.
    This routine prepares the parameters and calls the C routine _bin2d.
    The input array can be of type long, float or double.
    Last modified: Dec 15, 2003.
    Author: F.Rigaut
    SEE ALSO: _bin2d

    :parmeters:
        data_in: (np.ndarray) : data to binned

        binfact: (int) : binning factor

    """
    if(binfact < 1):
        raise ValueError("binfact has to be >= 1")

    cdef int nx, ny, fx, fy
    nx = data_in.shape[0]
    ny = data_in.shape[1]
    fx = int(np.ceil(nx / float(binfact)))
    fy = int(np.ceil(ny / float(binfact)))

    cdef np.ndarray data_out = np.zeros((fx, fy), dtype=data_in.dtype)

    cdef int i, j, i1, i2, j1, j2

    for i1 in range(fx):
        for j1 in range(fy):
            for i2 in range(binfact):
                for j2 in range(binfact):
                    i = i1 * binfact + i2
                    j = j1 * binfact + j2
                    if(i >= nx):
                        i = nx - 1
                    if(j >= ny):
                        j = ny - 1
                    data_out[i1, j1] += data_in[i, j]

    return data_out

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
    Sensors(nsensors,type_data,npup,nxsub,nvalid,nphase,pdiam,npix,nrebin,nfft,nftota,nphot,lgs,odevice,comm_size,rank)

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

    def __cinit__(self, int nsensors,
                  Telescope tel,
                  type_data,
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

        cdef char ** type_wfs = < char ** > malloc(len(type_data) * sizeof(char * ))
        cdef int i
        for i in range(nsensors):
            type_wfs[i] = type_data[i]
        cdef carma_context * context = carma_context.instance()
        if odevice < 0:
            odevice = context.get_activeDevice()
        if(type_data == "geo"):
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
                         np.ndarray[ndim=1, dtype=np.int64_t] seed):
        """Call the function sensors_initgs

    :parameters:
        xpos: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        ypos: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        Lambda: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        mag: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        zerop: (float) :

        size:  (np.ndarray[ndim=1,dtype=np.int64_t  ]) :

        noise: (np.ndarray[ndim=1,dtype=np.float32_t]) :

        seed:  (np.ndarray[ndim=1,dtype=np.int64_t  ]) :
        """

        if(noise.size == 0):
            self.sensors.sensors_initgs( < float * > xpos.data, < float * > ypos.data,
                                        < float * > Lambda.data, < float * > mag.data, zerop, < long * > size.data)

        elif(seed.size == 0):
            self.sensors.sensors_initgs( < float * > xpos.data, < float * > ypos.data,
                                        < float * > Lambda.data, < float * > mag.data, zerop, < long * > size.data,
                                        < float * > noise.data)
        else:
            self.sensors.sensors_initgs( < float * > xpos.data, < float * > ypos.data,
                                        < float * > Lambda.data, < float * > mag.data, zerop, < long * > size.data,
                                        < float * > noise.data, < long * > seed.data)

    cpdef sensors_initarr(self, int n, Param_wfs wfs):
        """Call the function wfs_initarrays from a sutra_wfs of the Sensors

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

        cdef np.ndarray[dtype = np.int32_t] phasemap_F = wfs._phasemap.flatten("F")
        cdef np.ndarray[dtype = np.int32_t] hrmap_F = wfs._hrmap.flatten("F")
        cdef np.ndarray[dtype = np.int32_t] binmap_F
        cdef int * binmap = NULL
        if(self.sensors.d_wfs[n].type == < bytes > "sh"):
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
        cdef int * cx = < int * > wfs._pyr_cx.data
        cdef int * cy = < int * > wfs._pyr_cy.data

        """
        #type depend on wfs type
        cdef np.ndarray halfxy_F=wfs._halfxy.flatten("F")
        # [dtype=np.float32_t] halfxy_F=
        cdef float *halfxy=<float*>halfxy_F.data
        #halfxy_F=wfs._halfxy.flatten("F")
        #halfxy=<float*>halfxy_F.data
        """
        cdef np.ndarray halfxy_F  # [dtype=np.float32_t] halfxy_F=wfs._halfxy.flatten("F")
        cdef float * halfxy  # =<float*>halfxy_F.data
        # if(self.sensors.d_wfs[n].type=="sh"):
        halfxy_F = wfs._halfxy.flatten("F")
        halfxy = < float * > halfxy_F.data

        cdef np.ndarray[dtype = np.float32_t] submask_F
        submask_F = wfs._submask.flatten("F").astype(np.float32)
        cdef np.ndarray[dtype = np.float32_t] sincar_F = wfs._sincar.flatten("F")

        cdef float * submask = < float * > submask_F.data
        cdef float * sincar = < float * > sincar_F.data

        cdef np.ndarray[dtype = np.complex64_t]ftkernel_F
        cdef float * ftkernel = NULL
        if(self.sensors.d_wfs[n].type == < bytes > "sh"):
            ftkernel_F = wfs._ftkernel.flatten("F")
            ftkernel = < float * > ftkernel_F.data

        cdef float * fluxPerSub = NULL
        cdef np.ndarray tmp
        cdef np.ndarray tmp_offset
        cdef cuFloatComplex * offset = NULL

        cdef np.ndarray tmp_halfxy
        cdef cuFloatComplex * pyr_halfxy = NULL

        # swap validx and validy due to transposition from yorick to python
        if(self.sensors.d_wfs[n].type == < bytes > "geo"):
            tmp = wfs._fluxPerSub.T[np.where(wfs._isvalid > 0)].copy()
            fluxPerSub = < float * > tmp.data
            wfs_geom = dynamic_cast_wfs_geom_ptr(self.sensors.d_wfs[n])
            wfs_geom.wfs_initarrays(phasemap, halfxy, fluxPerSub,
                                    validx, validy)

        elif(self.sensors.d_wfs[n].type== < bytes >"pyr"):
            tmp_offset=wfs._pyr_offsets.flatten("F").copy()
            offset=<cuFloatComplex*>tmp_offset.data
            wfs_pyr = dynamic_cast_wfs_pyr_pyr4_ptr(self.sensors.d_wfs[n])
            wfs_pyr.wfs_initarrays(<cuFloatComplex*>halfxy, offset,
                    submask, cx, cy,
                    sincar, phasemap, validx,validy)

        elif(self.sensors.d_wfs[n].type == < bytes > "pyrhr"):
            tmp = wfs._fluxPerSub.T[np.where(wfs._isvalid > 0)].copy()
            fluxPerSub = < float * > tmp.data
            tmp_halfxy = np.exp(
                1j * wfs._halfxy).astype(np.complex64).flatten("F").copy()
            pyr_halfxy = < cuFloatComplex * > tmp_halfxy.data
            wfs_pyrhr = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
            wfs_pyrhr.wfs_initarrays(pyr_halfxy, cx, cy, sincar, validx, validy, phasemap, fluxPerSub)

        elif(self.sensors.d_wfs[n].type == < bytes > "roof"):
            tmp_offset = wfs.__pyr_offsets.flatten("F").copy()
            offset = < cuFloatComplex * > tmp_offset.data
            wfs_roof = dynamic_cast_wfs_pyr_roof_ptr(self.sensors.d_wfs[n])
            wfs_roof.wfs_initarrays( < cuFloatComplex * > halfxy, offset,
                                    submask, cx, cy,
                                    sincar, phasemap, validx, validy)

        elif(self.sensors.d_wfs[n].type == < bytes > "sh"):
            tmp = wfs._fluxPerSub.T[np.where(wfs._isvalid > 0)].copy()
            fluxPerSub = < float * > tmp.data
            wfs_sh = dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[n])
            wfs_sh.wfs_initarrays(phasemap, hrmap, binmap, halfxy,
                                  fluxPerSub, validx, validy,
                                  istart, jstart, < cuFloatComplex * > ftkernel)

    cpdef sensors_addlayer(self, int i, bytes type_dm, float alt,
                           float xoff, float yoff):
        """Call function add_layer from the sutra_source of a sutra_wfs of the Sensors

        :parameters:
            i: (int) :

            type_dm: (string) :

            alt: (float) :

            xoff: (float) :

            yoff: (float) :
        """

        cdef carma_context * context = carma_context.instance()
        context.set_activeDevice(self.sensors.device, 1)
        self.sensors.d_wfs[i].d_gs.add_layer( < char * > type_dm, alt, xoff, yoff)

    def sensors_compimg(self, int n):
        """TODO doc

        :param n: (in) : index of the wfs
        """
        cdef carma_context * context = carma_context.instance()
        context.set_activeDeviceForCpy(self.sensors.device, 1)
        self.sensors.d_wfs[n].comp_image()

    def get_offsets(self, int n):
        """Return the 'offset' array of a given wfs

        :param n: (int) : number of the wfs to get the 'offset' from
        """
        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        img = self.sensors.d_wfs[n].d_offsets
        cdims = img.getDims()
        data = np.empty((cdims[1]), dtype=np.float32)
        img.device2host( < float * > data.data)
        return data

    def get_subsum(self, int n):
        """Return the 'subsum' array of a given wfs

        :param n: (int) : number of the wfs to get the 'subsum' from
        """
        cdef carma_obj[float] * subsum
        cdef const long * cdims
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data
        subsum = self.sensors.d_wfs[n].d_subsum
        cdims = subsum.getDims()
        data = np.empty((cdims[1]), dtype=np.float32)
        subsum.device2host( < float * > data.data)
        return data

    def get_imgtele(self, int n, Telescope tel=None, Atmos atmos=None, Dms dms=None):
        """Return the 'image_telemetry' array of a given wfs

        :param n: (int) : number of the wfs to get the 'image_telemetry' from

        :options for raw image computation
            tel (Telescope) : shesha telescope
            atmos (Atmos) : shesha atmos
            dms (Dms) : shesha dms
        """
        cdef carma_host_obj[float] * img
        cdef sutra_wfs_sh * wfs = dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[n])
        cdef const long * cdims
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data

        if (tel and atmos and dms):
            self.sensor_trace(n, "all", tel, atmos, dms)
            self.sensor_compimg(n)
        img = self.sensors.d_wfs[n].image_telemetry
        cdims = img.getDims()
        data = np.empty((cdims[1], cdims[2]), dtype=np.float32)

        wfs.fill_binimage(1)
        img.fill_into( < float * > data.data)
        data[np.where(data < 0)] = 0
        return data

    cpdef get_binimg(self, int n, Telescope tel=None, Atmos atmos=None, Dms dms=None):
        """Return the 'binimg' array of a given wfs

        :param
            n: (int) :number of the wfs to get the 'binimg' from
        :options for raw image computation
            tel (Telescope) : shesha telescope
            atmos (Atmos) : shesha atmos
            dms (Dms) : shesha dms
        """
        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        if (tel and atmos and dms):
            self.sensor_trace(n, "all", tel, atmos, dms)
            self.sensor_compimg(n)
        img = self.sensors.d_wfs[n].d_binimg
        cdims = img.getDims()
        data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
        data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
        dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[n]).fill_binimage(0)
        img.device2host( < float * > data_F.data)

        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        data[np.where(data < 0)] = 0
        return data

    def get_pyrimg(self, int n):
        """Return the image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from

        """

        return self._get_pyrimg(n)

    cdef _get_pyrimg(self, int n):

        """Return the image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from

        """
        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F

        cdef np.ndarray[ndim=3,dtype=np.float32_t] bincube
        cdef np.ndarray[ndim=2,dtype=np.float32_t] pyrimg
        cdef bytes type_wfs = < bytes > self.sensors.d_wfs[n].type
        cdef int npix

        if(type_wfs == "pyrhr"):
            img = self.sensors.d_wfs[n].d_binimg
            cdims = img.getDims()
            data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
            data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
            img.device2host( < float * > data_F.data)
            data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return data

        if(type_wfs == "pyr"):
            bincube = self.get_bincube(n)
            npix = bincube.shape[1]
            pyrimg = np.zeros((2*npix+3,2*npix+3),dtype=np.float32)
            pyrimg[1:npix+1,1:npix+1] = bincube[:,:,0]
            pyrimg[npix+2:2*npix+2,1:npix+1] = bincube[:,:,1]
            pyrimg[1:npix+1,npix+2:2*npix+2] = bincube[:,:,2]
            pyrimg[npix+2:2*npix+2,npix+2:2*npix+2] = bincube[:,:,3]

            return pyrimg

        else:
            raise TypeError("wfs should be a pyr")

    def get_pyrimghr(self, int n):
        """Return the high res image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from

        """

        return self._get_pyrimghr(n)

    cdef _get_pyrimghr(self, int n):

        """Return the high res image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from

        """
        cdef carma_obj[float] * img
        cdef const long * cdims
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F

        cdef bytes type_wfs = < bytes > self.sensors.d_wfs[n].type
        cdef sutra_wfs_pyr_pyrhr * wfs

        if(type_wfs == "pyrhr"):
            wfs = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
            img = wfs.d_hrimg
            cdims = img.getDims()
            data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
            data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
            img.device2host( < float * > data_F.data)
            data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
            return data

        else:
            raise TypeError("wfs should be a pyrhr")

    cpdef comp_modulation(self, int n, int cpt):

        """Return the high res image of a pyr wfs

        :param n: (int) : number of the wfs to get the image from

        """

        cdef bytes type_wfs = < bytes > self.sensors.d_wfs[n].type
        cdef sutra_wfs_pyr_pyrhr * wfs

        if(type_wfs == "pyrhr"):
            wfs = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
            wfs.comp_modulation(cpt)

        else:
            raise TypeError("wfs should be a pyrhr")

    cdef _get_bincube(self, int n):
        """Return the 'bincube' array of a given wfs

        :param n: (int) : number of the wfs to get the 'bincube' from
        """
        cdef carma_obj[float] * cube
        cdef const long * cdims
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] data_F
        cube = self.sensors.d_wfs[n].d_bincube
        cdims = cube.getDims()
        data = np.empty((cdims[1], cdims[2], cdims[3]), dtype=np.float32)
        data_F = np.empty((cdims[3], cdims[2], cdims[1]), dtype=np.float32)
        cube.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
        return data

    def get_bincube(self, int n):
        """Return the 'bincube' array of a given wfs

        :param n: (int) : number of the wfs to get the 'bincube' from
        """
        return self._get_bincube(n)

    cpdef set_bincube(self, int n, np.ndarray[ndim=3, dtype=np.float32_t] data):
        """ Set the bincube of the WFS numner n

        :parameters:
            n: (int) : WFS number
            data: (np.ndarray[ndim=3,dtype=np.float32_t]) : bincube to use
        """
        cdef carma_context * context = carma_context.instance()
        context.set_activeDeviceForCpy(self.sensors.device, 1)
        cdef np.ndarray[dtype = np.float32_t] data_F = data.flatten("F")

        self.sensors.d_wfs[n].d_bincube.host2device( < float * > data_F.data)

    cpdef get_bincubeNotNoisy(self, int n):
        """Return the 'bincube_not_noisy' array of a given wfs. It's the bincube
        before noise has been added

        :param n: (int) : number of the wfs to get the 'bincube_not_noisy' from
        """
        cdef carma_obj[float] * cube
        cdef const long * cdims
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] data_F
        if(self.sensors.error_budget):
            cube = self.sensors.d_wfs[n].d_bincube_notnoisy
            cdims = cube.getDims()
            data = np.empty((cdims[1], cdims[2], cdims[3]), dtype=np.float32)
            data_F = np.empty((cdims[3], cdims[2], cdims[1]), dtype=np.float32)
            cube.device2host( < float * > data_F.data)
            data = np.reshape(
                data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
            return data
        else:
            raise TypeError("the error budget analysis has to be enabled")

    def reset_phase(self, int n):
        """Reset the phase's array of a given wfs

        :param n: (int) : index of the given wfs
        """
        cdef carma_obj[float] * phase
        phase = self.sensors.d_wfs[n].d_gs.d_phase.d_screen
        phase.reset()

    def get_phase(self, int n):
        """Return the phase array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        """
        cdef carma_obj[float] * phase
        cdef const long * cdims
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        phase = self.sensors.d_wfs[n].d_gs.d_phase.d_screen
        cdims = phase.getDims()
        data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
        data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
        phase.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        return data

    def get_lgskern(self, int n):
        """Return the lgskern array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        """
        cdef carma_obj[float] * lgs_kern
        cdef const long * cdims
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 3, dtype = np.float32_t] data_F
        if(self.sensors.d_wfs[n].lgs):
            lgs_kern = self.sensors.d_wfs[n].d_gs.d_lgs.d_lgskern
            cdims = lgs_kern.getDims()
            data = np.empty((cdims[1], cdims[2], cdims[3]), dtype=np.float32)
            data_F = np.empty((cdims[3], cdims[2], cdims[1]), dtype=np.float32)
            lgs_kern.device2host( < float * > data_F.data)
            data = np.reshape(
                data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
            return data
        else:
            raise TypeError("the WFS should be a LGS")

    def get_ftlgskern(self, int n):
        """Return the ftlgskern array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        """
        cdef carma_obj[cuFloatComplex] * ftlgs_kern
        cdef const long * cdims
        cdef np.ndarray[ndim = 3, dtype = np.complex64_t] data
        cdef np.ndarray[ndim = 3, dtype = np.complex64_t] data_F
        if(self.sensors.d_wfs[n].lgs):
            ftlgs_kern = self.sensors.d_wfs[n].d_gs.d_lgs.d_ftlgskern
            cdims = ftlgs_kern.getDims()
            data = np.empty((cdims[1], cdims[2], cdims[3]), dtype=np.complex64)
            data_F = np.empty(
                (cdims[3], cdims[2], cdims[1]), dtype=np.complex64)
            ftlgs_kern.device2host( < cuFloatComplex * > data_F.data)
            data = np.reshape(
                data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
            return data
        else:
            raise TypeError("the WFS should be a LGS")

    def set_phase(self, int n, np.ndarray[ndim=2, dtype=np.float32_t] data):
        """Set the phase array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        :param data: (np.ndarray) : the phase to set
        """
        # self.context.set_activeDeviceForCpy(self.device)
        cdef sutra_source * src = self.sensors.d_wfs[n].d_gs

        cdef np.ndarray[dtype = np.float32_t] data_F = data.flatten("F")

        src.d_phase.d_screen.host2device( < float * > data_F.data)

    def get_camplipup(self, int n):
        """Return the 'camplipup' array of a given wfs

        :param n: (int) : number of the wfs to get the 'camplipup' from

        """
        cdef carma_obj[cuFloatComplex] * amplipup
        cdef const long * cdims

        cdef np.ndarray[ndim = 3, dtype = np.complex64_t] data
        cdef np.ndarray[ndim = 3, dtype = np.complex64_t] data_F
        amplipup = self.sensors.d_wfs[n].d_camplipup
        cdims = amplipup.getDims()

        data = np.zeros((cdims[1], cdims[2], cdims[3]), dtype=np.complex64)
        data_F = np.zeros((cdims[3], cdims[2], cdims[1]), dtype=np.complex64)
        amplipup.device2host( < cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))

        return data

    def get_camplipup_pyr(self, int n):
        """Return the 'camplipup' array of a given wfs in the pyr case

        :param n: (int) : number of the wfs to get the 'camplipup' from

        """
        cdef carma_obj[cuFloatComplex] * amplipup
        cdef const long * cdims

        cdef np.ndarray[ndim = 2, dtype = np.complex64_t] data
        cdef np.ndarray[ndim = 2, dtype = np.complex64_t] data_F
        amplipup = self.sensors.d_wfs[n].d_camplipup
        cdims = amplipup.getDims()

        data = np.zeros((cdims[1], cdims[2]), dtype=np.complex64)
        data_F = np.zeros((cdims[2], cdims[1]), dtype=np.complex64)
        amplipup.device2host( < cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))

        return data

    def get_amplifoc(self, int n):
        """Return the 'amplifoc' array of a given wfs

        :param n: (int) : number of the wfs to get the 'amplifoc' from
        """
        cdef carma_obj[cuFloatComplex] * amplifoc
        cdef const long * cdims
        cdef np.ndarray[ndim = 3, dtype = np.complex64_t] data
        cdef np.ndarray[ndim = 3, dtype = np.complex64_t] data_F
        amplifoc = self.sensors.d_wfs[n].d_camplifoc
        cdims = amplifoc.getDims()
        data = np.zeros((cdims[1], cdims[2], cdims[3]), dtype=np.complex64)
        data_F = np.zeros((cdims[3], cdims[2], cdims[1]), dtype=np.complex64)
        amplifoc.device2host( < cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2], cdims[3]))
        return data

    def get_amplifoc_pyr(self, int n):
        """Return the 'amplifoc' array of a given wfs

        :param n: (int) : number of the wfs to get the 'amplifoc' from
        """
        cdef carma_obj[cuFloatComplex] * amplifoc
        cdef const long * cdims
        cdef np.ndarray[ndim = 2, dtype = np.complex64_t] data
        cdef np.ndarray[ndim = 2, dtype = np.complex64_t] data_F
        amplifoc = self.sensors.d_wfs[n].d_camplifoc
        cdims = amplifoc.getDims()
        data = np.zeros((cdims[1], cdims[2]), dtype=np.complex64)
        data_F = np.zeros((cdims[2], cdims[1]), dtype=np.complex64)
        amplifoc.device2host( < cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        return data

    def get_fttotim_pyr(self, int n):
        """Return the 'fttotim' array of a given wfs

        :param n: (int) : number of the wfs to get the 'amplifoc' from
        """
        cdef carma_obj[cuFloatComplex] * fttotim
        cdef const long * cdims
        cdef np.ndarray[ndim = 2, dtype = np.complex64_t] data
        cdef np.ndarray[ndim = 2, dtype = np.complex64_t] data_F
        fttotim = self.sensors.d_wfs[n].d_fttotim
        cdims = fttotim.getDims()
        data = np.zeros((cdims[1], cdims[2]), dtype=np.complex64)
        data_F = np.zeros((cdims[2], cdims[1]), dtype=np.complex64)
        fttotim.device2host( < cuFloatComplex * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        return data

    def get_hrimg_pyr(self, int n):
        """Return the phase array of a given wfs

        :param n: (int) : number of the wfs to get the phase from
        """
        cdef carma_obj[float] * hrimg
        cdef const long * cdims
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data
        cdef np.ndarray[ndim = 2, dtype = np.float32_t] data_F
        wfs = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[n])
        hrimg = wfs.d_hrimg
        # hrimg=self.sensors.d_wfs[n].d_hrimg
        cdims = hrimg.getDims()
        data = np.empty((cdims[1], cdims[2]), dtype=np.float32)
        data_F = np.empty((cdims[2], cdims[1]), dtype=np.float32)
        hrimg.device2host( < float * > data_F.data)
        data = np.reshape(data_F.flatten("F"), (cdims[1], cdims[2]))
        return data

    cdef _get_slopesDims(self, int n):
        """return the dimension of the slopes array of a given wfs

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
        cdef carma_context * context = carma_context.instance()
        context.set_activeDevice(self.sensors.device, 1)
        cdims = self.sensors.d_wfs[n].d_slopes.getDims()
        dim_tot = cdims[1]

        IF USE_MPI:
            mpi.MPI_Allreduce(mpi.MPI_IN_PLACE, & dim_tot, 1, mpi.MPI_LONG, mpi.MPI_SUM, mpi.MPI_COMM_WORLD)
        return dim_tot

    def get_slopes(self, int n):
        """Return the 'slopes' array of a given wfs

        :param n: (int) : number of the wfs to get the 'slopes' from
        """

        return self._get_slopes(n)

    cdef _get_slopes(self, int n):
        """Return the 'slopes' array of a given wfs

        :param n: (int) : number of the wfs to get the 'slopes' from
        """
        cdef carma_obj[float] * slopes
        cdef const long * cdims
        cdef np.ndarray[ndim = 1, dtype = np.float32_t] data

        cdef carma_context * context = carma_context.instance()
        context.set_activeDevice(self.sensors.device, 1)

        slopes = self.sensors.d_wfs[n].d_slopes
        cdims = slopes.getDims()
        data = np.empty((cdims[1]), dtype=np.float32)
        slopes.device2host( < float * > data.data)

        IF USE_MPI:
            cdef int comm_size, rank
            mpi.MPI_Comm_size(mpi.MPI_COMM_WORLD, & comm_size)
            mpi.MPI_Comm_rank(mpi.MPI_COMM_WORLD, & rank)

            cdef int d = < int > (cdims[1] / 2)

            cdef int * count = < int * > malloc(comm_size * sizeof(int))
            mpi.MPI_Allgather( & d, 1, mpi.MPI_INT, count, 1, mpi.MPI_INT, mpi.MPI_COMM_WORLD)

            cdef int * disp = < int * > malloc((comm_size + 1) * sizeof(int))
            cdef int i, nvalid2
            disp[0] = 0
            for i in range(comm_size):
                disp[i + 1] = disp[i] + count[i]

            cdef np.ndarray[ndim = 1, dtype = np.float32_t] all_slopes = np.zeros(disp[comm_size] * 2,
                                                                                  dtype=np.float32)

            cdef float * send = < float * > data.data
            cdef float * recv = < float * > all_slopes.data

            mpi.MPI_Allgatherv(send, count[rank], mpi.MPI_FLOAT,
                               recv, count, disp,
                               mpi.MPI_FLOAT, mpi.MPI_COMM_WORLD)

            mpi.MPI_Allgatherv( & send[count[rank]], count[rank], mpi.MPI_FLOAT,
                               & recv[disp[comm_size]], count, disp,
                               mpi.MPI_FLOAT, mpi.MPI_COMM_WORLD)

            return all_slopes

        ELSE:
            return data

    cpdef slopes_geom(self, int nsensor, int t):
        """Compute the geometric slopes in a sutra_wfs object

        :parameters:
            nsensor: (int) : wfs number

            :param t: (int) : method (0 or 1)
        """
        cdef sutra_wfs_sh * wfs_sh = NULL
        cdef sutra_wfs_pyr_pyrhr * wfs_pyrhr = NULL

        if( < bytes > self.sensors.d_wfs[nsensor].type == "sh"):
            #raise TypeError("wfs should be a SH")
            wfs_sh = dynamic_cast_wfs_sh_ptr(self.sensors.d_wfs[nsensor])
            wfs_sh.slopes_geom(t)
        else: 
            if( < bytes > self.sensors.d_wfs[nsensor].type == "pyrhr"):
                wfs_pyrhr = dynamic_cast_wfs_pyr_pyrhr_ptr(self.sensors.d_wfs[nsensor])
                wfs_pyrhr.slopes_geom(t)
            else:
                raise TypeError("wfs should be a SH or PYRHR")


    cpdef sensors_trace(self, int n, str type_trace, Telescope tel=None, Atmos atmos=None, Dms dms=None, int rst=0):
        """ Does the raytracing for the wfs phase screen in sutra_wfs

        :parameters:
            n: (int) :

            type_trace: (str) : "all" : raytracing across atmos and dms seen
                                "dm"  : raytracing across dms seen only
                                "atmos" : raytracing across atmos only

            tel: (Telescope) :(optional) Telescope object

            atmos: (Atmos) :(optional) Atmos object

            dms: (Dms) : (optional) Dms object

            rst: (int) : (optional) reset before raytracing if rst = 1
        """

        cdef carma_context * context = carma_context.instance()
        cdef carma_obj[float] * d_screen = self.sensors.d_wfs[n].d_gs.d_phase.d_screen
        cdef carma_obj[float] * d_tel = tel.telescope.d_phase_ab_M1_m

        context.set_activeDeviceForce(self.sensors.device, 1)
        if(type_trace == "all"):
            self.sensors.d_wfs[n].sensor_trace(atmos.s_a, dms.dms)
            d_screen.axpy(1.0, d_tel, 1, 1)
        elif(type_trace == "atmos"):
            self.sensors.d_wfs[n].sensor_trace(atmos.s_a)
            d_screen.axpy(1.0, d_tel, 1, 1)
        elif(type_trace == "dm"):
            self.sensors.d_wfs[n].sensor_trace(dms.dms, rst)

    IF USE_MPI:
        cpdef Bcast_dscreen(self):
            """Broadcast the screen of every wfs on process 0 to all process

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
            """Broadcast the screen of every wfs on process 0 to all process

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
            """Gather the carma object 'bincube' of a wfs on the process 0

            :parameters:
                comm: (MPI.Intracomm) :communicator mpi

                :param n: (int) : number of the wfs where the gather will occured
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
            """Gather the carma object 'bincube' of a wfs on the process 0

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
        """Return the rank of one of the sensors wfs

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
        for i in range( < int > self.sensors.nsensors()):
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
        """Return the 'bincube' array of a given wfs

        :param n: (int) : number of the wfs to get the 'bincube' from
        """
        cdef carma_obj[int] * cube
        cdef const long * cdims
        cdef np.ndarray[ndim = 1, dtype = np.int32_t] data
        cube = self.sensors.d_wfs[n].d_hrmap
        cdims = cube.getDims()
        data = np.empty((cdims[1]), dtype=np.float32)
        cube.device2host( < int * > data.data)
        return data

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
