#!/usr/local/bin/python3.6
# encoding: utf-8
'''
Created on 3 aout 2017

@author: fferreira
'''
import numpy as np
from shesha_config import shesha_constants as scons
from shesha_config.PDMS import Param_dm
from shesha_config.PGEOM import Param_geom
from . import utilities as util
from Dms import Dms
from scipy.sparse import csr_matrix


def dim_dm_support(
        cent: float,
        extent: int,
        ssize: int):
    """
    Compute the DM support dimensions

    :parameters:
        cent : (float): center of the pupil
        extent: (float): size of the DM support
        ssize: (int): size of ipupil support
    """
    n1 = np.floor(cent - extent / 2)
    n2 = np.ceil(cent + extent / 2)
    if(n1 < 1):
        n1 = 1
    if(n2 > ssize):
        n2 = ssize

    return int(n1), int(n2)


def dim_dm_patch(
        pupdiam: int,
        diam: float,
        type_dm: bytes,
        alt: float,
        xpos_wfs: list,
        ypos_wfs: list):
    """ compute patchDiam for DM

    :parameters:
        pupdiam: (int) : pupil diameter

        diam: (float) : telescope diameter

        type_dm: (bytes) : type of dm

        alt: (float) : altitude of dm

        xpos_wfs: (list) : list of wfs xpos

        ypos_wfs: (list) : list of wfs ypos
    """

    norms = [
        np.linalg.norm([xpos_wfs[w], ypos_wfs[w]])
        for w in range(len(xpos_wfs))]
    if ((type_dm == scons.DmType.PZT) or (type_dm == scons.DmType.TT)):
        pp = (diam * pupdiam)
    elif (type_dm == scons.DmType.KL):
        pp = (pupdiam)
    else:
        raise TypeError("This type of DM doesn't exist ")

    patchDiam = int(
        pupdiam + 2 * np.max(norms) * scons.ARCSEC2RAD * np.abs(alt) /
        (pp))
    return patchDiam


def createSquarePattern(
        pitch: float,
        nxact: int):
    """
    Creates a list of M=nxact^2 actuator positions spread over an square grid.
    Coordinates are centred around (0,0).

    :parameters:
        pitch: (float) : distance in pixels between 2 adjacent actus
        nxact: (int) : number of actu across the pupil diameter
    :return:
        xy: (np.ndarray(dims=2,dtype=np.float32)) : xy[M,2] list of coodinates
    """

    xy = np.tile(np.arange(nxact) - (nxact - 1.) /
                 2., (nxact, 1)).astype(np.float32)
    xy = np.array([xy.flatten(), xy.T.flatten()]) * pitch
    xy = np.float32(xy)
    return xy


def createHexaPattern(
        pitch: float,
        supportSize: int):
    """
    Creates a list of M actuator positions spread over an hexagonal grid.
    The number M is the number of points of this grid, it cannot be
    known before the procedure is called.
    Coordinates are centred around (0,0).
    The support that limits the grid is a square [-n/2,n/2].

    :parameters:
        pitch: (float) : distance in pixels between 2 adjacent actus
        n: (float) : size in pixels of the support over which the coordinate list
             should be returned.
    :return:
        xy: (np.ndarray(dims=2,dtype=np.float32)) : xy[M,2] list of coodinates
    """
    V3 = np.sqrt(3)
    nx = int(np.ceil((supportSize / 2.0) / pitch) + 1)
    x = pitch * (np.arange(2 * nx + 1, dtype=np.float32) - nx)
    Nx = x.shape[0]
    ny = int(np.ceil((supportSize / 2.0) / pitch / V3) + 1)
    y = (V3 * pitch) * (np.arange(2 * ny + 1, dtype=np.float32) - ny)
    Ny = y.shape[0]
    x = np.tile(x, (Ny, 1)).flatten()
    y = np.tile(y, (Nx, 1)).T.flatten()
    x = np.append(x, x + pitch / 2.)
    y = np.append(y, y + pitch * V3 / 2.)
    xy = np.float32(np.array([y, x]))
    return xy


def createDoubleHexaPattern(
        pitch: float,
        supportSize: int):
    """
    Creates a list of M actuator positions spread over an hexagonal grid.
    The number M is the number of points of this grid, it cannot be
    known before the procedure is called.
    Coordinates are centred around (0,0).
    The support that limits the grid is a square [-n/2,n/2].

    :parameters:
        pitch: (float) : distance in pixels between 2 adjacent actus
        n: (float) : size in pixels of the support over which the coordinate list
             should be returned.
    :return:
        xy: (np.ndarray(dims=2,dtype=np.float32)) : xy[M,2] list of coodinates
    """
    V3 = np.sqrt(3)
    pi = np.pi
    nx = int(np.ceil((supportSize / 2.0) / pitch) + 1)
    x = pitch * (np.arange(2 * nx + 1, dtype=np.float32) - nx)
    Nx = x.shape[0]
    ny = int(np.ceil((supportSize / 2.0) / pitch / V3) + 1)
    y = (V3 * pitch) * (np.arange(2 * ny + 1, dtype=np.float32) - ny) + pitch
    Ny = y.shape[0]
    x = np.tile(x, (Ny, 1)).flatten()
    y = np.tile(y, (Nx, 1)).T.flatten()
    x = np.append(x, x + pitch / 2.)
    y = np.append(y, y + pitch * V3 / 2.)
    xy = np.float32(np.array([x, y]))

    th = np.arctan2(y, x)
    nn = np.where(((th > pi / 3) & (th < 2 * pi / 3)))
    x = x[nn]
    y = y[nn]
    X = np.array([])
    Y = np.array([])
    for k in range(6):
        xx = np.cos(k * pi / 3) * x + np.sin(k * pi / 3) * y
        yy = -np.sin(k * pi / 3) * x + np.cos(k * pi / 3) * y
        X = np.r_[X, xx]
        Y = np.r_[Y, yy]
    return np.float32(np.array([Y, X]))


def select_actuators(
        xc: np.ndarray,
        yc: np.ndarray,
        nxact: int,
        pitch: int,
        cobs: float,
        margin_in: float,
        margin_out: float,
        N=None):
    """
    Select the "valid" actuators according to the system geometry
    :parameters:
        p_dm: (Param_dm) : dm settings
        cobs: telescope cobs
        xc: actuators x positions (origine in center of mirror)
        yc: actuators y positions (origine in center of mirror)

    :return:
        liste_fin: actuator indice selection for xpos/ypos


    """
    # the following determine if an actuator is to be considered or not
    # relative to the pitchmargin parameter.
    dis = np.sqrt(xc**2 + yc**2)

    # test Margin_in
    rad_in = (((nxact - 1) / 2) * cobs - margin_in) * pitch

    if N is None:
        if(margin_out < 0):
            margin_out = 1.44
        rad_out = ((nxact - 1.) / 2. + margin_out) * pitch

        valid_actus = np.where((dis <= rad_out) * (dis >= rad_in))[0]

    else:
        valid_actus = np.where(dis >= rad_in)[0]
        indsort = np.argsort(dis[valid_actus])

        if(N > valid_actus.size):
            print('Too many actuators wanted, restricted to ', valid_actus.size)
        else:
            valid_actus = np.sort(indsort[:N])

    return valid_actus


def make_zernike(
        nzer: int,
        size: int,
        diameter: int,
        xc=-1.,
        yc=-1.,
        ext=0):
    """Compute the zernike modes

    :parameters:
        nzer: (int) : number of modes

        size: (int) : size of the screen

        diameter: (int) : pupil diameter

        xc: (float) : (optional) x-position of the center


        yc: (float) : (optional) y-position of the center

        ext: (int) : (optional) extension
    :return:
        z : (np.ndarray(ndims=3,dtype=np.float64)) : zernikes modes
    """
    m = 0
    n = 0

    if(xc == -1):
        xc = size / 2
    if(yc == -1):
        yc = size / 2

    radius = (diameter + 1.) / 2.
    zr = util.dist(size, xc, yc).astype(np.float32).T / radius
    zmask = np.zeros((zr.shape[0], zr.shape[1], nzer), dtype=np.float32)
    zmaskmod = np.zeros((zr.shape[0], zr.shape[1], nzer), dtype=np.float32)

    zmask[:, :, 0] = (zr <= 1).astype(np.float32)
    zmaskmod[:, :, 0] = (zr <= 1.2).astype(np.float32)

    for i in range(1, nzer):
        zmask[:, :, i] = zmask[:, :, 0]
        zmaskmod[:, :, i] = zmaskmod[:, :, 0]

    zrmod = zr * zmaskmod[:, :, 0]

    zr = zr * zmask[:, :, 0]

    x = np.tile(np.linspace(1, size, size).astype(np.float32), (size, 1))
    zteta = np.arctan2(x - yc, x.T - xc).astype(np.float32)

    z = np.zeros((size, size, nzer), dtype=np.float32)

    for zn in range(nzer):
        n, m = zernumero(zn + 1)

        if ext:
            for i in range((n - m) // 2 + 1):
                z[:, :, zn] = z[:, :, zn] + (-1.) ** i * zrmod ** (n - 2. * i) * float(np.math.factorial(n - i)) / \
                    float(np.math.factorial(i) * np.math.factorial((n + m) / 2 - i) *
                          np.math.factorial((n - m) / 2 - i))
        else:
            for i in range((n - m) // 2 + 1):
                z[:, :, zn] = z[:, :, zn] + (-1.) ** i * zr ** (n - 2. * i) * float(np.math.factorial(n - i)) / \
                    float(np.math.factorial(i) * np.math.factorial((n + m) / 2 - i) *
                          np.math.factorial((n - m) / 2 - i))

        if((zn + 1) % 2 == 1):
            if(m == 0):
                z[:, :, zn] = z[:, :, zn] * np.sqrt(n + 1.)
            else:
                z[:, :, zn] = z[:, :, zn] * \
                    np.sqrt(2. * (n + 1)) * np.sin(m * zteta)
        else:
            if(m == 0):
                z[:, :, zn] = z[:, :, zn] * np.sqrt(n + 1.)
            else:
                z[:, :, zn] = z[:, :, zn] * \
                    np.sqrt(2. * (n + 1)) * np.cos(m * zteta)

    if(ext):
        return z * zmaskmod
    else:
        return z * zmask


def zernumero(zn: int):
    """
    Returns the radial degree and the azimuthal number of zernike
    number zn, according to Noll numbering (Noll, JOSA, 1976)

    :parameters:
        zn: (int) : zernike number

    :returns:
        rd: (int) : radial degrees

        an: (int) : azimuthal numbers

    """
    j = 0
    for n in range(101):
        for m in range(n + 1):
            if((n - m) % 2 == 0):
                j = j + 1
                if(j == zn):
                    return n, m
                if(m != 0):
                    j = j + 1
                    if(j == zn):
                        return n, m


def compute_KLbasis(
        g_dm: Dms,
        p_dm: Param_dm,
        p_geom: Param_geom,
        r0: float,
        diam: float):
    """Compute a Karhunen-Loeve basis for the dm:
            - compute the phase covariance matrix on the actuators using Kolmogorov
            - compute the geometric covariance matrix
            - double diagonalisation to obtain KL basis

    :parameters:
        g_dm: (Dms) : Dms object

        p_dm: (Param_dm) : dm settings

        p_geom: (Param_geom) : geom settings

        r0: (float) : atmos r0 in meter

        diam: (float) : telescope diameter
    """

    if(p_dm.type_dm == scons.DmType.PZT):
        tmp = (p_geom._ipupil.shape[0] - (p_dm._n2 - p_dm._n1 + 1)) // 2
        tmp_e0 = p_geom._ipupil.shape[0] - tmp
        tmp_e1 = p_geom._ipupil.shape[1] - tmp
        pup = p_geom._ipupil[tmp:tmp_e0, tmp:tmp_e1]
        indx_valid = np.where(pup.flatten("F") > 0)[0].astype(np.int32)
        p2m = diam / p_geom.pupdiam
        norm = -(p2m * diam / (2 * r0)) ** (5. / 3)

        g_dm.compute_KLbasis(scons.DmType.PZT, p_dm.alt, p_dm._xpos,
                             p_dm._ypos, indx_valid, indx_valid.size, norm, 1.0)
        KLbasis = np.fliplr(g_dm.get_KLbasis(scons.DmType.PZT, p_dm.alt))
    else:
        raise TypeError("DM must be pzt type")

    return KLbasis


def compute_DMbasis(
        g_dm: Dms,
        p_dm: Param_dm,
        p_geom: Param_geom):
    """Compute a the DM basis as a sparse matrix :
            - push on each actuator
            - get the corresponding dm shape
            - apply pupil mask and store in a column

    :parameters:
        g_dm: (Dms) : Dms object

        p_dm: (Param_dm) : dm settings

        p_geom: (Param_geom) : geom settings
    :return:
        IFbasis = (csr_matrix) : DM IF basis
    """
    tmp = (p_geom._ipupil.shape[0] - (p_dm._n2 - p_dm._n1 + 1)) // 2
    tmp_e0 = p_geom._ipupil.shape[0] - tmp
    tmp_e1 = p_geom._ipupil.shape[1] - tmp
    pup = p_geom._ipupil[tmp:tmp_e0, tmp:tmp_e1]
    indx_valid = np.where(pup.flatten("F") > 0)[0].astype(np.int32)

    #IFbasis = np.ndarray((indx_valid.size, p_dm._ntotact), dtype=np.float32)
    for i in range(p_dm._ntotact):
        g_dm.resetdm(p_dm.type_dm, p_dm.alt)
        g_dm.comp_oneactu(p_dm.type_dm, p_dm.alt, i, 1.0)
        shape = g_dm.get_dm(p_dm.type_dm, p_dm.alt)
        IFvec = csr_matrix(shape.flatten("F")[indx_valid])
        if(i == 0):
            val = IFvec.data
            col = IFvec.indices
            row = np.append(0, IFvec.getnnz())
        else:
            val = np.append(val, IFvec.data)
            col = np.append(col, IFvec.indices)
            row = np.append(row, row[-1] + IFvec.getnnz())
    g_dm.resetdm(p_dm.type_dm, p_dm.alt)
    IFbasis = csr_matrix((val, col, row))
    return IFbasis


def compute_IFsparse(
        g_dm: Dms,
        p_dms: list,
        p_geom: Param_geom):
    """Compute the influence functions of all DMs as a sparse matrix :
            - push on each actuator
            - get the corresponding dm shape
            - apply pupil mask and store in a column

    :parameters:
        g_dm: (Dms) : Dms object

        p_dms: (Param_dms) : dms settings

        p_geom: (Param_geom) : geom settings
    :return:
        IFbasis = (csr_matrix) : DM IF basis
    """
    ndm = len(p_dms)
    for i in range(ndm):
        IFi = computeDMbasis(g_dm, p_dms[i], p_geom)
        if(i == 0):
            val = IFi.data
            col = IFi.indices
            row = IFi.indptr
        else:
            val = np.append(val, IFi.data)
            col = np.append(col, IFi.indices)
            row = np.append(row, row[-1] + IFi.indptr[1:])
    IFsparse = csr_matrix((val, col, row))
    return IFsparse
