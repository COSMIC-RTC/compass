import numpy as np
import h5py
from Groot import groot_init
import time
import sys
import os

from guardian import gamora
from guardian.tools import Dphi
from guardian.tools import roket_exploitation as rexp
import matplotlib.pyplot as plt
plt.ion()


def compute_Cerr(filename, modal=True, ctype="float"):
    """ Returns the residual error covariance matrix using GROOT from a ROKET file
    :parameter:
        filename : (string) : full path to the ROKET file
        modal : (bool) : if True (default), Cerr is returned in the Btt modal basis,
                         in the actuator basis if False
        ctype : (string) : "float" or "double"
    :return:
        Cerr : (np.ndarray(dim=2, dtype=np.float32)) : residual error covariance matrix
    """
    f = h5py.File(filename, 'r')
    Lambda_tar = f.attrs["_Param_target__Lambda"][0]
    Lambda_wfs = f.attrs["_Param_wfs__Lambda"]
    dt = f.attrs["_Param_loop__ittime"]
    gain = f.attrs["_Param_controller__gain"]
    wxpos = f.attrs["_Param_wfs__xpos"][0]
    wypos = f.attrs["_Param_wfs__ypos"][0]
    r0 = f.attrs["_Param_atmos__r0"] * (Lambda_tar / Lambda_wfs)**(6. / 5.)
    RASC = 180. / np.pi * 3600.
    xpos = f["dm.xpos"][:]
    ypos = f["dm.ypos"][:]
    p2m = f.attrs["_Param_tel__diam"] / f.attrs["_Param_geom__pupdiam"]
    pupshape = int(2**np.ceil(np.log2(f.attrs["_Param_geom__pupdiam"]) + 1))
    xactu = (xpos - pupshape / 2) * p2m
    yactu = (ypos - pupshape / 2) * p2m
    H = f.attrs["_Param_atmos__alt"]
    L0 = f.attrs["_Param_atmos__L0"]
    speed = f.attrs["_Param_atmos__windspeed"]
    theta = f.attrs["_Param_atmos__winddir"] * np.pi / 180.
    frac = f.attrs["_Param_atmos__frac"]

    Htheta = np.linalg.norm([wxpos, wypos]) / RASC * H
    vdt = speed * dt / gain
    angleht = np.arctan2(wypos, wxpos)
    fc = 1 / (2 * (xactu[1] - xactu[0]))
    scale = (1 / r0)**(5 / 3.) * frac * (Lambda_tar / (2 * np.pi))**2
    Nact = f["Nact"][:]
    Nact = np.linalg.inv(Nact)
    P = f["P"][:]
    Btt = f["Btt"][:]
    Tf = Btt[:-2, :-2].dot(P[:-2, :-2])
    IF, T = rexp.get_IF(filename)
    IF = IF.T
    T = T.T
    N = IF.shape[0]
    deltaTT = T.T.dot(T) / N
    deltaF = IF.T.dot(T) / N
    pzt2tt = np.linalg.inv(deltaTT).dot(deltaF.T)

    if (ctype == "float"):
        groot = groot_init(Nact.shape[0],
                           int(f.attrs["_Param_atmos__nscreens"]), angleht, fc,
                           vdt.astype(np.float32),
                           Htheta.astype(np.float32), f.attrs["_Param_atmos__L0"], theta,
                           scale.astype(np.float32),
                           xactu.astype(np.float32),
                           yactu.astype(np.float32),
                           pzt2tt.astype(np.float32),
                           Tf.astype(np.float32), Nact.astype(np.float32))
    elif (ctype == "double"):
        groot = groot_initD(Nact.shape[0],
                            int(f.attrs["_Param_atmos__nscreens"]), angleht, fc,
                            vdt.astype(np.float64),
                            Htheta.astype(np.float64),
                            f.attrs["_Param_atmos__L0"].astype(np.float64),
                            theta.astype(np.float64),
                            scale.astype(np.float64),
                            xactu.astype(np.float64),
                            yactu.astype(np.float64),
                            pzt2tt.astype(np.float64),
                            Tf.astype(np.float64), Nact.astype(np.float64))
    else:
        raise TypeError("Unknown ctype : must be float or double")
    tic = time.time()
    groot.compute_Cerr()
    Cerr = groot.get_Cerr()
    cov_err_groot = np.zeros((Nact.shape[0] + 2, Nact.shape[0] + 2))
    cov_err_groot[:-2, :-2] = Cerr
    cov_err_groot[-2:, -2:] = groot.get_TTcomp()
    tac = time.time()
    print("Cee computed in : %.2f seconds" % (tac - tic))
    if (modal):
        cov_err_groot = P.dot(cov_err_groot).dot(P.T)

    f.close()
    return cov_err_groot


def compute_Cerr_cpu(filename, modal=True):
    """ Returns the residual error covariance matrix using CPU version of GROOT
    from a ROKET file
    :parameter:
        filename : (string) : full path to the ROKET file
        modal : (bool) : if True (default), Cerr is returned in the Btt modal basis,
                         in the actuator basis if False
    :return:
        Cerr : (np.ndarray(dim=2, dtype=np.float32)) : residual error covariance matrix
    """
    f = h5py.File(filename, 'r')

    tabx, taby = Dphi.tabulateIj0()
    Lambda_tar = f.attrs["_Param_target__Lambda"][0]
    Lambda_wfs = f.attrs["_Param_wfs__Lambda"]
    dt = f.attrs["_Param_loop__ittime"]
    gain = f.attrs["_Param_controller__gain"]
    wxpos = f.attrs["_Param_wfs__xpos"][0]
    wypos = f.attrs["_Param_wfs__ypos"][0]
    r0 = f.attrs["_Param_atmos__r0"] * (Lambda_tar / Lambda_wfs)**(6. / 5.)
    RASC = 180. / np.pi * 3600.
    xpos = f["dm.xpos"][:]
    ypos = f["dm.ypos"][:]
    p2m = f.attrs["_Param_tel__diam"] / f.attrs["_Param_geom__pupdiam"]
    pupshape = int(2**np.ceil(np.log2(f.attrs["_Param_geom__pupdiam"]) + 1))
    xactu = (xpos - pupshape / 2) * p2m
    yactu = (ypos - pupshape / 2) * p2m
    Ccov = np.zeros((xpos.size, xpos.size))
    Caniso = np.zeros((xpos.size, xpos.size))
    Cbp = np.zeros((xpos.size, xpos.size))
    xx = np.tile(xactu, (xactu.shape[0], 1))
    yy = np.tile(yactu, (yactu.shape[0], 1))
    xij = xx - xx.T
    yij = yy - yy.T

    for l in range(f.attrs["_Param_atmos__nscreens"]):
        H = f.attrs["_Param_atmos__alt"][l]
        L0 = f.attrs["_Param_atmos__L0"][l]
        speed = f.attrs["_Param_atmos__windspeed"][l]
        theta = f.attrs["_Param_atmos__winddir"][l] * np.pi / 180.
        frac = f.attrs["_Param_atmos__frac"][l]

        Htheta = np.linalg.norm([wxpos, wypos]) / RASC * H
        vdt = speed * dt / gain
        # Covariance matrices models on actuators space
        M = np.zeros((xpos.size, xpos.size))
        Mvdt = M.copy()
        Mht = M.copy()
        Mhvdt = M.copy()
        angleht = np.arctan2(wypos, wxpos)
        fc = xactu[1] - xactu[0]

        M = np.linalg.norm([xij, yij], axis=0)
        Mvdt = np.linalg.norm([xij - vdt * np.cos(theta), yij - vdt * np.sin(theta)],
                              axis=0)
        Mht = np.linalg.norm([
                xij - Htheta * np.cos(angleht), yij - Htheta * np.sin(angleht)
        ], axis=0)
        Mhvdt = np.linalg.norm([
                xij - vdt * np.cos(theta) - Htheta * np.cos(angleht),
                yij - vdt * np.sin(theta) - Htheta * np.sin(angleht)
        ], axis=0)

        Ccov += 0.5 * (Dphi.dphi_lowpass(Mhvdt, fc, L0, tabx, taby) - Dphi.dphi_lowpass(
                Mht, fc, L0, tabx, taby) - Dphi.dphi_lowpass(Mvdt, fc, L0, tabx, taby) +
                       Dphi.dphi_lowpass(M, fc, L0, tabx, taby)) * (1. / r0)**(
                               5. / 3.) * frac

        Caniso += 0.5 * (Dphi.dphi_lowpass(Mht, fc, L0, tabx, taby) - Dphi.dphi_lowpass(
                M, fc, L0, tabx, taby)) * (1. / r0)**(5. / 3.) * frac
        Cbp += 0.5 * (Dphi.dphi_lowpass(Mvdt, fc, L0, tabx, taby) - Dphi.dphi_lowpass(
                M, fc, L0, tabx, taby)) * (1. / r0)**(5. / 3.) * frac

    Sp = (Lambda_tar / (2 * np.pi))**2
    Ctt = (Caniso + Caniso.T) * Sp
    Ctt += ((Cbp + Cbp.T) * Sp)
    Ctt += ((Ccov + Ccov.T) * Sp)

    P = f["P"][:]
    Btt = f["Btt"][:]
    Tf = Btt[:-2, :-2].dot(P[:-2, :-2])

    IF, T = rexp.get_IF(filename)
    IF = IF.T
    T = T.T
    N = IF.shape[0]
    deltaTT = T.T.dot(T) / N
    deltaF = IF.T.dot(T) / N
    pzt2tt = np.linalg.inv(deltaTT).dot(deltaF.T)

    Nact = f["Nact"][:]
    N1 = np.linalg.inv(Nact)
    Ctt = N1.dot(Ctt).dot(N1)
    ttcomp = pzt2tt.dot(Ctt).dot(pzt2tt.T)
    Ctt = Tf.dot(Ctt).dot(Tf.T)
    cov_err = np.zeros((Ctt.shape[0] + 2, Ctt.shape[0] + 2))
    cov_err[:-2, :-2] = Ctt
    cov_err[-2:, -2:] = ttcomp
    if (modal):
        cov_err = P.dot(cov_err).dot(P.T)
    f.close()

    return cov_err


def compare_GPU_vs_CPU(filename):
    """ Compare results of GROOT vs its CPU version in terms of execution time
    and precision on the PSF renconstruction
    :parameter:
        filename : (string) : full path to the ROKET file

    """
    timer = ch.naga_timer()

    ch.threadSync()
    timer.start()
    ch.threadSync()
    synctime = timer.stop()
    timer.reset()

    timer.start()
    cov_err_gpu_s = compute_Cerr(filename)
    ch.threadSync()
    gpu_time_s = timer.stop() - synctime
    timer.reset()

    timer.start()
    cov_err_gpu_d = compute_Cerr(filename, ctype="double")
    ch.threadSync()
    gpu_time_d = timer.stop() - synctime
    timer.reset()

    tic = time.time()
    cov_err_cpu = compute_Cerr_cpu(filename)
    tac = time.time()
    cpu_time = tac - tic

    otftel, otf2, psf_cpu, gpu = gamora.psf_rec_Vii(filename, fitting=False,
                                                    cov=cov_err_cpu.astype(np.float32))
    otftel, otf2, psf_gpu_s, gpu = gamora.psf_rec_Vii(
            filename, fitting=False, cov=cov_err_gpu_s.astype(np.float32))
    otftel, otf2, psf_gpu_d, gpu = gamora.psf_rec_Vii(
            filename, fitting=False, cov=cov_err_gpu_d.astype(np.float32))

    print("-----------------------------------------")
    print("CPU time : ", cpu_time, " s ")
    print("GPU time simple precision : ", gpu_time_s, " s ")
    print("GPU time double precision : ", gpu_time_d, " s ")
    print("Max absolute difference in PSFs simple precision : ",
          np.abs(psf_cpu - psf_gpu_s).max())
    print("Max absolute difference in PSFs double precision : ",
          np.abs(psf_cpu - psf_gpu_d).max())
    gamora.cutsPSF(filename, psf_cpu, psf_gpu_s)
    gamora.cutsPSF(filename, psf_cpu, psf_gpu_d)


def compute_Ca_cpu(filename, modal=True):
    """ Returns the aliasing error covariance matrix using CPU version of GROOT
    from a ROKET file
    :parameter:
        filename : (string) : full path to the ROKET file
        modal : (bool) : if True (default), Ca is returned in the Btt modal basis,
                         in the actuator basis if False
    :return:
        Ca : (np.ndarray(dim=2, dtype=np.float32)) : aliasing error covariance matrix
    """
    f = h5py.File(filename, 'r')
    nsub = f["R"][:].shape[1] // 2
    nssp = f.attrs["_Param_wfs__nxsub"][0]
    validint = f.attrs["_Param_tel__cobs"]
    x = np.linspace(-1, 1, nssp)
    x, y = np.meshgrid(x, x)
    r = np.sqrt(x * x + y * y)
    rorder = np.sort(r.reshape(nssp * nssp))
    ncentral = nssp * nssp - np.sum(r >= validint)
    validext = rorder[ncentral + nsub]
    valid = (r < validext) & (r >= validint)
    ivalid = np.where(valid)
    xvalid = ivalid[0] + 1
    yvalid = ivalid[1] + 1
    ivalid = (xvalid, yvalid)
    d = f.attrs["_Param_tel__diam"] / (f.attrs["_Param_dm__nact"][0] - 1)
    r0 = f.attrs["_Param_atmos__r0"] * (f.attrs["_Param_target__Lambda"] / 0.5)**(
            6. / 5.)
    RASC = 180 / np.pi * 3600.

    scale = 0.23 * (d / r0)**(5 / 3.) * \
        (f.attrs["_Param_target__Lambda"] * 1e-6 / (2 * np.pi * d))**2 * RASC**2

    mask = np.zeros((nssp + 2, nssp + 2))
    Ca = np.identity(nsub * 2)

    for k in range(nsub):
        mask *= 0
        mask[xvalid[k], yvalid[k]] = 1
        mask[xvalid[k], yvalid[k] - 1] = -0.5
        mask[xvalid[k], yvalid[k] + 1] = -0.5
        Ca[k, :nsub] = mask[ivalid].flatten()
        mask *= 0
        mask[xvalid[k], yvalid[k]] = 1
        mask[xvalid[k] - 1, yvalid[k]] = -0.5
        mask[xvalid[k] + 1, yvalid[k]] = -0.5
        Ca[k + nsub, nsub:] = mask[ivalid].flatten()

    R = f["R"][:]
    Ca = R.dot(Ca * scale).dot(R.T)
    if (modal):
        P = f["P"][:]
        Ca = P.dot(Ca).dot(P.T)
    f.close()
    return Ca


def compute_Cn_cpu(filename, model="data", modal=True):
    """ Returns the noise error covariance matrix using CPU version of GROOT
    from a ROKET file
    :parameter:
        filename : (string) : full path to the ROKET file
        modal : (bool) : if True (default), Cn is returned in the Btt modal basis,
                         in the actuator basis if False
    :return:
        Cn : (np.ndarray(dim=2, dtype=np.float32)) : noise error covariance matrix
    """
    f = h5py.File(filename, 'r')
    if (model == "data"):
        N = f["noise"][:]
        Cn = N.dot(N.T) / N.shape[1]
        if modal:
            P = f["P"][:]
            Cn = P.dot(Cn).dot(P.T)
    else:
        nslopes = f["R"][:].shape[1]
        Cn = np.zeros(nslopes)
        noise = f.attrs["_Param_wfs__noise"][0]
        RASC = 180 / np.pi * 3600.
        if (noise >= 0):
            Nph = f.attrs["_Param_wfs__zerop"] * 10 ** (-0.4 * f.attrs["_Param_wfs__gsmag"]) * \
                f.attrs["_Param_wfs__optthroughput"] * \
                (f.attrs["_Param_tel__diam"] / f.attrs["_Param_wfs__nxsub"]
                 ) ** 2. * f.attrs["_Param_loop__ittime"]

            r0 = (f.attrs["_Param_wfs__Lambda"] / 0.5)**(
                    6.0 / 5.0) * f.attrs["_Param_atmos__r0"]

            sig = (np.pi ** 2 / 2) * (1 / Nph) * \
                (1. / r0) ** 2  # Photon noise in m^-2
            # Noise variance in arcsec^2
            sig = sig * ((f.attrs["_Param_wfs__Lambda"] * 1e-6) /
                         (2 * np.pi))**2 * RASC**2

            Ns = f.attrs["_Param_wfs__npix"]  # Number of pixel
            Nd = (f.attrs["_Param_wfs__Lambda"] * 1e-6
                  ) * RASC / f.attrs["_Param_wfs__pixsize"]
            sigphi = (np.pi ** 2 / 3.0) * (1 / Nph ** 2) * (f.attrs["_Param_wfs__noise"]) ** 2 * \
                Ns ** 2 * (Ns / Nd) ** 2  # Phase variance in m^-2
            # Noise variance in arcsec^2
            sigsh = sigphi * \
                ((f.attrs["_Param_wfs__Lambda"] * 1e-6) / (2 * np.pi)) ** 2 * RASC ** 2

            Cn[:len(sig)] = sig + sigsh
            Cn[len(sig):] = sig + sigsh

        Cn = np.diag(Cn)
        R = f["R"][:]
        Cn = R.dot(Cn).dot(R.T)
        if (modal):
            P = f["P"][:]
            Cn = P.dot(Cn).dot(P.T)
    f.close()
    return Cn


def compute_OTF_fitting(filename, otftel):
    """
    Modelize the OTF due to the fitting using dphi_highpass

    :parameters:
        filename: (str) : ROKET hdf5 file path
        otftel: (np.ndarray) : Telescope OTF
    :return:
        otf_fit: (np.ndarray) : Fitting OTF
        psf_fit (np.ndarray) : Fitting PSF
    """
    f = h5py.File(filename, 'r')
    r0 = f.attrs["_Param_atmos__r0"] * (f.attrs["_Param_target__Lambda"][0] / 0.5)**(
            6. / 5.)
    ratio_lambda = 2 * np.pi / f.attrs["_Param_target__Lambda"][0]
    # Telescope OTF
    spup = rexp.get_pup(filename)
    mradix = 2
    fft_size = mradix**int((np.log(2 * spup.shape[0]) / np.log(mradix)) + 1)
    mask = np.ones((fft_size, fft_size))
    mask[np.where(otftel < 1e-5)] = 0

    x = np.arange(fft_size) - fft_size / 2
    pixsize = f.attrs["_Param_tel__diam"] / f.attrs["_Param_geom__pupdiam"]
    x = x * pixsize
    r = np.sqrt(x[:, None] * x[:, None] + x[None, :] * x[None, :])
    tabx, taby = Dphi.tabulateIj0()
    dphi = np.fft.fftshift(
            Dphi.dphi_highpass(r, f.attrs["_Param_tel__diam"] /
                               (f.attrs["_Param_dm__nact"][0] - 1), tabx, taby) *
            (1 / r0)**(5 / 3.))  # * den * ratio_lambda**2 * mask
    otf_fit = np.exp(-0.5 * dphi) * mask
    otf_fit = otf_fit / otf_fit.max()

    psf_fit = np.fft.fftshift(np.real(np.fft.ifft2(otftel * otf_fit)))
    psf_fit *= (fft_size * fft_size / float(np.where(spup)[0].shape[0]))

    f.close()
    return otf_fit, psf_fit


def compute_PSF(filename):
    """
    Modelize the PSF using GROOT model for aniso and bandwidth, Gendron model for aliasing,
    dphi_highpass for fitting, noise extracted from datas. Non linearity not taken into account

    :parameters:
        filename: (str) : ROKET hdf5 file path
    :return:
        psf: (np.ndarray) : PSF
    """
    tic = time.time()
    spup = rexp.get_pup(filename)
    Cab = compute_Cerr(filename)
    Cn = compute_Cn_cpu(filename)
    Ca = compute_Ca_cpu(filename)
    Cee = Cab + Cn + Ca
    otftel, otf2, psf, gpu = gamora.psf_rec_Vii(filename, fitting=False,
                                                cov=(Cee).astype(np.float32))
    otf_fit, psf_fit = compute_OTF_fitting(filename, otftel)
    psf = np.fft.fftshift(np.real(np.fft.ifft2(otf_fit * otf2 * otftel)))
    psf *= (psf.shape[0] * psf.shape[0] / float(np.where(spup)[0].shape[0]))
    tac = time.time()
    print("PSF computed in ", tac - tic, " seconds")

    return psf


def test_Calias(filename):
    f = h5py.File(filename, 'r')
    tabx, taby = Dphi.tabulateIj0()
    nsub = f["Cmm"][:].shape[0] // 2
    nssp = f.attrs["_Param_wfs__nxsub"][0]
    npix = f.attrs["_Param_wfs__npix"][0]
    validint = f.attrs["_Param_tel__cobs"]
    x = np.linspace(-1, 1, nssp)
    x, y = np.meshgrid(x, x)
    r = np.sqrt(x * x + y * y)
    rorder = np.sort(r.reshape(nssp * nssp))
    ncentral = nssp * nssp - np.sum(r >= validint, dtype=np.int32)
    validext = rorder[ncentral + nsub]
    valid = (r < validext) & (r >= validint)
    ivalid = np.where(valid)
    r0 = f.attrs["_Param_atmos__r0"]
    Lambda_wfs = f.attrs["_Param_wfs__Lambda"][0]
    d = f.attrs["_Param_tel__diam"] / nssp
    RASC = 180 / np.pi * 3600
    scale = (RASC * Lambda_wfs * 1e-6 / 2 / np.pi)**2 / d**2
    x = (np.arange(nssp) - nssp / 2) * d
    x, y = np.meshgrid(x, x)
    x = x[ivalid]
    y = y[ivalid]
    fc = d  #/ npix
    xx = np.tile(x, (nsub, 1))
    yy = np.tile(y, (nsub, 1))

    xx = xx - xx.T
    yy = yy - yy.T

    AB = np.linalg.norm([xx, yy], axis=0)
    Ab = np.linalg.norm([xx - d, yy], axis=0)
    aB = np.linalg.norm([xx + d, yy], axis=0)
    ab = AB

    Cmm = np.zeros((2 * nsub, 2 * nsub))
    Cmm[:nsub, :nsub] = 0.5 * (
            Dphi.dphi_highpass(Ab, fc, tabx, taby) + Dphi.dphi_highpass(
                    aB, fc, tabx, taby) - 2 * Dphi.dphi_highpass(AB, fc, tabx, taby)) * (
                            1 / r0)**(5. / 3.)
    CD = AB
    Cd = np.linalg.norm([xx, yy - d], axis=0)
    cD = np.linalg.norm([xx, yy + d], axis=0)
    cd = CD

    Cmm[nsub:, nsub:] = 0.5 * (
            Dphi.dphi_highpass(Cd, fc, tabx, taby) + Dphi.dphi_highpass(
                    cD, fc, tabx, taby) - 2 * Dphi.dphi_highpass(CD, fc, tabx, taby)) * (
                            1 / r0)**(5. / 3.)

    # aD = np.linalg.norm([xx + d/2, yy + d/2], axis=0)
    # ad = np.linalg.norm([xx + d/2, yy - d/2], axis=0)
    # Ad = np.linalg.norm([xx - d/2, yy - d/2], axis=0)
    # AD = np.linalg.norm([xx - d/2, yy + d/2], axis=0)
    #
    # Cmm[nsub:,:nsub] = 0.25 * (Dphi.dphi_highpass(Ad, d, tabx, taby)
    #                 + Dphi.dphi_highpass(aD, d, tabx, taby)
    #                 - Dphi.dphi_highpass(AD, d, tabx, taby)
    #                 - Dphi.dphi_highpass(ad, d, tabx, taby)) * (1 / r0)**(5. / 3.)
    # Cmm[:nsub,nsub:] = Cmm[nsub:,:nsub].copy()

    a = f["alias_meas"][:]
    Calias = a.dot(a.T) / a.shape[1]

    return Cmm * scale, Calias
