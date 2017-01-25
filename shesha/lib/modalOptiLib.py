'''
Created on 23 janv. 2017

@author: vdeo
'''

import numpy as np
import time
import matplotlib.pyplot as plt



def applyVoltGetSlopes(wao):
    wao.rtc.applycontrol(0, wao.dms)
    for w in range(len(wao.config.p_wfss)):
        wao.wfs.reset_phase(w)
        wao.wfs.sensors_trace(w, "dm", wao.tel, wao.atm, wao.dms)
        wao.wfs.sensors_compimg(w)
    wao.rtc.docentroids(0)
    return wao.rtc.getCentroids(0)


def measureIMatKL(wao, ampliVec, KL2V, Nslopes):
    iMatKL = np.zeros((KL2V.shape[1], Nslopes))
    wao.openLoop()
    st = time.time()
    for kl in range(KL2V.shape[1]):
        v = ampliVec[kl] * KL2V[:, kl:kl + 1].T.copy()
        wao.rtc.set_perturbcom(0, v)
        iMatKL[kl, :] = applyVoltGetSlopes(wao)
        print "Doing KL interaction matrix on mode: #%d\r" % kl,
    print "Modal interaction matrix done in %3.1f seconds" % (time.time() - st)
    return iMatKL


def normalizeKL2V(KL2V):
    KL2VNorm = KL2V * 0.
    for kl in range(KL2V.shape[1]):
        klmaxVal = np.abs(KL2V[:, kl]).max()
        # Norm du max des modes a 1
        KL2VNorm[:, kl] = KL2V[:, kl] / klmaxVal
    return KL2VNorm


def cropKL2V(KL2V):
    nkl = KL2V.shape[1]
    KL2Vmod = np.zeros((KL2V.shape[0], KL2V.shape[1] - 2))
    KL2Vmod[:, 0:-2] = KL2V[:, 0:-4]
    KL2Vmod[:, -2:] = KL2V[:, -2:]
    return KL2Vmod


def plotVDm(wao, numdm, V, size = 16, fignum = False):
    """
    plotVDm(wao, 0, KL2V[:,0][:-2], size=25, fignum=15)
    """
    if(wao.config.p_dms[numdm]._j1.shape[0] != V.shape[0]):
        print "Error V should have %d dimension not %d " % (wao.config.p_dms[numdm]._j1.shape[0], V.shape[0])
    else:
        if(fignum):
            plt.figure(fignum)
        plt.clf()
        plt.scatter(wao.config.p_dms[numdm]._j1, wao.config.p_dms[
                    numdm]._i1, c = V, marker = "s", s = size ** 2)
        plt.colorbar()


def computeKLModesImat(wao, pushDMMic, pushTTArcsec, KL2V, Nslopes):
    modesAmpli = np.ones(KL2V.shape[0])
    # normalisation de la valeur max des KL:
    # Amplitude des modes DM en microns et arc sec pour le TT
    KL2VNorm = normalizeKL2V(KL2V)
    NKL = KL2VNorm.shape[1]
    modesAmpli[0:NKL - 2] = pushDMMic  # DM en microns DM en microns
    modesAmpli[NKL - 2:] = pushTTArcsec  # arc sec pour le TT
    imat = measureIMatKL(wao, modesAmpli, KL2VNorm, Nslopes)
    imat[:-2, :] /= pushDMMic
    imat[-2:, :] /= pushTTArcsec
    return imat, KL2VNorm


def computeCmatKL(DKL, KL2V, nfilt, gains, gainTT):
    nmodes = (DKL.shape[0] - nfilt)
    KL2V2Filt = np.zeros((KL2V.shape[0], KL2V.shape[1] - nfilt))
    DKLfilt = np.zeros((nmodes, DKL.shape[1]))
    # Filter the nfilt modes
    DKLfilt[0:nmodes, :] = DKL[0:nmodes, :]
    DKLfilt[-2:, :] = DKL[-2:, :]  # Concatenating the TT (last 2 columns)
    KL2V2Filt[:, 0:nmodes] = KL2V[:, 0:nmodes]
    KL2V2Filt[:, -2:] = KL2V[:, -2:]
    # Direct inversion
    Dmp = np.linalg.inv(DKLfilt.dot(DKLfilt.T)).dot(DKLfilt)
    S2M = Dmp
    for i in range(nmodes - 2):
        Dmp[i, :] *= gains[i]
    Dmp[-2:, :] *= gainTT
    # Command matrix
    cmat = KL2V2Filt.dot(Dmp).astype(np.float32)
    return (cmat, S2M)  # Dmp : slopes -> modes

"""
ENTRY POINT
"""

def setupStandardCmat(wao):

    Nactu = sum(wao.config.p_rtc.controllers[0].nactu)
    Nslopes = wao.rtc.getCentroids(0).shape[0]
    wao.setIntegratorLaw()
    gain = 0.7  # FIXME

    wao.openLoop()  # openLoop

    # wao.rtc.set_gain(0, gain) # Useless for generic controller

    decay = np.ones(Nslopes, dtype = (np.float32))
    wao.rtc.set_decayFactor(0, decay)

    mgain = gain * np.ones(Nactu, dtype = (np.float32))
    wao.rtc.set_mgain(0, mgain)

    wao.closeLoop()  # closeLoop

    KL2V = wao.returnkl2V()

    pushDMMic = 0.01  # 10nm
    pushTTArcsec = 0.005  # 5 mas
    wao.rtc.do_centroids_ref(0)
    miKL, KL2VN = computeKLModesImat(wao, pushDMMic, pushTTArcsec, KL2V, Nslopes)
    nfilt = 0  # 450

#     magicModalGain = np.linspace(1, 0.3, Nactu - 2 - nfilt)
    magicModalGain = np.ones((Nactu - 2 - nfilt))
    (cmat, S2KL) = computeCmatKL(miKL, KL2VN, nfilt, magicModalGain, 1)

    wao.rtc.set_cmat(0, cmat.astype(np.float32).copy())

    gDM = 1.
    gTT = 1.
    gains = np.ones(Nactu, dtype = np.float32) * gDM
    gains[-2:] = gTT
    wao.rtc.set_mgain(0, gains)
    wao.closeLoop()

    return miKL.T, S2KL, KL2VN


def updateToOptimalCmat(wao, slopes, volts, miKL, S2KL, M2V):

    nfilt = 0

    mia = np.dot(miKL, np.linalg.pinv(M2V))

    gainVector = modalControlOptimizationClosedLoopData(
            slopes, volts, mia, S2KL, M2V)


    optiCmat = computeSystemOptimizedControlMatrix(S2KL, gainVector, M2V)

    wao.openLoop()
    wao.rtc.set_cmat(0, optiCmat.astype(np.float32).copy())

    wao.closeLoop()

    return gainVector

def main(wao, nFrames):
    (miKL, S2KL, KL2VN) = setupStandardCmat(wao)
    
    nSlopes = S2KL.shape[1]
    nActu = KL2VN.shape[0]
    
    wao.closeLoop()
    wao.set_atmos(True)

    wao.aoLoopClicked(False)

    slopes = np.zeros((nSlopes, nFrames + 1), dtype=np.float32)
    volts = np.zeros((nActu, nFrames + 1), dtype=np.float32)

    for i in range(nFrames + 1):
        wao.mainLoop()
        slopes[:, i] = wao.getSlopes()
        volts[:, i] = wao.getVolts()

    gainVector = updateToOptimalCmat(wao, slopes, volts, miKL, S2KL, KL2VN)
        
    return gainVector



def hcor(freq, Fs, latency, G, BP):
    """
    Input arguments :
    <freq> is a 1D array of frequencies (usually 1024 points ranging from Fs/2048 to Fs/2).
    <Fs> is the sampling frequency
    <latency> is the latency, in seconds, between the BEGINNING of the integration and the
              start of the command.
    <G> is a scalar, it's the loop gain
    <BP> is a scalar, the cutoff frequency of the DM (seen as a 1st order filter)

    On output, returns the square modulus of the correction transfer function of the
    system
    """
    Te = 1. / Fs  # sampling period
    p = 1j * 2 * np.pi * freq  # Laplace variable
    Hint = 1. / (1 - np.exp(-p * Te))  # numeric integrator
    Hccd = (1. - np.exp(-p * Te)) / (p * Te)  # zero-order hold with 1/2 frame delay
    Hdac = Hccd  # well, same.
    tdelay = latency - Te  # time between END of the integratino and start of command
    Hret = np.exp(-p * tdelay)  # latency transfer function
    # transfer func of the DM, as a 1st order filter
    Hmir = 1. / (1. + 1j * freq / BP)
    # open-loop transfer function
    Hbo = Hint * Hccd * Hdac * Hret * Hmir
    # correction transfer function
    Hcor = 1. / np.abs(1 + Hbo * G) ** 2
    return Hcor







def modalControlOptimizationOpenLoopData(slopes, S2M, M2V, gmax = 1, Fs = 500, latency = 1, BP = 1e12):
    """
    Input arguments :
    <slopes>     : set of N frames of slopes (i.e. centroids), (nSlopes x N)
    <S2M>        : modal control matrix
    <M2V>        : matrix of system modes
    <gmax>       : scalar floating point value, maximum gain.
    <Fs>         : sampling frequency (Hertz)
    <latency>    : latency, in seconds, i.e. time between the start of WFS integration and time
                  when the actuator reaches 50% of its command.
    <BP>         : scalar, cutoff frequency of the DM (seen as a 1st order filter)

    Returned value :
    Array of optimal gains on the nModes modes.

    """
    # determines how many frames are recorded in the data set,
    # will be used later (could also be passed as an argument..)
    (_, nrec) = slopes.shape

    # conversion of slopes into modal coefficients
    modes = np.dot(S2M, slopes)
    (nmod, _) = modes.shape  # just to have the number of modes here

    # Produces an array of gains
    ngain = 100  # value to be adjusted during AITs
    gmin = 0.01  # TBC, maybe we'll keep gmin=0.001 to allow for static aberration compensation, if any. TBC during AIT.
    G = np.linspace(np.sqrt(gmin), np.sqrt(gmax), ngain) ** 2  # 1D array from gmin to gmax in ngain points

    npfft = nrec // 2  # number of useful points in the FFT
    # create a 1D array of frequency ranging from Fs/nrec to Fs/2.0 in npfft points
    freq = np.linspace(Fs / nrec, Fs / 2.0, npfft)
    optimumGain = np.zeros(nmod)  # memory alloc for 1D array of optimum gains

    # Initialization of all transfer functions computed for all the gains
    Hcor = np.zeros((ngain, npfft))
    for j in range(ngain):
        Hcor[j, :] = hcor(freq, Fs, latency, G[j], BP)

    # Fourier transform of modes coefficients.
    # Compute the square modulus of Fourier transform along time, for each mode.
    # Then multiplies by transfer function for different gains, sum everything
    # over frequency to get the error, and selects the minimum value.
    for i in range(nmod):
        # square modulus of Fast Fourier Transform
        tmp = np.abs(np.fft.fft(modes[i, :])) ** 2
        # take only half of the points, and reject 1st point
        fftmodes = (2.0 / nrec) * tmp[1:npfft + 1]

        for j in range(ngain):
            phaseError = np.sum(Hcor[j, :] * fftmodes)

            # Initializes things at startup
            if(j == 0):
                minPhaseError = phaseError
                jmin = j

            # Search for the minimum value, and the associated gain
            if(phaseError < minPhaseError):
                minPhaseError = phaseError
                jmin = j


        optimumGain[i] = G[jmin]

    return optimumGain


def modalControlOptimizationClosedLoopData(slopesClosed, voltsData, mia, S2M, M2V,
                                           gmax = 0.5, Fs = 500, latency = 1, BP = 1e12):
    """
    Input arguments :
    <slopesClosed> : set of N frames of slopes (i.e. centroids), (nSlopes x (N + 1))
    <voltsData>   : set of N frames of volts (i.e. DM commands) synchronized with slopes, (nActu x (N + 1))
    <mia>         : measured interaction matrix (nSlopes x nActu)
    <S2M>         : modal control matrix
    <M2V>         : matrix of system modes
    <gmax>        : scalar floating point value, maximum gain.
    <Fs>          : sampling frequency (Hertz)
    <latency>     : latency, in seconds, i.e. time between the start of WFS integration and time
                    when the actuator reaches 50% of its command.
    <BP>          : scalar, cutoff frequency of the DM (seen as a 1st order filter)

    Returned value :
    Array of optimal gains on the 50 modes.

    """
    # determines how many frames are recorded in the data set,
    # will be used later (could also be passed as an argument..)
    (nrec, ns) = slopesClosed.shape

    # Reconstruction of pseudo open-loop slopes
    slopesOpen = slopesClosed[:, 1:] - np.dot(mia, voltsData[:, :-1])

    return modalControlOptimizationOpenLoopData(slopesOpen, S2M, M2V,
                                                gmax = gmax, Fs = Fs, latency = latency, BP = BP)



def computeSystemOptimizedControlMatrix(S2M, gainVector, M2V):
    """

    Input parameters :
    <S2M>        is the modal control matrix as computed by function computeModalControlMatrix()
    <gainVector> a 1D vector of modal gains, supposed to be provided by function modalControlOptimization()
    <M2V>        is the matrix of the modes computed by the function createModalBasis()

    Returned parameter : the command matrix of the system.
    """
    # let's allocate memory space
    tmp = np.zeros(S2M.shape)
    # let's multiply by the gains (dimensions of gainVector
    # and S2M must comply)
    ngain = gainVector.shape[0]
    for i in range(ngain):
        tmp[i, :] = S2M[i, :] * gainVector[i]

    # matrix multiply  M2V * (gains*S2M)
    mc = np.dot(M2V, tmp)

    return mc
