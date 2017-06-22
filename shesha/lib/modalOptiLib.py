'''
Created on 23 janv. 2017

@author: vdeo
'''

import numpy as np
import time
import matplotlib.pyplot as plt
import os, sys

try:
    HRAAPATH = os.environ["HRAAPATH"]
    PYRDATAPATH = HRAAPATH + "/PYRCADO/data"
    sys.path.insert(0, HRAAPATH + '/PYRCADO/Python')
except:
    print("ERROR COULD NOT FIND HRAAPATH environment variable. Please set it up in the .bashrc in re-start")
    pass





def applyVoltGetSlopes(wao, volts, withAtm, extPyrc = None):
    '''
        Applies voltages
        on the dms and reads out wfs slopes
    :param wao: AO Widget
    :param withAtm: Show atmosphere to wfs or only DM phase
    '''
    # FIXME add a refSlopes somewhere ?
    wao.dms.set_full_comm(volts.astype(np.float32).copy())
    for w in range(len(wao.config.p_wfss)):
        wao.wfs.reset_phase(w)
        if withAtm:
            wao.wfs.sensors_trace(w, "all", wao.tel, wao.atm, wao.dms)
        else:
            wao.wfs.sensors_trace(w, "dm", wao.tel, wao.atm, wao.dms)
        wao.wfs.sensors_compimg(w, noise = False)

    if extPyrc is None:
        wao.rtc.docentroids(0)
        slopes = wao.rtc.getCentroids(0)
    else:
        import controlRoutines as cr
        pyrhr = wao.wfs.get_pyrimghr(0)
        _, _, _, grad = cr.getPyramidData(pyrhr, extPyrc)
        slopes = cr.gradToSlopes(grad, extPyrc["validZone"])
    return slopes


def applyTiltsGetFlatSlopes(wao, TTpush, extPyrc):
    pyrhr = None
    for i in [-1, 1]:
        for j in [-1, 1]:
            wao.dms.set_comm('tt', 0, np.array([i, j], np.float32) * TTpush,
                             shape_dm = True)
            
            for w in range(len(wao.config.p_wfss)):
                wao.wfs.reset_phase(w)
                wao.wfs.sensors_trace(w, "dm", wao.tel, wao.atm, wao.dms)
                wao.wfs.sensors_compimg(w, noise = False)

            import controlRoutines as cr
            if pyrhr is None:
                pyrhr = wao.wfs.get_pyrimghr(0) / 4.
            else:
                pyrhr += wao.wfs.get_pyrimghr(0) / 4.

    _, _, _, grad = cr.getPyramidData(pyrhr, extPyrc)
    slopes = cr.gradToSlopes(grad, extPyrc["validZone"])
    return slopes



def measureIMatKLPP(wao, ampliVec, KL2V, nSlopes, withAtm, extPyrc = None):
    '''
        Make modal interaction matrix using push-pull normalized difference
    :param wao: AO Widget
    :param ampliVec: nActu vector of push for imat values
    :param KL2V: Matrix of modes in the DM space
    :param nSlopes: WFS output dimension
    :param withAtm: Make interaction matrix around currently shown atmosphere (True) or around null phase (False)
    '''
    currentVolts = wao.rtc.getVoltage(0)
    if not withAtm:
        currentVolts[:] = 0.

    iMatKL = np.zeros((nSlopes, KL2V.shape[1]))

    st = time.time()

    ref = applyVoltGetSlopes(wao, currentVolts, withAtm, extPyrc = extPyrc)
    

    for kl in range(KL2V.shape[1]):
        v = ampliVec[kl] * KL2V[:, kl]

        iMatKL[:, kl] = (applyVoltGetSlopes(wao, v + currentVolts, withAtm, extPyrc = extPyrc) - ref) / ampliVec[kl]

        print("Doing KL interaction matrix on mode: #%d\r" % kl, end=' ')
        os.sys.stdout.flush()

    print("Modal interaction matrix done in %3.1f seconds" % (time.time() - st))
    return iMatKL


def measureIMatKLSine(wao, ampliVec, KL2V, nSlopes, withAtm, extPyrc = None):

    if withAtm:
        currentVolts = wao.rtc.getVoltage(0)
    else:
        currentVolts = 0.0

    primeFreq = nPrimes(KL2V.shape[1])
    primeFreq = np.array(primeFreq[2:] + primeFreq[:2])  # Assign low frequencies to tip-tilt
    nFrames = 2 * (np.max(primeFreq) + 1)
    freqArr = 1.0 * primeFreq * np.pi / (np.max(primeFreq) + 1)

    frames = np.zeros((nSlopes, nFrames))

    ampliKL = ampliVec * KL2V
    voltsToSet = np.dot(ampliKL, np.cos(np.dot(freqArr[..., np.newaxis], np.arange(nFrames)[np.newaxis, ...])))

    st = time.time()
    for f in range(nFrames):
        v = voltsToSet[:, f]
        frames[:, f] = applyVoltGetSlopes(wao, v + currentVolts, withAtm, extPyrc = extPyrc)
        print("Doing sine KL interaction matrix on frame: #%d / %d\r" % (f + 1, nFrames), end=' ')
        os.sys.stdout.flush()

    print("\nModal sine interaction recording done in %3.1f seconds" % (time.time() - st))

    FTImat = np.fft.fft(frames, axis = 1)
    imatKL = np.real(FTImat[:, primeFreq]) / ((np.max(primeFreq) + 1) * ampliVec)

    return imatKL


def nPrimes(n):
    def prime(i, primes):
        for prime in primes:
            if not (i == prime or i % prime):
                return False
        primes.append(i)
        return i

    primes = []
    i, p = 2, 0
    while True:
        if prime(i, primes):
            p += 1
            if p == n:
                return primes
        i += 1


def normalizeKL2V(KL2V, mode = 'linf'):
    '''
        L-Inf normalization of mirror modes
    :param KL2V: Matrix of modes in the DM Space
    '''
    if mode == 'linf':
        return KL2V / np.amax(np.abs(KL2V), axis = 0)[np.newaxis, ...]
    elif mode == 'l2':
        return KL2V / np.sqrt(np.sum(KL2V ** 2, axis = 0))[np.newaxis, ...]
    else:
        print("Invalid normalize KL2V mode - required optional mode = 'l2' or 'linf'")
        return None


def cropKL2V(KL2V):
    '''
        Remove the least 2 KL modes of the pzt mirror (columns -4:-2) to remove Tip-Tilt redundancy
    :param KL2V: Matrix of modes in the DM Space
    '''
    return KL2V[:, list(range(KL2V.shape[1] - 4)) + [-2, -1]]


def plotVDm(wao, numdm, V, size = 16, fignum = False):
    if(wao.config.p_dms[numdm]._j1.shape[0] != V.shape[0]):
        print("Error V should have %d dimension not %d " % (wao.config.p_dms[numdm]._j1.shape[0], V.shape[0]))
    else:
        if(fignum):
            plt.figure(fignum)
        plt.clf()
        plt.scatter(wao.config.p_dms[numdm]._j1, wao.config.p_dms[
                    numdm]._i1, c = V, marker = "s", s = size ** 2)
        plt.colorbar()


def computeImatKL(wao, pushPZT, pushTT, KL2V, nSlopes, withAtm, extPyrc = None):
    '''
        Compute the modal interaction matrix with given input parameters
    :param wao: AO Widget
    :param pushPZT: push for Imat value for PZT DM -- MICRON --
    :param pushTT: push for Imat value for TT mirror -- ARCSEC --
    :param KL2V: Matrix of modes in the DM space
    :param nSlopes: Dimension of WFS slope space
    :param withAtm: Make interaction matrix around currently shown atmosphere (True) or around null phase (False)
    '''
    NKL = KL2V.shape[1]
    modesAmpli = np.array([pushPZT] * (NKL - 2) + [pushTT] * 2)
    iMat = measureIMatKLPP(wao, modesAmpli, KL2V, nSlopes, withAtm, extPyrc = extPyrc)
    return iMat


def computeCmatKL(iMatKL, nFilt):
    '''
        Invert the interaction matrix into a mode filtered control matrix
    :param iMatKL: interaction matrix nSlopes x nModes
    :param nFilt: number of modes to filter before inversion
    '''
    nModes = (iMatKL.shape[1] - nFilt)
    iMatKLFilt = iMatKL[:, list(range(0, nModes - 2)) + [-2, -1]]

    cMatKL = np.linalg.inv(iMatKLFilt.T.dot(iMatKLFilt)).dot(iMatKLFilt.T)

    return cMatKL


def cMatFromcMatKL(cMatKL, KL2V, gainVector = None):
    '''
        Compute the slope -> volts control matrix from the slope -> KL modal control matrix
    :param cMatKL: Modal control matrix (nModes - filtered) x nSlopes
    :param KL2V: Matrix of modes in the DM space
    :param gainVector: Vector of modal gains to apply on the modal control matrix rows
    '''
    if gainVector is None:
        gainVector = np.ones((cMatKL.shape[0]))

    cMat = np.dot(KL2V, gainVector[..., np.newaxis] * cMatKL)

    return cMat


class ModalGainOptimizer:
    '''
        ModalGainOptimizer class
        Utility class to compute custom iMat/cMat around zero point or around current atmosphere
        within a Compass/Canapass instance
    '''


    def __init__(self, wao, extPyrc = None):
        '''
            ModalGainOptimizer constructor.
            Resets AO session's RTC to some defaults.
            Gets modal basis and sets up some storage
        :param wao: Handle to Compass/Canapass session
        '''
        self.wao = wao

        self.nActu = sum(wao.config.p_rtc.controllers[0].nactu)
        self.waoSlopes = wao.rtc.getCentroids(0).shape[0]
        if extPyrc is None:
            self.nSlopes = self.waoSlopes
        else:
            self.nSlopes = applyVoltGetSlopes(wao, np.zeros(self.nActu), False, extPyrc).shape[0]
        self.pushPzt = 0.01  # 10 nm
        self.pushTT = 0.005  # 5 mas

        self.gain = 0.7
        self.nFilter = 10

        self.pyrc = extPyrc

        # Clear WAO and do some controller setup
        self.wao.aoLoopClicked(False)
        self.wao.set_atmos(False)
        self.wao.openLoop()

        self.wao.setIntegratorLaw()
        self.wao.rtc.set_decayFactor(0, np.ones(self.waoSlopes, dtype = (np.float32)))
        self.wao.rtc.set_mgain(0, self.gain * np.ones(self.nActu, dtype = (np.float32)))
        
        self.initBasis(wao.returnkl2V())


    
    def initBasis(self, rawKL2V, norm = False):
        '''
            Sets and normalizes a DM basis
        :param rawKL2V: Self-explanatory.
        '''
        if norm:
            self.KL2V = normalizeKL2V(rawKL2V, mode = 'l2')
        else:
            self.KL2V = rawKL2V.copy()

        self.KL2VFilt = \
            self.KL2V[:, list(range(0, self.KL2V.shape[1] - self.nFilter - 2)) + [-2, -1]]

        self.iMatKLRef = None
        self.cMatKLRef = None

        self.gainValues = None

    def computeRef(self):
        '''
            Compute modal interaction and control matrices around 0 atmosphere point.
            Store them within ModalGainOptimizer instance
        '''
        self.wao.closeLoop()
        self.iMatKLRef = computeImatKL(self.wao, self.pushPzt, self.pushTT,
                                       self.KL2V, self.nSlopes, withAtm = False, extPyrc = self.pyrc)
        self.cMatKLRef = computeCmatKL(self.iMatKLRef, self.nFilter)

        self.gainValues = np.ones(self.cMatKLRef.shape[0])

    def setRef(self):
        '''
            Sets the reference control matrix into the AO session RTC.
        '''
        cMat = cMatFromcMatKL(self.cMatKLRef, self.KL2VFilt)
        self.wao.rtc.set_cmat(0, cMat.astype(np.float32).copy())


    def setGeomRef(self):
        '''
            Sets the geom control matrix into the AO session RTC
        '''
        cMat = cMatFromcMatKL(self.cMatKLGeom, self.KL2VFilt)
        self.wao.rtc.set_cmat(0, cMat.astype(np.float32).copy())


    def computeAtmGains(self):
        '''
            Compute a new modal interaction matrix around the current residual atmosphere setpoint.
            Apply the sensitivity compensation coefficients found to the zero-point modal control matrix
            Set the updated control matrix to the AO session.
        '''
        self.wao.closeLoop()
        iMatKL = computeImatKL(self.wao, self.pushPzt, self.pushTT,
                               self.KL2V, self.nSlopes, withAtm = True, extPyrc = self.pyrc)
        ksiVal = np.sqrt(np.diag(np.dot(self.iMatKLRef.T, self.iMatKLRef)) /
                         np.diag(np.dot(iMatKL.T, iMatKL)))

        ksiFilter = ksiVal[list(range(iMatKL.shape[1] - self.nFilter - 2)) + [-2, -1]]
#         ksiFilter = ksiFilter * np.sqrt(np.shape(ksiFilter)[0] / np.sum(ksiFilter ** 2))

        return ksiFilter, iMatKL

    def setCmat(self, ksis = None):
        if ksis is None:
            ksis = self.gainValues
        else:
            self.gainValues = ksis

        cMat = cMatFromcMatKL(self.cMatKLRef, self.KL2VFilt, gainVector = ksis)
        self.wao.rtc.set_cmat(0, cMat.astype(np.float32).copy())


def refreshAtmos(wao, nSteps = 10000):
    for i in range(nSteps):
        wao.atm.move_atmos()

def loopStaticAtmos(wao, nIter):
    for i in range(nIter):
        for w in range(len(wao.config.p_wfss)):
            wao.wfs.sensors_trace(w, "all", wao.tel, wao.atm, wao.dms)
            wao.wfs.sensors_compimg(w)

        wao.rtc.docentroids(0)
        wao.rtc.docontrol(0)
        wao.rtc.applycontrol(0, wao.dms)

    for t in range(wao.config.p_target.ntargets):
        wao.tar.atmos_trace(t, wao.atm, wao.tel)
        wao.tar.dmtrace(t, wao.dms)

        
def AOLoop(nIter, wao, cmat = None, pyrc = None,
           reset = None, gain = 0.7, clip = -1, refSlopes = None,
           moveAtmos = True, showAtmos = True, computeLE = None,
           computeResidualsM2PH = None):
    '''
    
    :param nIter:
    :param cmat: None to use wao.rtc, or a control matrix to externalize computation
    :param wao:
    :param reset: None to reset or an initial command array otherwise
    :param sleep_time:
    :param refSlopes: None for zero refslopes or slopes array
    :param moveAtmos:
    :param showAtmos:
    :param showOffset:
    :param computeLE: None or [target numbers]
    :param computeResidualsM2PH: M2PH phase cube for residual, or None for no residuals
    '''
    
    if cmat is None:  # Internal WFS/RTC
        nactu = 0
        for p_dm in wao.config.p_dms:
            nactu += p_dm._ntotact
    else:  # Externalize WFS/RTC
        nactu = cmat.shape[0]
        
    # MACROS
    def trace_wfs(showAtmos):
        for w in range(len(wao.config.p_wfss)):
            if showAtmos:
                wao.wfs.sensors_trace(w, "all", wao.tel, wao.atm, wao.dms, rst = 1)
            else:
                wao.wfs.sensors_trace(w, "dm", wao.tel, wao.atm, wao.dms, rst = 1)
            wao.wfs.sensors_compimg(w)
    def trace_tar(showAtmos):
        for t in range(wao.config.p_target.ntargets):
            wao.tar.reset_phase(t)
            if showAtmos:
                wao.tar.atmos_trace(t, wao.atm, wao.tel)
            wao.tar.dmtrace(t, wao.dms)

    if reset is None:
        wao.dms.resetdm('pzt', 0)
        wao.dms.resetdm('tt', 0)
        currCommand = np.zeros(nactu)
    else:
        currCommand = reset.copy()
        wao.dms.set_full_comm(currCommand.copy().astype(np.float32))
    trace_wfs(showAtmos)
    trace_tar(showAtmos)

    if refSlopes is None:
        if cmat is None:
            refSlopes = np.zeros(wao.rtc.getCentroids(0).shape[0])
        else:
            refSlopes = np.zeros(cmat.shape[1])

    longExp = []
    if computeLE is not None:
        for i in range(len(computeLE)):
            longExp += [wao.tar.get_image(computeLE[i], "se")]


    if computeResidualsM2PH is not None:
        residualTot = np.zeros((nIter, computeResidualsM2PH.shape[0] + 1))

    for i in range(nIter):

        if moveAtmos:
            wao.atm.move_atmos()

        if cmat is None:  # Internal WFS/RTC
            wao.rtc.applycontrol(0, wao.dms)
            trace_wfs(showAtmos)
            wao.rtc.docentroids(0)
            wao.rtc.docontrol(0)
            currCommand = np.r_[wao.dms.getComm('pzt', 0), wao.dms.getComm('tt', 0)]

        else:  # Externalize WFS/RTC
            slopes = applyVoltGetSlopes(wao, currCommand, showAtmos, extPyrc = pyrc)

            currCommand -= gain * cmat.dot(slopes - refSlopes)
            if clip > 0:
                currCommand = np.clip(currCommand, -clip, clip)

        if computeResidualsM2PH is not None:
            residualTot[i, :] = residuals(computeResidualsM2PH, phMap(wao, showAtmos))

        if computeLE is not None:
            trace_tar(showAtmos)
            for i in range(len(computeLE)):
                longExp[i] += wao.tar.get_image(computeLE[i], "se")
        
    arg = 2 * (computeLE is not None) + (computeResidualsM2PH is not None)
    if arg == 0:
        return currCommand
    elif arg == 1:
        return currCommand, residualTot
    elif arg == 2:
        return currCommand, longExp
    else:
        return currCommand, longExp, residualTot
        
    for t in range(wao.config.p_target.ntargets):
        wao.tar.atmos_trace(t, wao.atm, wao.tel)
        wao.tar.dmtrace(t, wao.dms)
    

def phMap(wao, showAtmos):
    '''
        Get a WFS phase map
    '''
    wao.wfs.reset_phase(0)
    if showAtmos:
        wao.wfs.sensors_trace(0, "all", wao.tel, wao.atm, wao.dms)
    else:
        wao.wfs.sensors_trace(0, "dm", wao.tel, wao.atm, wao.dms)

    pupil = wao.config.p_geom.get_mpupil().astype(bool)
    phase = wao.wfs.get_phase(0) * pupil
    phase[pupil] -= np.mean(phase[pupil])
    phase /= np.sqrt(np.sum(pupil)) * wao.config.p_wfs0.Lambda

    return phase

def residuals(M2PH, phaseMap):
    '''
        Project a phase map on the M2PH tensor
    '''
    # M2PH : nModes * npup * npup
    # phaseMap : npup * npup

    coefficients = np.tensordot(M2PH, phaseMap, axes = 2)
    fitting = np.sqrt(np.sum(phaseMap ** 2) - np.sum(coefficients ** 2))
    return np.r_[coefficients, fitting]




def simuAndRecord(wao, nIter):
    
    nSlopes = wao.rtc.getCentroids(0).shape[0]
    
    slopes = np.zeros((nSlopes, nIter))

    for i in range(nIter):
        slopes[:, i] = doubleWFSLoop(wao)
    
    return slopes

def doubleWFSLoop(wao):
    wao.atm.move_atmos()

    for t in range(wao.config.p_target.ntargets):
        wao.tar.atmos_trace(t, wao.atm, wao.tel)
        wao.tar.dmtrace(t, wao.dms)
    for w in range(len(wao.config.p_wfss)):
        wao.wfs.sensors_trace(w, "all", wao.tel, wao.atm, wao.dms)
        wao.wfs.sensors_compimg(w)

    # Pyramid slopes

    wao.rtc.docentroids(0)
    slopes = wao.rtc.getCentroids(0)

    wao.rtc.docontrol(0)
    wao.rtc.applycontrol(0, wao.dms)
#     wao.tar.get_strehl(0)[0]

    return slopes


def mode2ph(KL2V, wao):
    nModes = KL2V.shape[1]
    pupil = wao.config.p_geom.get_mpupil().astype(bool)
    
    M2PH = np.zeros((nModes, pupil.shape[0], pupil.shape[1]))
    S = np.sum(pupil)

    for i in range(nModes):
        wao.wfs.reset_phase(0)
        wao.dms.set_full_comm(KL2V[:, i].astype(np.float32).copy())
        wao.wfs.sensors_trace(0, "dm", wao.tel, wao.atm, wao.dms)
        M2PH[i] = wao.wfs.get_phase(0) * pupil
        M2PH[i][pupil] -= np.mean(M2PH[i][pupil])
        M2PH[i] /= np.sqrt(np.sum(M2PH[i][pupil] ** 2))
        print(str(i + 1) + "/" + str(nModes))

    return M2PH




def updateToOptimalCmat(wao, slopes, volts, miKL, S2KL, M2V, offset):

    nfilt = 20

    mia = np.dot(miKL, np.linalg.pinv(M2V))

    gainVector = modalControlOptimizationClosedLoopData(
        slopes, volts, mia, S2KL, M2V, offset, gmax = 0.5)

    if nfilt == 0:
        optiCmat = cMatFromcMatKL(S2KL, M2V, gainVector = gainVector)
    else:
        optiCmat = cMatFromcMatKL(
            S2KL, M2V[:, list(range(0, (M2V.shape[1] - nfilt - 2))) + [-2, -1]], gainVector = gainVector)

    wao.openLoop()
    wao.rtc.set_cmat(0, optiCmat.astype(np.float32).copy())

    wao.closeLoop()

    return gainVector


def main(wao, withAtm, nFrames):

    wao.aoLoopClicked(False)
    wao.set_atmos(False)

    (miKL, S2KL, KL2VN) = setupRef(wao, withAtm)

    nSlopes = S2KL.shape[1]
    nActu = KL2VN.shape[0]

    wao.closeLoop()
    wao.set_atmos(False)

    return miKL, np.dot(miKL, np.linalg.pinv(KL2VN))

    offset = 2
    slopes = np.zeros((nSlopes, nFrames + offset), dtype = np.float32)
    volts = np.zeros((nActu, nFrames + offset), dtype = np.float32)

    for i in range(nFrames + 1):
        wao.mainLoop()
        slopes[:, i] = wao.getSlopes()
        volts[:, i] = wao.getVolts()

    gainVector = updateToOptimalCmat(
        wao, slopes, volts, miKL, S2KL, KL2VN, offset)

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
    # zero-order hold with 1/2 frame delay
    Hccd = (1. - np.exp(-p * Te)) / (p * Te)
    Hdac = Hccd  # well, same.
    tdelay = latency - Te  # time between END of the integration and start of command
    Hret = np.exp(-p * tdelay)  # latency transfer function
    # transfer func of the DM, as a 1st order filter
    Hmir = 1. / (1. + 1j * freq / BP)
    # open-loop transfer function
    Hbo = Hint * Hccd * Hdac * Hret * Hmir
    # correction transfer function
    Hcor = 1. / np.abs(1 + G * Hbo) ** 2
    return Hcor


def modalControlOptimizationOpenLoopData(slopes, S2M, M2V, gmax = 0.5, Fs = 500, latency = 1, BP = 1e12):
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
    # TBC, maybe we'll keep gmin=0.001 to allow for static aberration
    # compensation, if any. TBC during AIT.
    gmin = 0.0001
    # 1D array from gmin to gmax in ngain points
    G = np.exp(np.linspace(np.log(gmin), np.log(gmax), ngain))

    npfft = nrec // 2  # number of useful points in the FFT
    # create a 1D array of frequency ranging from Fs/nrec to Fs/2.0 in npfft
    # points
    freq = np.linspace(Fs * 1.0 / nrec, Fs / 2.0, npfft)
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
    # STOP
    return optimumGain


def modalControlOptimizationClosedLoopData(slopesClosed, voltsData, mia, S2M, M2V, offset,
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
    slopesOpen = slopesClosed[:, offset:] - np.dot(mia, voltsData[:, :-offset])

    return modalControlOptimizationOpenLoopData(slopesOpen, S2M, M2V,
                                                gmax = gmax, Fs = Fs, latency = latency, BP = BP)











