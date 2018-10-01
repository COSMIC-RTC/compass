"""Widget to simulate a closed loop

Usage:
  canapassSupervisor.py [<parameters_filename>]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
"""

import os, sys
import numpy as np
import time
from collections import OrderedDict

from tqdm import trange
from tqdm import tqdm

import astropy.io.fits as pfits
from threading import Thread
from subprocess import Popen, PIPE

import shesha.ao as ao
import shesha.constants as scons

from typing import Any, Dict, Tuple, Callable, List
from .compassSupervisor import CompassSupervisor

# from naga.obj import obj_Double2D
# from naga.magma import syevd_Double, svd_host_Double
# from naga.context import context as naga_context

# from naga.host_obj import host_obj_Double1D, host_obj_Double2D


class CanapassSupervisor(CompassSupervisor):

    def __init__(self, configFile: str=None, BRAHMA: bool=True) -> None:
        CompassSupervisor.__init__(self, configFile=configFile, BRAHMA=BRAHMA)

        #############################################################
        #                 CONNECTED BUTTONS                         #
        #############################################################
        # Default path for config files

        self.ph2modes = None
        self.KL2V = None
        self.P = None
        self.currentBuffer = 1
        #############################################################
        #                       METHODS                             #
        #############################################################

    """
          ____ ___  __  __ ____   _    ____ ____
         / ___/ _ \|  \/  |  _ \ / \  / ___/ ___|
        | |  | | | | |\/| | |_) / _ \ \___ \___ \
        | |__| |_| | |  | |  __/ ___ \ ___) |__) |
         \____\___/|_|  |_|_| /_/   \_\____/____/
         ____  _   _ ____  _____ ______     _____ ____   ___  ____
        / ___|| | | |  _ \| ____|  _ \ \   / /_ _/ ___| / _ \|  _ \
        \___ \| | | | |_) |  _| | |_) \ \ / / | |\___ \| | | | |_) |
         ___) | |_| |  __/| |___|  _ < \ V /  | | ___) | |_| |  _ <
        |____/ \___/|_|   |_____|_| \_\ \_/  |___|____/ \___/|_| \_\
         __  __ _____ _____ _   _  ___  ____  ____
        |  \/  | ____|_   _| | | |/ _ \|  _ \/ ___|
        | |\/| |  _|   | | | |_| | | | | | | \___ \
        | |  | | |___  | | |  _  | |_| | |_| |___) |
        |_|  |_|_____| |_| |_| |_|\___/|____/|____/
    """

    def getConfig(self, path=None):
        ''' Returns the configuration in use, in a supervisor specific format '''
        if path:
            self.writeConfigOnFile(path)
            return
        return self._sim.config

    def loadConfig(self, configFile: str=None, sim=None) -> None:
        ''' Load the configuration for the compass supervisor'''
        CompassSupervisor.loadConfig(self, configFile=configFile, sim=sim)
        print("switching to a generic controller")
        self._sim.config.p_controllers[0].type = scons.ControllerType.GENERIC

    """
          ____    _    _   _    _    ____   _    ____ ____
         / ___|  / \  | \ | |  / \  |  _ \ / \  / ___/ ___|
        | |     / _ \ |  \| | / _ \ | |_) / _ \ \___ \___ \
        | |___ / ___ \| |\  |/ ___ \|  __/ ___ \ ___) |__) |
         \____/_/   \_\_| \_/_/   \_\_| /_/   \_\____/____/
         __  __ _____ _____ _   _  ___  ____  ____
        |  \/  | ____|_   _| | | |/ _ \|  _ \/ ___|
        | |\/| |  _|   | | | |_| | | | | | | \___ \
        | |  | | |___  | | |  _  | |_| | |_| |___) |
        |_|  |_|_____| |_| |_| |_|\___/|____/|____/
    """

    def computePh2Modes(self):
        oldnoise = self._sim.config.p_wfs0.noise
        self.setNoise(-1)

        if (self.modalBasis is None):
            self.modalBasis, self.P = self.getModes2VBasis("Btt")
        nbmode = self.modalBasis.shape[1]
        pup = self._sim.config.p_geom._spupil
        ph = self._sim.tar.get_phase(0)
        ph2KL = np.zeros((nbmode, ph.shape[0], ph.shape[1]))
        S = np.sum(pup)
        for i in trange(nbmode):
            self.resetTarPhase(0)
            self._sim.dms.set_full_comm(
                    (self.modalBasis[:, i]).astype(np.float32).copy())
            self._sim.next(see_atmos=False)
            ph = self.getTarPhase(0) * pup
            # Normalisation pour les unites rms en microns !!!
            norm = np.sqrt(np.sum((ph)**2) / S)
            ph2KL[i] = ph / norm
        self.ph2modes = ph2KL
        self.setNoise(oldnoise)
        return ph2KL

    def next(self, nbiters, see_atmos=True):
        for i in trange(nbiters):
            self._sim.next(see_atmos=see_atmos)

    def loop(self, n: int=1, monitoring_freq: int=100, **kwargs):
        """
        Perform the AO loop for n iterations

        :parameters:
            n: (int): (optional) Number of iteration that will be done
            monitoring_freq: (int): (optional) Monitoring frequency [frames]
        """
        self._sim.loop(n, monitoring_freq=monitoring_freq)

    def computePh2ModesFits(self, fullpath):
        data = self.computePh2Modes()
        self.writeDataInFits(data, fullpath)

    def getModes2VBasis(self, ModalBasisType, merged=False, nbpairs=None):
        if (ModalBasisType == "KL2V"):
            print("Computing KL2V basis...")
            self.modalBasis, _ = self.returnkl2V()
            return self.modalBasis, 0
        elif (ModalBasisType == "Btt"):
            print("Computing Btt basis...")
            self.modalBasis, self.P = self.compute_Btt2(inv_method="cpu_svd",
                                                        merged=merged, nbpairs=nbpairs)
            return self.modalBasis, self.P

    def returnkl2V(self):
        """
        KL2V = ao.compute_KL2V(self._sim.config.p_controllers[
                               0], self._sim.dms, self._sim.config.p_dms, self._sim.config.p_geom, self._sim.config.p_atmos, self._sim.config.p_tel)
        """
        if (self.KL2V is None):
            print("Computing KL2V...")
            KL2V = ao.compute_KL2V(self._sim.config.p_controllers[0], self._sim.dms,
                                   self._sim.config.p_dms, self._sim.config.p_geom,
                                   self._sim.config.p_atmos, self._sim.config.p_tel)
            print("KL2V Done!")
            self.KL2V = KL2V
            return KL2V, 0
        else:
            return self.KL2V, 0

    def setGain(self, gain: float) -> None:
        CompassSupervisor.setGain(self, gain)

    def computeMerged(self, nbpairs=None):
        p_geom = self._sim.config.p_geom
        import shesha.util.make_pupil as mkP
        import shesha.util.utilities as util
        import scipy.ndimage

        cent = p_geom.pupdiam / 2. + 0.5
        p_tel = self._sim.config.p_tel
        p_tel.t_spiders = 0.51
        spup = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent,
                              cent).astype(np.float32).T

        p_tel.t_spiders = 0.
        spup2 = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent,
                               cent).astype(np.float32).T

        spiders = spup2 - spup

        (spidersID, k) = scipy.ndimage.label(spiders)
        spidersi = util.pad_array(spidersID, p_geom.ssize).astype(np.float32)
        pxListSpider = [np.where(spidersi == i) for i in range(1, k + 1)]

        # DM positions in iPupil:
        dmposx = self._sim.config.p_dm0._xpos - 0.5
        dmposy = self._sim.config.p_dm0._ypos - 0.5
        dmposMat = np.c_[dmposx, dmposy].T  # one actu per column

        pitch = self._sim.config.p_dm0._pitch
        DISCARD = np.zeros(len(dmposx), dtype=np.bool)
        PAIRS = []

        # For each of the k pieces of the spider
        for k, pxList in enumerate(pxListSpider):
            pts = np.c_[pxList[1], pxList[0]]  # x,y coord of pixels of the spider piece
            # lineEq = [a, b]
            # Which minimizes leqst squares of aa*x + bb*y = 1
            lineEq = np.linalg.pinv(pts).dot(np.ones(pts.shape[0]))
            aa, bb = lineEq[0], lineEq[1]

            # Find any point of the fitted line.
            # For simplicity, the intercept with one of the axes x = 0 / y = 0
            if np.abs(bb) < np.abs(aa):  # near vertical
                onePoint = np.array([1 / aa, 0.])
            else:  # otherwise
                onePoint = np.array([0., 1 / bb])

            # Rotation that aligns the spider piece to the horizontal
            rotation = np.array([[-bb, aa], [-aa, -bb]]) / (aa**2 + bb**2)**.5

            # Rotated the spider mask
            rotatedPx = rotation.dot(pts.T - onePoint[:, None])
            # Min and max coordinates along the spider length - to filter actuators that are on
            # 'This' side of the pupil and not the other side
            minU, maxU = rotatedPx[0].min() - 5. * pitch, rotatedPx[0].max() + 5. * pitch

            # Rotate the actuators
            rotatedActus = rotation.dot(dmposMat - onePoint[:, None])
            selGoodSide = (rotatedActus[0] > minU) & (rotatedActus[0] < maxU)
            seuil = 0.05
            # Actuators below this piece of spider
            selDiscard = (np.abs(rotatedActus[1]) < seuil * pitch) & selGoodSide
            DISCARD |= selDiscard

            # Actuator 'near' this piece of spider
            selPairable = (np.abs(rotatedActus[1]) > seuil  * pitch) & \
                            (np.abs(rotatedActus[1]) < 1. * pitch) & \
                            selGoodSide

            pairableIdx = np.where(selPairable)[0]  # Indices of these actuators
            uCoord = rotatedActus[
                    0, selPairable]  # Their linear coord along the spider major axis

            order = np.sort(uCoord)  # Sort by linear coordinate
            orderIdx = pairableIdx[np.argsort(
                    uCoord)]  # And keep track of original indexes

            # i = 0
            # while i < len(order) - 1:
            if (nbpairs is None):
                i = 0
                ii = len(order) - 1
            else:
                i = len(order) // 2 - nbpairs
                ii = len(order) // 2 + nbpairs
            while (i < ii):
                # Check if next actu in sorted order is very close
                # Some lonely actuators may be hanging in this list
                if np.abs(order[i] - order[i + 1]) < .2 * pitch:
                    PAIRS += [(orderIdx[i], orderIdx[i + 1])]
                    i += 2
                else:
                    i += 1
        print('To discard: %u actu' % np.sum(DISCARD))
        print('%u pairs to slave' % len(PAIRS))
        if np.sum(DISCARD) == 0:
            DISCARD = []
        else:
            list(np.where(DISCARD)[0])
        return np.asarray(PAIRS), list(np.where(DISCARD)[0])

    def computeActusPetals(self):
        import shesha.util.make_pupil as mkP
        import shesha.util.utilities as util
        N = self._sim.config.p_geom.pupdiam
        i0 = j0 = self._sim.config.p_geom.pupdiam / 2. - 0.5
        pixscale = self._sim.config.p_tel.diam / self._sim.config.p_geom.pupdiam
        dspider = self._sim.config.p_tel.t_spiders
        #spup = self.getSpupil()
        oldspidersize = self._sim.config.p_tel.t_spiders
        self._sim.config.p_tel.set_t_spiders = 0.
        spup = mkP.make_pupil(self._sim.config.p_geom.pupdiam,
                              self._sim.config.p_geom.pupdiam, self._sim.config.p_tel,
                              i0, j0)
        self._sim.config.p_tel.set_t_spiders = oldspidersize
        petalsPups = mkP.compute6Segments(spup, N, pixscale, 0.51 / 2, i0, j0)
        ipup = util.pad_array(spup, self._sim.config.p_geom.ssize).astype(np.float32)
        dmpup = np.zeros_like(ipup)
        dmposx = self._sim.config.p_dm0._xpos
        dmposy = self._sim.config.p_dm0._ypos
        dmpup[dmposx.astype(int), dmposy.astype(int)] = 1
        dmpupPetals = np.zeros((6, dmpup.shape[0], dmpup.shape[1]))
        indActusPetals = []  # actuators to disable because under spider
        for i in range(6):
            dmpupPetals[i, :, :] = dmpup * util.pad_array(petalsPups[i],
                                                          self._sim.config.p_geom.ssize)
            ipetal = []
            for j in range(len(dmposx)):
                if (dmpupPetals[i, :, :][int(dmposx[j]), int(dmposy[j])] == 1):
                    ipetal.append(j)
            indActusPetals.append(ipetal)
        return indActusPetals, petalsPups

    def computeMerged_fab(self):
        print("Computing merged actuators...")
        p_geom = self._sim.config.p_geom
        import shesha.util.make_pupil as mkP
        import shesha.util.utilities as util

        cent = p_geom.pupdiam / 2. + 0.5
        p_tel = self._sim.config.p_tel
        # spider + pupil
        p_tel.t_spiders = 0.51
        spup = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent,
                              cent).astype(np.float32)

        # Without spider pupil
        p_tel.t_spiders = 0.
        spup2 = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent,
                               cent).astype(np.float32)

        # doubling spider size  + l eft and right spiders only
        p_tel.t_spiders = 0.51 * 2.5
        pp = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent,
                            cent).astype(np.float32)
        l, r = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent, cent,
                              halfSpider=1).astype(np.float32)

        spiders = spup2 - spup
        spidersLeft = l - pp
        spidersRight = r - pp
        dmposx = self._sim.config.p_dm0._xpos
        dmposy = self._sim.config.p_dm0._ypos

        # ----------------------------------------
        # Under spider actuators
        # ----------------------------------------
        ipup = util.pad_array(spup, p_geom.ssize).astype(np.float32)

        spidersi = util.pad_array(spiders, p_geom.ssize).astype(np.float32)
        dmpup = np.zeros_like(spidersi)

        dmpup[dmposy.astype(int), dmposx.astype(int)] = 1
        dmpupUnderSpiders = dmpup * spidersi
        indUnderSpiders = []  # actuators to disable because under spider
        for i in range(len(dmposx)):
            if (dmpupUnderSpiders[int(dmposy[i]), int(dmposx[i])] == 1):
                indUnderSpiders.append(i)

        # Left sided spiders actuators
        spidersiLeft = util.pad_array(spidersLeft, p_geom.ssize).astype(np.float32)
        spidersiRight = util.pad_array(spidersRight, p_geom.ssize).astype(np.float32)

        dmpupLeft = np.zeros_like(spidersiLeft)
        dmpupLeft[dmposy.astype(int), dmposx.astype(int)] = 1

        dmpupRight = np.zeros_like(spidersiRight)
        dmpupRight[dmposy.astype(int), dmposx.astype(int)] = 1

        dmpupLeftSpiders = dmpupLeft * spidersiLeft
        dmpupRightSpiders = dmpupRight * spidersiRight
        indLeftSpiders = []  # actuators to disable because under spider
        indRightSpiders = []  # actuators to disable because under spider
        for i in range(len(dmposx)):
            if (i not in indUnderSpiders):
                if (dmpupLeftSpiders[int(dmposy[i]), int(dmposx[i])] == 1):
                    indLeftSpiders.append(i)
                if (dmpupRightSpiders[int(dmposy[i]), int(dmposx[i])] == 1):
                    indRightSpiders.append(i)

        couplesActus = np.zeros((len(indLeftSpiders), 2)).astype(int)
        couplesActusSelect = np.zeros((len(indLeftSpiders), 2)).astype(int)
        for i in range(len(indRightSpiders)):
            d = np.sqrt((dmposx[indRightSpiders[i]] - dmposx[indLeftSpiders])**2 +
                        (dmposy[indRightSpiders[i]] - dmposy[indLeftSpiders])**2)
            indproche = np.where(d == min(d))[0][0]
            couplesActus[i, 0] = indRightSpiders[i]
            couplesActus[i, 1] = indLeftSpiders[indproche]
            couplesActusSelect[i, 0] = i
            couplesActusSelect[i, 1] = indproche
        print("Done")

        return couplesActus, indUnderSpiders

    def compute_Btt2(self, inv_method: str="cpu_svd", merged=False, nbpairs=None):

        IF = self.getIFsparse(1)
        if (merged):
            couplesActus, indUnderSpiders = self.computeMerged(nbpairs=nbpairs)
            IF2 = IF.copy()
            indremoveTmp = indUnderSpiders.copy()
            indremoveTmp += list(couplesActus[:, 1])
            print("Pairing Actuators...")
            for i in tqdm(range(couplesActus.shape[0])):
                IF2[couplesActus[i, 0], :] += IF2[couplesActus[i, 1], :]
            print("Pairing Done")
            boolarray = np.zeros(IF2.shape[0], dtype=np.bool)
            boolarray[indremoveTmp] = True
            self.slavedActus = boolarray
            self.selectedActus = ~boolarray
            self.couplesActus = couplesActus
            self.indUnderSpiders = indUnderSpiders
            IF2 = IF2[~boolarray, :]
            IF = IF2
        else:
            self.slavedActus = None
            self.selectedActus = None
            self.couplesActus = None
            self.indUnderSpiders = None
        n = IF.shape[0]
        N = IF.shape[1]
        T = self.getIFtt(1)

        delta = IF.dot(IF.T).toarray() / N

        # Tip-tilt + piston
        Tp = np.ones((T.shape[0], T.shape[1] + 1))
        Tp[:, :2] = T.copy()  #.toarray()
        deltaT = IF.dot(Tp) / N
        # Tip tilt projection on the pzt dm
        tau = np.linalg.inv(delta).dot(deltaT)

        # Famille generatrice sans tip tilt
        G = np.identity(n)
        tdt = tau.T.dot(delta).dot(tau)
        subTT = tau.dot(np.linalg.inv(tdt)).dot(tau.T).dot(delta)
        G -= subTT

        # Base orthonormee sans TT
        gdg = G.T.dot(delta).dot(G)

        startTimer = time.time()
        if inv_method == "cpu_svd":
            print("Doing SVD (CPU)")
            U, s, _ = np.linalg.svd(gdg)
        # elif inv_method == "gpu_svd":
        #     print("Doing SVD on CPU of a matrix...")
        #     m = gdg.shape[0]
        #     h_mat = host_obj_Double2D(data=gdg, mallocType="pagelock")
        #     h_eig = host_obj_Double1D(data=np.zeros([m], dtype=np.float64),
        #                               mallocType="pagelock")
        #     h_U = host_obj_Double2D(data=np.zeros((m, m), dtype=np.float64),
        #                             mallocType="pagelock")
        #     h_VT = host_obj_Double2D(data=np.zeros((m, m), dtype=np.float64),
        #                              mallocType="pagelock")
        #     svd_host_Double(h_mat, h_eig, h_U, h_VT)
        #     U = h_U.getData().T.copy()
        #     s = h_eig.getData()[::-1].copy()
        # elif inv_method == "gpu_evd":
        #     print("Doing EVD on GPU of a matrix...")
        #     c = naga_context()
        #     m = gdg.shape[0]
        #     d_mat = obj_Double2D(c, data=gdg)
        #     d_U = obj_Double2D(c, data=np.zeros([m, m], dtype=np.float64))
        #     h_s = np.zeros(m, dtype=np.float64)
        #     syevd_Double(d_mat, h_s, d_U)
        #     U = d_U.device2host().T.copy()
        #     s = h_s[::-1].copy()
        else:
            raise "ERROR cannot recognize inv_method"
        print("Done in %fs" % (time.time() - startTimer))
        U = U[:, :U.shape[1] - 3]
        s = s[:s.size - 3]
        L = np.identity(s.size) / np.sqrt(s)
        B = G.dot(U).dot(L)

        # Rajout du TT
        TT = T.T.dot(T) / N  #.toarray()/N
        Btt = np.zeros((n + 2, n - 1))
        Btt[:B.shape[0], :B.shape[1]] = B
        mini = 1. / np.sqrt(np.abs(TT))
        mini[0, 1] = 0
        mini[1, 0] = 0
        Btt[n:, n - 3:] = mini

        # Calcul du projecteur actus-->modes
        delta = np.zeros((n + T.shape[1], n + T.shape[1]))
        delta[:-2, :-2] = IF.dot(IF.T).toarray() / N
        delta[-2:, -2:] = T.T.dot(T) / N
        P = Btt.T.dot(delta)
        if (merged):
            Btt2 = np.zeros((len(boolarray) + 2, Btt.shape[1]))
            Btt2[np.r_[~boolarray, True, True], :] = Btt
            Btt2[couplesActus[:, 1], :] = Btt2[couplesActus[:, 0], :]

            P2 = np.zeros((Btt.shape[1], len(boolarray) + 2))
            P2[:, np.r_[~boolarray, True, True]] = P
            P2[:, couplesActus[:, 1]] = P2[:, couplesActus[:, 0]]
            return Btt2.astype(np.float32), P.astype(np.float32)
        else:
            return Btt.astype(np.float32), P.astype(np.float32)

    def doImatModal(self, ampliVec, KL2V, Nslopes, noise=False, nmodesMax=0,
                    withTurbu=False, pushPull=False):
        iMatKL = np.zeros((KL2V.shape[1], Nslopes))
        #currentVolts = self._sim.rtc.get_voltage(0).copy()[None,:]

        if (nmodesMax):
            KLMax = nmodesMax
        else:
            KLMax = KL2V.shape[1]
        for kl in trange(KLMax, desc="Modal IM"):
            v = ampliVec[kl] * KL2V[:, kl:kl + 1].T.copy()
            if ((pushPull is True) or
                (withTurbu is True)):  # with turbulence/aberrations => push/pull
                self.setPerturbationVoltage(
                        0, v)  # Adding Perturbation voltage on current iteration
                devpos = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                self.setPerturbationVoltage(0, -v)
                devmin = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                iMatKL[kl, :] = (devpos - devmin) / (2. * ampliVec[kl])
                #imat[:-2, :] /= pushDMMic
                #if(nmodesMax == 0):# i.e we measured all modes including TT
                #imat[-2:, :] /= pushTTArcsec
            else:  # No turbulence => push only
                self.openLoop()  # openLoop
                self.setPerturbationVoltage(0, v)
                iMatKL[kl, :] = self.applyVoltGetSlopes(noise=noise) / ampliVec[kl]
        self.setPerturbationVoltage(0, v * 0.)  # removing perturbvoltage...
        # print("Modal interaction matrix done in %3.0f seconds" % (time.time() - st))

        return iMatKL

    def doImatPhase(self, cubePhase, Nslopes, noise=False, nmodesMax=0, withTurbu=False,
                    pushPull=False, wfsnum=0):
        iMatPhase = np.zeros((cubePhase.shape[0], Nslopes))
        for nphase in trange(cubePhase.shape[0], desc="Phase IM"):
            if ((pushPull is True) or
                (withTurbu is True)):  # with turbulence/aberrations => push/pull
                self.setNcpaWfs(cubePhase[nphase, :, :], wfsnum=wfsnum)
                devpos = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                self.setNcpaWfs(-cubePhase[nphase, :, :], wfsnum=wfsnum)
                devmin = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                iMatPhase[nphase, :] = (devpos - devmin) / 2
            else:  # No turbulence => push only
                self.openLoop()  # openLoop
                self.setNcpaWfs(cubePhase[nphase, :, :], wfsnum=wfsnum)
                iMatPhase[nphase, :] = self.applyVoltGetSlopes(noise=noise)
        self.setNcpaWfs(cubePhase[nphase, :, :] * 0.,
                        wfsnum=wfsnum)  # Remove the Phase on WFS
        _ = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
        # print("Phase interaction matrix done in %3.0f seconds" % (time.time() - st))

        return iMatPhase

    def applyVoltGetSlopes(self, noise=False, turbu=False, reset=1):
        self._sim.rtc.apply_control(0, self._sim.dms)
        for w in range(len(self._sim.wfs.d_wfs)):

            if (turbu):
                self._sim.raytraceWfs(w, "all")
                # w.d_gs.raytrace(self._sim.atm)
                # w.d_gs.raytrace(self._sim.tel)
                # w.d_gs.raytrace(self._sim.dms)
                # w.d_gs.raytrace()
            else:
                self._sim.raytraceWfs(w, ["dm", "ncpa"], rst=reset)
                # w.d_gs.raytrace(self._sim.dms, rst=reset)
                # w.d_gs.raytrace()

            self._sim.compWfsImage(w, noise=noise)
        self._sim.rtc.do_centroids(0)
        c = self.getCentroids(0)
        return c

    def computeModalResiduals(self):
        self._sim.doControl(1, 0)
        v = self.getCom(
                1
        )  # We compute here the residual phase on the DM modes. Gives the Equivalent volts to apply/
        if (self.P is None):
            self.modalBasis, self.P = self.getModes2VBasis("Btt")
        if (self.selectedActus is None):
            ai = self.P.dot(v) * 1000.  # np rms units
        else:  # Slaving actus case
            v2 = v[:-2][list(
                    self.selectedActus)]  # If actus are slaved then we select them.
            v3 = v[-2:]
            ai = self.P.dot(np.concatenate((v2, v3))) * 1000.
        return ai

    def writeConfigOnFile(self,
                          filepath=os.environ["SHESHA_ROOT"] + "/widgets/canapass.conf"):
        aodict = OrderedDict()

        aodict.update({"Fe": 1 / self._sim.config.p_loop.ittime})
        aodict.update({"teldiam": self._sim.config.p_tel.diam})
        aodict.update({"telobs": self._sim.config.p_tel.cobs})

        # WFS
        aodict.update({"nbWfs": len(self._sim.config.p_wfss)})
        aodict.update({"nbTargets": len(self._sim.config.p_targets)})
        aodict.update({"nbCam": aodict["nbWfs"]})
        aodict.update({"nbOffaxis": 0})
        aodict.update({"nbNgsWFS": 1})
        aodict.update({"nbLgsWFS": 0})
        aodict.update({"nbFigSensor": 0})
        aodict.update({"nbSkyWfs": aodict["nbWfs"]})
        aodict.update({"nbOffNgs": 0})

        # DMS
        aodict.update({"nbDms": len(self._sim.config.p_dms)})
        aodict.update({"Nactu": self._sim.rtc.d_control[0].nactu})

        # List of things
        aodict.update({"list_NgsOffAxis": []})
        aodict.update({"list_Fig": []})
        aodict.update({"list_Cam": [0]})
        aodict.update({"list_SkyWfs": [0]})
        aodict.update({"list_ITS": []})
        aodict.update({"list_Woofer": []})
        aodict.update({"list_Tweeter": []})
        aodict.update({"list_Steering": []})

        # fct of Nb of wfss
        NslopesList = []
        NsubapList = []
        listWfsType = []
        pyrModulationList = []
        pyr_npts = []
        pyr_pupsep = []
        pixsize = []
        xPosList = []
        yPosList = []
        fstopsize = []
        fstoptype = []
        npixPerSub = []
        nxsubList = []
        nysubList = []
        lambdaList = []
        dms_seen = []
        colTmpList = []
        noise = []
        new_hduwfsl = pfits.HDUList()
        new_hduwfsSubapXY = pfits.HDUList()
        for i in range(aodict["nbWfs"]):
            new_hduwfsl.append(pfits.ImageHDU(
                    self._sim.config.p_wfss[i]._isvalid))  # Valid subap array
            new_hduwfsl[i].header["DATATYPE"] = "valid_wfs%d" % i

            xytab = np.zeros((2, self._sim.config.p_wfss[i]._validsubsx.shape[0]))
            xytab[0, :] = self._sim.config.p_wfss[i]._validsubsx
            xytab[1, :] = self._sim.config.p_wfss[i]._validsubsy

            new_hduwfsSubapXY.append(
                    pfits.ImageHDU(xytab))  # Valid subap array inXx Y on the detector
            new_hduwfsSubapXY[i].header["DATATYPE"] = "validXY_wfs%d" % i

            pixsize.append(self._sim.config.p_wfss[i].pixsize)
            NslopesList.append(self._sim.config.p_wfss[i]._nvalid * 2)  # slopes per wfs
            NsubapList.append(self._sim.config.p_wfss[i]._nvalid)  # subap per wfs
            listWfsType.append(self._sim.config.p_wfss[i].type)
            xPosList.append(self._sim.config.p_wfss[i].xpos)
            yPosList.append(self._sim.config.p_wfss[i].ypos)
            fstopsize.append(self._sim.config.p_wfss[i].fssize)
            fstoptype.append(self._sim.config.p_wfss[i].fstop)
            nxsubList.append(self._sim.config.p_wfss[i].nxsub)
            nysubList.append(self._sim.config.p_wfss[i].nxsub)
            lambdaList.append(self._sim.config.p_wfss[i].Lambda)
            dms_seen.append(list(self._sim.config.p_wfss[i].dms_seen))
            noise.append(self._sim.config.p_wfss[i].noise)

            if (self._sim.config.p_wfss[i].type == "pyrhr"):
                pyrModulationList.append(self._sim.config.p_wfss[i].pyr_ampl)
                pyr_npts.append(self._sim.config.p_wfss[i].pyr_npts)
                pyr_pupsep.append(self._sim.config.p_wfss[i].pyr_pup_sep)
                npixPerSub.append(1)
            else:
                pyrModulationList.append(0)
                pyr_npts.append(0)
                pyr_pupsep.append(0)
                npixPerSub.append(self._sim.config.p_wfss[i].npix)
        confname = filepath.split("/")[-1].split('.conf')[0]
        new_hduwfsl.writeto(
                filepath.split(".conf")[0] + '_wfsConfig.fits', overwrite=True)
        new_hduwfsSubapXY.writeto(
                filepath.split(".conf")[0] + '_wfsValidXYConfig.fits', overwrite=True)
        aodict.update({"listWFS_NslopesList": NslopesList})
        aodict.update({"listWFS_NsubapList": NsubapList})
        aodict.update({"listWFS_WfsType": listWfsType})
        aodict.update({"listWFS_pixarc": pixsize})
        aodict.update({"listWFS_pyrModRadius": pyrModulationList})
        aodict.update({"listWFS_pyrModNPts": pyr_npts})
        aodict.update({"listWFS_pyrPupSep": pyr_pupsep})
        aodict.update({"listWFS_fstopsize": fstopsize})
        aodict.update({"listWFS_fstoptype": fstoptype})
        aodict.update({"listWFS_dms_seen": dms_seen})
        aodict.update({"listWFS_NsubX": nxsubList})
        aodict.update({"listWFS_NsubY": nysubList})
        aodict.update({"listWFS_Nsub": nysubList})
        aodict.update({"listWFS_NpixPerSub": npixPerSub})
        aodict.update({"listWFS_Lambda": lambdaList})
        aodict.update({"listWFS_noise": noise})

        listDmsType = []
        NactuX = []
        unitPerVolt = []
        push4imat = []
        coupling = []
        push4iMatArcSec = []
        new_hdudmsl = pfits.HDUList()

        for j in range(aodict["nbDms"]):
            listDmsType.append(self._sim.config.p_dms[j].type)
            NactuX.append(self._sim.config.p_dms[j].nact)
            unitPerVolt.append(self._sim.config.p_dms[j].unitpervolt)
            push4imat.append(self._sim.config.p_dms[j].push4imat)
            coupling.append(self._sim.config.p_dms[j].coupling)
            tmp = []
            if (self._sim.config.p_dms[j].type != 'tt'):
                tmpdata = np.zeros((2, len(self._sim.config.p_dm0._i1)))
                tmpdata[0, :] = self._sim.config.p_dm0._j1
                tmpdata[1, :] = self._sim.config.p_dm0._i1
                new_hdudmsl.append(pfits.ImageHDU(tmpdata))  # Valid subap array
                new_hdudmsl[j].header["DATATYPE"] = "valid_dm%d" % j
            #for k in range(aodict["nbWfs"]):
            #    tmp.append(self._sim.computeDMrange(j, k))

            push4iMatArcSec.append(tmp)
        new_hdudmsl.writeto(
                filepath.split(".conf")[0] + '_dmsConfig.fits', overwrite=True)
        aodict.update({"listDMS_push4iMatArcSec": push4iMatArcSec})
        aodict.update({"listDMS_push4iMat": push4imat})
        aodict.update({"listDMS_unitPerVolt": unitPerVolt})
        aodict.update({"listDMS_Nxactu": NactuX})
        aodict.update({"listDMS_Nyactu": NactuX})
        aodict.update({"listDMS_type": listDmsType})
        aodict.update({"listDMS_coupling": coupling})

        listTargetsLambda = []
        listTargetsXpos = []
        listTargetsYpos = []
        listTargetsDmsSeen = []
        listTargetsMag = []
        for k in range(aodict["nbTargets"]):
            listTargetsLambda.append(self._sim.config.p_targets[k].Lambda)
            listTargetsXpos.append(self._sim.config.p_targets[k].xpos)
            listTargetsYpos.append(self._sim.config.p_targets[k].ypos)
            listTargetsMag.append(self._sim.config.p_targets[k].mag)
            listTargetsDmsSeen.append(list(self._sim.config.p_targets[k].dms_seen))

        aodict.update({"listTARGETS_Lambda": listTargetsLambda})
        aodict.update({"listTARGETS_Xpos": listTargetsXpos})
        aodict.update({"listTARGETS_Ypos": listTargetsYpos})
        aodict.update({"listTARGETS_Mag": listTargetsMag})
        aodict.update({"listTARGETS_DmsSeen": listTargetsDmsSeen})

        listDmsType = []
        Nslopes = sum(NslopesList)
        Nsubap = sum(NsubapList)
        aodict.update({"Nslopes": Nslopes})
        aodict.update({"Nsubap": Nsubap})
        f = open(filepath, 'w+')
        for dictval in aodict:
            f.write(dictval + ":" + str(aodict[dictval]) + "\n")
        f.close()
        print("OK: Config File wrote in:" + filepath)
        #return aodict

    def setPyrModulation(self, pyrmod):
        CompassSupervisor.setPyrModulation(self, pyrmod)
        self._sim.rtc.do_centroids(0)  # To be ready for the next getSlopes

    def setPyrMethod(self, pyrmethod):
        CompassSupervisor.setPyrMethod(self, pyrmethod)
        self._sim.rtc.do_centroids(0)  # To be ready for the next getSlopes

    # def setNoise(self, noise, numwfs=0):
    #     CompassSupervisor.setNoise(self, noise, numwfs)

    def getTargetPhase(self, tarnum):
        pup = self.getSpupil()
        ph = self.getTarPhase(tarnum) * pup
        return ph

    def getNcpaWfs(self, wfsnum):
        return np.array(self._sim.wfs.d_wfs[wfsnum].d_gs.d_ncpa_phase)

    def getNcpaTar(self, tarnum):
        return np.array(self._sim.tar.d_targets[tarnum].d_ncpa_phase)

    #def getVolts(self):
    #    return self._sim.rtc.get_voltage(0)
    #
    #def getSlopes(self):
    #     return self._sim.rtc.get_centroids(0)

    def writeDataInFits(self, data, fullpath):
        pfits.writeto(fullpath, data, overwrite=True)

    def getFrameCounter(self):
        return self._sim.iter

    def recordCB(self, CBcount, subSample=1, tarnum=0, seeAtmos=True,
                 tarPhaseFilePath="", NCPA=False, ncpawfs=None, refSlopes=None):
        slopesdata = None
        voltsdata = None
        tarPhaseData = None
        aiData = None
        k = 0
        srseList = []
        srleList = []
        gNPCAList = []
        # Resets the target so that the PSF LE is synchro with the data
        for i in range(len(self._sim.config.p_targets)):
            self.resetStrehl(i)

        # Starting CB loop...
        for j in trange(CBcount, desc="recording"):
            if (NCPA):
                if (j % NCPA == 0):
                    ncpaDiff = refSlopes[None, :]
                    ncpaturbu = self.doImatPhase(-ncpawfs[None, :, :],
                                                 refSlopes.shape[0], noise=False,
                                                 withTurbu=True)
                    gNCPA = float(
                            np.sqrt(
                                    np.dot(ncpaDiff, ncpaDiff.T) / np.dot(
                                            ncpaturbu, ncpaturbu.T)))
                    if (gNCPA > 1e18):
                        gNCPA = 0
                        print('Warning NCPA ref slopes gain too high!')
                        gNPCAList.append(gNCPA)
                        self.setRefSlopes(-refSlopes * gNCPA)
                    else:
                        gNPCAList.append(gNCPA)
                        print('NCPA ref slopes gain: %4.3f' % gNCPA)
                        self.setRefSlopes(-refSlopes / gNCPA)

            self._sim.next(see_atmos=seeAtmos)
            for t in range(len(self._sim.config.p_targets)):
                self._sim.compTarImage(t)

            srse, srle, _, _ = self.getStrehl(tarnum)
            srseList.append(srse)
            srleList.append(srle)
            if (j % subSample == 0):
                aiVector = self.computeModalResiduals()
                if (aiData is None):
                    aiData = np.zeros((len(aiVector), int(CBcount / subSample)))
                aiData[:, k] = aiVector

                slopesVector = self.getCentroids(0)
                if (slopesdata is None):
                    slopesdata = np.zeros((len(slopesVector), int(CBcount / subSample)))
                slopesdata[:, k] = slopesVector

                voltsVector = self.getCom(0)
                if (voltsdata is None):
                    voltsdata = np.zeros((len(voltsVector), int(CBcount / subSample)))
                voltsdata[:, k] = voltsVector

                if (tarPhaseFilePath != ""):
                    tarPhaseArray = self.getTargetPhase(tarnum)
                    if (tarPhaseData is None):
                        tarPhaseData = np.zeros((*tarPhaseArray.shape,
                                                 int(CBcount / subSample)))
                    tarPhaseData[:, :, k] = tarPhaseArray
                k += 1
        if (tarPhaseFilePath != ""):
            print("Saving tarPhase cube at: ", tarPhaseFilePath)
            pfits.writeto(tarPhaseFilePath, tarPhaseData, overwrite=True)
        psfLE = self.getTarImage(tarnum, "le")
        return slopesdata, voltsdata, aiData, psfLE, srseList, srleList, gNPCAList

        #wao.sim.config.p_geom._ipupil
        """
        ------------------------------------------------
        With a SH:
        ------------------------------------------------
        plt.clf()
        plt.matshow(wao.config.p_geom._ipupil, fignum=1)

        # DM positions in iPupil:
        dmposx = wao.config.p_dm0._xpos
        dmposy = wao.config.p_dm0._ypos
        plt.scatter(dmposy, dmposx, color="blue", label="DM actuators")
        #plt.scatter(dmposy[list(wao.supervisor.slavedActus)], dmposx[list(wao.supervisor.slavedActus)], color="orange", label="Slaved")

        #WFS position in ipupil
        ipup = wao.config.p_geom._ipupil
        spup = wao.config.p_geom._spupil
        s2ipup = (ipup.shape[0] - spup.shape[0]) / 2.
        posx = wao.config.p_wfss[0]._validsubsx * wao.config.p_wfs0._pdiam + s2ipup
        posy = wao.config.p_wfss[0]._validsubsy * wao.config.p_wfs0._pdiam + s2ipup

        #center of ssp position in ipupil
        demissp = wao.config.p_wfs0._pdiam / 2 - 0.5
        sspx = posx+demissp
        sspy = posy+demissp
        plt.scatter(sspy, sspx, color="red", label="WFS SSP")

        ------------------------------------------------
        With a PYRAMID:
        ------------------------------------------------

        plt.matshow(wao.config.p_geom._ipupil, fignum=1)
        size=10


        # DM positions in iPupil:
        dmposx = wao.config.p_dm0._xpos
        dmposy = wao.config.p_dm0._ypos
        plt.scatter(dmposy, dmposx,s=size, color="blue", label="DM actuators")
        #plt.scatter(dmposy[list(wao.supervisor.couplesActus)], dmposx[list(wao.supervisor.couplesActus)], color="orange", label="Slaved")

        #WFS position in ipupil
        ipup = wao.config.p_geom._ipupil
        spup = wao.config.p_geom._spupil
        s2ipup = (ipup.shape[0] - spup.shape[0]) / 2.
        posx = wao.config.p_wfss[0]._istart + s2ipup
        posx = np.tile(posx,(posx.size,1))
        posy = posx.T.copy()
        posx = posx *  wao.config.p_wfss[0]._isvalid
        posx = posx[np.where(posx > 0)] - ipup.shape[0] / 2 -wao.config.p_wfs0.npix
        posy = posy * wao.config.p_wfss[0]._isvalid
        posy = posy[np.where(posy > 0)] - ipup.shape[0] / 2 -wao.config.p_wfs0.npix

        #center of ssp position in ipupil
        demissp = wao.config.p_wfs0._pdiam / 2 -0.5
        sspx = posx+ipup.shape[0]/2+demissp
        sspy = posy+ipup.shape[0]/2+demissp
        plt.scatter(sspy, sspx, color="red", s=size, label="WFS SSP")




        """


if __name__ == '__main__':
    from docopt import docopt
    arguments = docopt(__doc__)
    supervisor = CanapassSupervisor(arguments["<parameters_filename>"], True)
    supervisor.initConfig()

    try:
        from subprocess import Popen, PIPE
        from hraa.server.pyroServer import PyroServer

        p = Popen("whoami", shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if (err != b''):
            print(err)
            raise ValueError("ERROR CANNOT RECOGNIZE USER")
        else:
            user = out.split(b"\n")[0].decode("utf-8")
            print("User is " + user)
        server = PyroServer()
        server.add_device(supervisor, "waoconfig_" + user)
        server.start()
    except:
        raise EnvironmentError(
                "Missing dependencies (code HRAA or Pyro4 or Dill Serializer)")
