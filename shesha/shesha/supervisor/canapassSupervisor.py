## @package   shesha.supervisor.canapassSupervisor
## @brief     Initialization and execution of a CANAPASS supervisor
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   4.3.2
## @date      2011/01/28
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
#  All rights reserved.
#  Distributed under GNU - LGPL
#
#  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
#  General Public License as published by the Free Software Foundation, either version 3 of the License,
#  or any later version.
#
#  COMPASS: End-to-end AO simulation tool using GPU acceleration
#  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
#
#  The final product includes a software package for simulating all the critical subcomponents of AO,
#  particularly in the context of the ELT and a real-time core based on several control approaches,
#  with performances consistent with its integration into an instrument. Taking advantage of the specific
#  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
#  conduct large simulation campaigns called to the ELT.
#
#  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
#  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
#  various systems configurations such as multi-conjugate AO.
#
#  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
#  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
#  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
"""
Initialization and execution of a CANAPASS supervisor

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

from tqdm import tqdm
import astropy.io.fits as pfits
from threading import Thread
from subprocess import Popen, PIPE

import shesha.ao as ao
import shesha.constants as scons
from shesha.constants import CentroiderType, WFSType

from typing import Any, Dict, Tuple, Callable, List
from .compassSupervisor import CompassSupervisor

# from carmaWrap.obj import obj_Double2D
# from carmaWrap.magma import syevd_Double, svd_host_Double
# from carmaWrap.context import context as carmaWrap_context

# from carmaWrap.host_obj import host_obj_Double1D, host_obj_Double2D


class CanapassSupervisor(CompassSupervisor):

    def __init__(self, configFile: str = None, cacao: bool = True) -> None:
        CompassSupervisor.__init__(self, configFile=configFile, cacao=cacao)

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

    def setDmCommand(self, numdm, volts):
        """
        Allows to by-pass the RTC for sending a command to the
        specified DM <numdm>.
        This command comes in addition to the RTC computation.
        It allows a direct access the DM without using the RTC.

        <numdm> : number of the DM
        <volts> : voltage vector to be applied on the DM.
        """
        ntotDm = len(self._sim.config.p_dms)
        if (numdm < ntotDm):
            self._sim.dms.d_dms[numdm].set_com(volts)
        else:
            print("ERROR !!!!\nRequested DM (", numdm,
                  ") conflicts with number of available DMs (", ntotDm, ").")

    def getConfig(self, path=None):
        ''' Returns the configuration in use, in a supervisor specific format '''
        if path:
            return self.writeConfigOnFile(path)
        else:
            return CompassSupervisor.getConfig(self)

    def loadConfig(self, configFile: str = None, sim=None) -> None:
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
        for i in range(nbmode):
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

    def computePh2ModesFits(self, fullpath):
        data = self.computePh2Modes()
        self.writeDataInFits(data, fullpath)

    """
    def getModes2VBasis(self, ModalBasisType, merged=False, nbpairs=None):
        if (ModalBasisType == "KL2V"):
            print("Computing KL2V basis...")
            self.modalBasis, _ = self.returnkl2V()
            self.modalBasis *= np.sign(self.modalBasis[0,:])[None,:]
            return self.modalBasis, 0
        elif (ModalBasisType == "Btt"):
            print("Computing Btt basis...")
            self.modalBasis, self.P = self.compute_Btt2(inv_method="cpu_svd",
                                                        merged=merged, nbpairs=nbpairs)
            self.modalBasis *= np.sign(self.modalBasis[0,:])[None,:]
            return self.modalBasis, self.P
    """

    def first_nonzero(self, arr, axis, invalid_val=-1):
        """
        find the first non zero element of an array.
        """
        mask = arr != 0
        return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)

    def getModes2VBasis(self, ModalBasisType, merged=False, nbpairs=None,
                        returnDelta=False):
        """
        Pos signifies the sign of the first non zero element of the eigen vector is forced to be +
        """

        if (ModalBasisType == "KL2V"):
            print("Computing KL2V basis...")
            self.modalBasis, _ = self.returnkl2V()
            fnz = self.first_nonzero(self.modalBasis, axis=0)
            # Computing the sign of the first non zero element
            #sig = np.sign(self.modalBasis[[fnz, np.arange(self.modalBasis.shape[1])]])
            sig = np.sign(self.modalBasis[tuple([
                    fnz, np.arange(self.modalBasis.shape[1])
            ])])  # pour remove le future warning!
            self.modalBasis *= sig[None, :]
            return self.modalBasis, 0
        elif (ModalBasisType == "Btt"):
            print("Computing Btt basis...")
            self.modalBasis, self.P = self.compute_Btt2(inv_method="cpu_svd",
                                                        merged=merged, nbpairs=nbpairs,
                                                        returnDelta=returnDelta)
            fnz = self.first_nonzero(self.modalBasis, axis=0)
            # Computing the sign of the first non zero element
            #sig = np.sign(self.modalBasis[[fnz, np.arange(self.modalBasis.shape[1])]])
            sig = np.sign(self.modalBasis[tuple([
                    fnz, np.arange(self.modalBasis.shape[1])
            ])])  # pour remove le future warning!
            self.modalBasis *= sig[None, :]
            return self.modalBasis, self.P
        elif (ModalBasisType == "Btt_petal"):
            print("Computing Btt with a Petal basis...")
            self.modalBasis, self.P = self.compute_btt_petal()
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

    def compute_Btt2(self, inv_method: str = "cpu_svd", merged=False, nbpairs=None,
                     returnDelta=False):

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
        #     c = carmaWrap_context()
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
            Btt = Btt2
            P = P2
        if (returnDelta):
            P = delta
        return Btt, P

    def compute_btt_petal(self):
        """
        Done
        """

        # Tip-tilt + piston + petal modes
        IF = self.getIFsparse(1)
        IFpetal = self.getIFdm(1)
        IFtt = self.getIFdm(2)

        n = IF.shape[0]  # number of points (pixels) over the pupil
        N = IF.shape[1]  # number of influence functions (nb of actuators)

        # Compute matrix delta (geometric covariance of actus)
        delta = IF.dot(IF.T).toarray() / N

        # Petal basis generation (orthogonal to global piston)
        nseg = IFpetal.toarray().shape[0]

        petal_modes = -1 / (nseg - 1) * np.ones((nseg, (nseg - 1)))
        petal_modes += nseg / (nseg - 1) * np.eye(nseg)[:, 0:(
                nseg - 1)]  # petal modes within the petal dm space
        tau_petal = np.dot(IF.toarray(), IFpetal.toarray().T).dot(petal_modes)

        Tp = np.concatenate((IFtt.toarray(), np.ones((1, N))),
                            axis=0)  # Matrice contenant Petal Basis + Tip/Tilt + Piston
        deltaT = IF.dot(Tp.T) / N

        # Tip tilt + petals projection on the pzt dm
        tau = np.concatenate((tau_petal, np.linalg.inv(delta).dot(deltaT)), axis=1)

        # Famille generatrice sans tip tilt ni pétales
        G = np.identity(n)
        tdt = tau.T.dot(delta).dot(tau)
        subTT = tau.dot(np.linalg.inv(tdt)).dot(tau.T).dot(delta)
        G -= subTT

        # Base orthonormee sans Tip, Tilp, Piston, Pétales
        gdg = G.T.dot(delta).dot(G)
        U, s, V = np.linalg.svd(gdg)
        U = U[:, :U.shape[1] - 8]
        s = s[:s.size - 8]
        L = np.identity(s.size) / np.sqrt(s)
        B = G.dot(U).dot(L)

        # Rajout du TT et Pétales
        TT = IFtt.toarray().dot(IFtt.toarray().T) / N  # .toarray()/N
        Btt = np.zeros((n + 2, n - 1))
        Btt[:n, :B.shape[1]] = B
        mini = 1. / np.sqrt(np.abs(TT))
        mini[0, 1] = 0
        mini[1, 0] = 0
        Btt[n:, -2:] = mini  # ajout du tip tilt sur le miroir tip tilt
        Btt[:n, -7:-2] = tau_petal  # ajout des modes pétales sur le miroir M4

        # Calcul du projecteur actus-->modes
        delta = np.zeros((n + IFtt.shape[0], n + IFtt.shape[0]))
        delta[:-2, :-2] = IF.dot(IF.T).toarray() / N
        delta[-2:, -2:] = IFtt.toarray().dot(IFtt.toarray().T) / N
        P = Btt.T.dot(delta)

        return Btt.astype(np.float32), P.astype(np.float32)

    def doImatModal(self, ampliVec, KL2V, Nslopes, noise=False, nmodesMax=0,
                    withTurbu=False, pushPull=False):
        """
        KL2V is the btt2V matrix
        """
        iMatKL = np.zeros((Nslopes, KL2V.shape[1]))

        if (nmodesMax):
            KLMax = nmodesMax
        else:
            KLMax = KL2V.shape[1]
        vold = self.getCom(0)
        self.openLoop(rst=False)
        for kl in range(KLMax):
            # v = ampliVec[kl] * KL2V[:, kl:kl + 1].T.copy()
            v = ampliVec[kl] * KL2V[:, kl]
            if ((pushPull is True) or
                (withTurbu is True)):  # with turbulence/aberrations => push/pull
                self.setPerturbationVoltage(
                        0, "imatModal",
                        vold + v)  # Adding Perturbation voltage on current iteration
                devpos = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                self.setPerturbationVoltage(0, "imatModal", vold - v)
                devmin = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                iMatKL[:, kl] = (devpos - devmin) / (2. * ampliVec[kl])
                #imat[:-2, :] /= pushDMMic
                #if(nmodesMax == 0):# i.e we measured all modes including TT
                #imat[-2:, :] /= pushTTArcsec
            else:  # No turbulence => push only
                self.openLoop()  # openLoop
                self.setPerturbationVoltage(0, "imatModal", v)
                iMatKL[:, kl] = self.applyVoltGetSlopes(noise=noise) / ampliVec[kl]
        self.removePerturbationVoltage(0, "imatModal")
        if ((pushPull is True) or (withTurbu is True)):
            self.closeLoop()  # We are supposed to be in close loop now
        return iMatKL

    def doImatPhase(self, cubePhase, Nslopes, noise=False, nmodesMax=0, withTurbu=False,
                    pushPull=False, wfsnum=0):
        iMatPhase = np.zeros((cubePhase.shape[0], Nslopes))
        for nphase in range(cubePhase.shape[0]):
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
        self._sim.rtc.apply_control(0)
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
        dataDict = {}
        root = self._sim

        if (root.config.p_tel is not None):
            aodict.update({"teldiam": root.config.p_tel.diam})
            aodict.update({"telobs": root.config.p_tel.cobs})
            # TURBU
            aodict.update({"r0": root.config.p_atmos.r0})
            aodict.update({"Fe": 1 / root.config.p_loop.ittime})

        # WFS
        aodict.update({"nbWfs": len(root.config.p_wfss)})
        aodict.update({"nbTargets": len(root.config.p_targets)})
        aodict.update({"nbCam": aodict["nbWfs"]})
        aodict.update({"nbOffaxis": 0})
        aodict.update({"nbNgsWFS": 1})
        aodict.update({"nbLgsWFS": 0})
        aodict.update({"nbFigSensor": 0})
        aodict.update({"nbSkyWfs": aodict["nbWfs"]})
        aodict.update({"nbOffNgs": 0})

        # DMS
        aodict.update({"nbDms": len(root.config.p_dms)})
        aodict.update({"Nactu": root.rtc.d_control[0].nactu})
        # List of things
        aodict.update({"list_NgsOffAxis": []})
        aodict.update({"list_Fig": []})
        aodict.update({"list_Cam": [0]})
        aodict.update({"list_SkyWfs": [0]})
        aodict.update({"list_ITS": []})
        aodict.update({"list_Woofer": []})
        aodict.update({"list_Tweeter": []})
        aodict.update({"list_Steering": []})

        listOfNstatesPerController = []
        listOfcontrolLawTypePerController = []
        for control in self.config.p_controllers:
            listOfNstatesPerController.append(control.nstates)
            listOfcontrolLawTypePerController.append(control.type)
        aodict.update({"list_nstatesPerController": listOfNstatesPerController})
        aodict.update({"list_controllerType": listOfcontrolLawTypePerController})

        # fct of Nb of wfss
        NslopesList = []
        NsubapList = []
        listWfsType = []
        listCentroType = []

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
                    root.config.p_wfss[i]._isvalid))  # Valid subap array
            new_hduwfsl[i].header["DATATYPE"] = "valid_wfs%d" % i
            dataDict["wfsValid_" + str(i)] = root.config.p_wfss[i]._isvalid

            xytab = np.zeros((2, root.config.p_wfss[i]._validsubsx.shape[0]))
            xytab[0, :] = root.config.p_wfss[i]._validsubsx
            xytab[1, :] = root.config.p_wfss[i]._validsubsy
            dataDict["wfsValidXY_" + str(i)] = xytab

            new_hduwfsSubapXY.append(
                    pfits.ImageHDU(xytab))  # Valid subap array inXx Y on the detector
            new_hduwfsSubapXY[i].header["DATATYPE"] = "validXY_wfs%d" % i
            pixsize.append(root.config.p_wfss[i].pixsize)
            """
            if (root.config.p_centroiders[i].type == "maskedpix"):
                factor = 4
            else:
                factor = 2
            NslopesList.append(
                    root.config.p_wfss[i]._nvalid * factor)  # slopes per wfs
            """
            listCentroType.append(
                    root.config.p_centroiders[i].
                    type)  # assumes that there is the same number of centroiders and wfs
            NsubapList.append(root.config.p_wfss[i]._nvalid)  # subap per wfs
            listWfsType.append(root.config.p_wfss[i].type)
            xPosList.append(root.config.p_wfss[i].xpos)
            yPosList.append(root.config.p_wfss[i].ypos)
            fstopsize.append(root.config.p_wfss[i].fssize)
            fstoptype.append(root.config.p_wfss[i].fstop)
            nxsubList.append(root.config.p_wfss[i].nxsub)
            nysubList.append(root.config.p_wfss[i].nxsub)
            lambdaList.append(root.config.p_wfss[i].Lambda)
            dms_seen.append(list(root.config.p_wfss[i].dms_seen))
            noise.append(root.config.p_wfss[i].noise)

            if (root.config.p_centroiders[i].type == CentroiderType.MASKEDPIX):
                NslopesList.append(root.config.p_wfss[i]._nvalid * 4)  # slopes per wfs
            else:
                NslopesList.append(root.config.p_wfss[i]._nvalid * 2)  # slopes per wfs

            if (root.config.p_wfss[i].type == "pyrhr"):
                pyrModulationList.append(root.config.p_wfss[i].pyr_ampl)
                pyr_npts.append(root.config.p_wfss[i].pyr_npts)
                pyr_pupsep.append(root.config.p_wfss[i].pyr_pup_sep)
                npixPerSub.append(1)
            else:
                pyrModulationList.append(0)
                pyr_npts.append(0)
                pyr_pupsep.append(0)
                npixPerSub.append(root.config.p_wfss[i].npix)
        confname = filepath.split("/")[-1].split('.conf')[0]
        new_hduwfsl.writeto(
                filepath.split(".conf")[0] + '_wfsConfig.fits', overwrite=True)
        new_hduwfsSubapXY.writeto(
                filepath.split(".conf")[0] + '_wfsValidXYConfig.fits', overwrite=True)

        aodict.update({"listWFS_NslopesList": NslopesList})
        aodict.update({"listWFS_NsubapList": NsubapList})
        aodict.update({"listWFS_CentroType": listCentroType})
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
        Nactu = []
        unitPerVolt = []
        push4imat = []
        coupling = []
        push4iMatArcSec = []
        new_hdudmsl = pfits.HDUList()

        for j in range(aodict["nbDms"]):
            listDmsType.append(root.config.p_dms[j].type)
            NactuX.append(
                    root.config.p_dms[j].nact)  # nb of actuators across the diameter !!
            Nactu.append(root.config.p_dms[j]._ntotact)  # nb of actuators in total
            unitPerVolt.append(root.config.p_dms[j].unitpervolt)
            push4imat.append(root.config.p_dms[j].push4imat)
            coupling.append(root.config.p_dms[j].coupling)
            tmp = []
            if (root.config.p_dms[j]._i1 is
                        not None):  # Simu Case where i1 j1 is known (simulated)
                if (root.config.p_dms[j].type != 'tt'):
                    tmpdata = np.zeros((4, len(root.config.p_dms[j]._i1)))
                    tmpdata[0, :] = root.config.p_dms[j]._j1
                    tmpdata[1, :] = root.config.p_dms[j]._i1
                    tmpdata[2, :] = root.config.p_dms[j]._xpos
                    tmpdata[3, :] = root.config.p_dms[j]._ypos
                else:
                    tmpdata = np.zeros((4, 2))

                dataDict["dmData" + str(j)] = tmpdata
                new_hdudmsl.append(pfits.ImageHDU(tmpdata))  # Valid subap array
                new_hdudmsl[j].header["DATATYPE"] = "valid_dm%d" % j
                #for k in range(aodict["nbWfs"]):
                #    tmp.append(root.computeDMrange(j, k))

                push4iMatArcSec.append(tmp)
        new_hdudmsl.writeto(
                filepath.split(".conf")[0] + '_dmsConfig.fits', overwrite=True)
        aodict.update({"listDMS_push4iMatArcSec": push4iMatArcSec})
        aodict.update({"listDMS_push4iMat": push4imat})
        aodict.update({"listDMS_unitPerVolt": unitPerVolt})
        aodict.update({"listDMS_Nxactu": NactuX})
        aodict.update({"listDMS_Nyactu": NactuX})
        aodict.update({"listDMS_Nactu": Nactu})

        aodict.update({"listDMS_type": listDmsType})
        aodict.update({"listDMS_coupling": coupling})

        if (root.config.p_targets is not None):  # simu case
            listTargetsLambda = []
            listTargetsXpos = []
            listTargetsYpos = []
            listTargetsDmsSeen = []
            listTargetsMag = []
            for k in range(aodict["nbTargets"]):
                listTargetsLambda.append(root.config.p_targets[k].Lambda)
                listTargetsXpos.append(root.config.p_targets[k].xpos)
                listTargetsYpos.append(root.config.p_targets[k].ypos)
                listTargetsMag.append(root.config.p_targets[k].mag)
                listTargetsDmsSeen.append(list(root.config.p_targets[k].dms_seen))

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
        return aodict, dataDict

    def setPyrModulation(self, pyrmod):
        CompassSupervisor.setPyrModulation(self, pyrmod)
        self._sim.rtc.do_centroids(0)  # To be ready for the next getSlopes

    def setPyrSourceArray(self, cx, cy, nwfs=0):
        """
        Sets the Pyramid source Array
        cx, cy must be in arcseconds units

        """
        pyr_npts = len(cx)
        wfs = self._sim.wfs
        pwfs = self._sim.config.p_wfss[nwfs]
        pwfs.set_pyr_npts(pyr_npts)
        pwfs.set_pyr_cx(cx)
        pwfs.set_pyr_cy(cy)
        wfs.d_wfs[nwfs].set_pyr_modulation(cx, cy, pyr_npts)

        # RTC scale units to be updated ????
        #scale = pwfs.Lambda * 1e-6 / p_tel.diam * ampli * 180. / np.pi * 3600.
        #rtc.d_centro[nwfs].set_scale(scale)

    def setPyrMultipleStarsSource(self, coords, weights=None, pyrmod=3., niters=None,
                                  nwfs=0):
        """
        Sets the Pyramid source Array with a multiple star system
        coords is a list of couples of length n, coordinates of the n stars in lambda/D
        pyrmod is the modulation of the pyramid in lambda/D
        niters is the number of iteration

        """
        if niters == None:
            perim = pyrmod * 2 * np.pi
            niters = int((perim // 4 + 1) * 4)
            print(niters)
        nstars = len(coords)
        pyr_npts = niters * nstars
        wfs = self._sim.wfs
        pwfs = self._sim.config.p_wfss[nwfs]
        ptel = self._sim.config.p_tel
        #Computes the positions of the stars during the modulation
        pyrsize = pwfs._Nfft
        pixsize = (np.pi * pwfs._qpixsize) / (3600 * 180)
        scale_circ = 2 * np.pi / pyrsize * \
            (pwfs.Lambda * 1e-6 / ptel.diam) / pixsize * pyrmod
        scale_pos = 2 * np.pi / pyrsize * \
            (pwfs.Lambda * 1e-6 / ptel.diam) / pixsize
        temp_cx = []
        temp_cy = []
        for k in coords:
            temp_cx.append(scale_circ * \
                np.sin((np.arange(niters)) * 2. * np.pi / niters) + \
                k[0] * scale_pos)
            temp_cy.append(scale_circ * \
                np.cos((np.arange(niters)) * 2. * np.pi / niters) + \
                k[1] * scale_pos)
        cx = np.concatenate(np.array(temp_cx))
        cy = np.concatenate(np.array(temp_cy))
        #Gives the arguments to the simulation
        pwfs.set_pyr_npts(pyr_npts)
        pwfs.set_pyr_cx(cx)
        pwfs.set_pyr_cy(cy)
        if weights == None:
            wfs.d_wfs[nwfs].set_pyr_modulation(cx, cy, pyr_npts)
        else:
            w = []
            for k in weights:
                w += niters * [k]
            wfs.d_wfs[nwfs].set_pyr_modulation(cx, cy, w, pyr_npts)

        # RTC scale units to be updated ????
        #scale = pwfs.Lambda * 1e-6 / p_tel.diam * ampli * 180. / np.pi * 3600.
        #rtc.d_centro[nwfs].set_scale(scale)

    def setPyrDiskSourceHP(self, radius, density=1., nwfs=0):
        """
        radius is the radius of the disk object in lambda/D
        density is the spacing between the packed PSF in the disk object, in lambda/D

        create disk object by packing PSF in a given radius, using hexagonal packing
        /!\ There is no modulation
        """
        wfs = self._sim.wfs
        pwfs = self._sim.config.p_wfss[nwfs]
        ptel = self._sim.config.p_tel
        pyrsize = pwfs._Nfft
        pixsize = (np.pi * pwfs._qpixsize) / (3600 * 180)
        scale_pos = 2 * np.pi / pyrsize * \
            (pwfs.Lambda * 1e-6 / ptel.diam) / pixsize
        #Vectors used to generate the hexagonal paving
        gen_xp, gen_yp = np.array([1,
                                   0.]), np.array([np.cos(np.pi / 3),
                                                   np.sin(np.pi / 3)])
        n = 1 + int(1.2 * radius)
        mat_circ = []
        for k in range(-n, n):
            for l in range(-n, n):
                coord = k * gen_xp + l * gen_yp
                if np.sqrt(coord[0]**2 + coord[1]**2) <= radius:
                    mat_circ.append(coord)
        mat_circ = np.array(mat_circ)
        cx, cy = mat_circ[:, 0], mat_circ[:, 1]
        pyr_npts = len(cx)
        pwfs.set_pyr_npts(pyr_npts)
        pwfs.set_pyr_cx(cx)
        pwfs.set_pyr_cy(cy)
        wfs.d_wfs[nwfs].set_pyr_modulation(cx, cy, pyr_npts)

    def setPyrDiskSourceSP(self, radius, density=1., nwfs=0):
        """
        radius is the radius of the disk object in lambda/D
        density is the spacing between the packed PSF in the disk object, in lambda/D

        create disk object by packing PSF in a given radius, using square packing
        /!\ There is no modulation
        """

        def generate_square_circ(radius, density=1.):
            x = np.linspace(-radius, radius, 1 + 2 * int(radius / density))
            cx, cy = np.meshgrid(x, x, indexing='ij')
            cx = cx.flatten()
            cy = cy.flatten()
            r = cx * cx + cy * cy <= radius**2
            return (cx[r], cy[r])

        wfs = self._sim.wfs
        pwfs = self._sim.config.p_wfss[nwfs]
        ptel = self._sim.config.p_tel
        pyrsize = pwfs._Nfft
        pixsize = (np.pi * pwfs._qpixsize) / (3600 * 180)
        scale_pos = 2 * np.pi / pyrsize * \
            (pwfs.Lambda * 1e-6 / ptel.diam) / pixsize

        cx, cy = generate_square_circ(radius, density)
        cx = cx.flatten() * scale_pos
        cy = cy.flatten() * scale_pos
        pyr_npts = len(cx)
        pwfs.set_pyr_npts(pyr_npts)
        pwfs.set_pyr_cx(cx)
        pwfs.set_pyr_cy(cy)
        wfs.d_wfs[nwfs].set_pyr_modulation(cx, cy, pyr_npts)

    def setPyrSquareSource(self, radius, density=1., nwfs=0):
        """
        radius is half of the side of the object in lambda/D
        density is the spacing between the packed PSF in the square object, in lambda/D

        create square object by packing PSF in a given radius, using square packing
        /!\ There is no modulation
        """
        wfs = self._sim.wfs
        pwfs = self._sim.config.p_wfss[nwfs]
        ptel = self._sim.config.p_tel
        pyrsize = pwfs._Nfft
        pixsize = (np.pi * pwfs._qpixsize) / (3600 * 180)
        scale_pos = 2 * np.pi / pyrsize * \
            (pwfs.Lambda * 1e-6 / ptel.diam) / pixsize
        x = np.linspace(-radius, radius, 1 + 2 * int(radius / density)) * scale_pos
        cx, cy = np.meshgrid(x, x, indexing='ij')
        cx = cx.flatten()
        cy = cy.flatten()
        pyr_npts = len(cx)
        pwfs.set_pyr_npts(pyr_npts)
        pwfs.set_pyr_cx(cx)
        pwfs.set_pyr_cy(cy)
        wfs.d_wfs[nwfs].set_pyr_modulation(cx, cy, pyr_npts)

    def setPyrPseudoSource(self, radius, additional_psf=0, density=1., nwfs=0):

        def generate_square(radius, density=1.):
            """
            radius is half the length of a side in lambda/D
            density is the number of psf per lambda/D
            """
            x = np.linspace(-radius, radius, 1 + 2 * int(radius / density))
            cx, cy = np.meshgrid(x, x, indexing='ij')
            cx = cx.flatten()
            cy = cy.flatten()
            return (cx, cy)

        def generate_square_circ(radius, density=1.):
            x = np.linspace(-radius, radius, 1 + 2 * int(radius / density))
            cx, cy = np.meshgrid(x, x, indexing='ij')
            cx = cx.flatten()
            cy = cy.flatten()
            r = cx * cx + cy * cy <= radius**2
            return (cx[r], cy[r])

        def generate_pseudo_source(radius, additional_psf=0, density=1.):
            struct_size = (1 + 2 * additional_psf)**2
            center_x, center_y = generate_square(additional_psf, density)
            center_weight = (1 + 2 * int(additional_psf / density))**2 * [1]
            center_size = 1 + 2 * int(additional_psf / density)

            weight_edge = [(1 + 2 * int(radius / density) - center_size) // 2]
            xc, yc = generate_square_circ(radius, density)
            for k in range(additional_psf):
                line_length = np.sum(yc == (k + 1))
                print(line_length)
                weight_edge.append((line_length - center_size) // 2)

            edge_dist = (radius + additional_psf) // 2
            V_edge_x = []
            V_edge_y = []
            V_edge_weight = []
            for m in [-1, 1]:
                V_edge_x.append(0)
                V_edge_y.append(m * edge_dist)
                V_edge_weight.append(weight_edge[0])
            for k, val in enumerate(weight_edge[1:]):
                for l in [-1, 1]:
                    for m in [-1, 1]:
                        V_edge_x.append(l * (k + 1) * density)
                        V_edge_y.append(m * edge_dist)
                        V_edge_weight.append(val)
            H_edge_x = []
            H_edge_y = []
            H_edge_weight = []
            for m in [-1, 1]:
                H_edge_x.append(m * edge_dist)
                H_edge_y.append(0)
                H_edge_weight.append(weight_edge[0])
            for k, val in enumerate(weight_edge[1:]):
                for l in [-1, 1]:
                    for m in [-1, 1]:
                        H_edge_x.append(m * edge_dist)
                        H_edge_y.append(l * (k + 1) * density)
                        H_edge_weight.append(val)
            pup_cent_x = []
            pup_cent_y = []
            pup_cent_weight = 4 * [
                    (len(xc) - 2 * np.sum(H_edge_weight) - struct_size) / 4
            ]
            pup_cent_dist = int(edge_dist // np.sqrt(2))
            for l in [-1, 1]:
                for m in [-1, 1]:
                    pup_cent_x.append(l * pup_cent_dist)
                    pup_cent_y.append(m * pup_cent_dist)
            ox = np.concatenate((center_x, V_edge_x, H_edge_x, pup_cent_x))
            oy = np.concatenate((center_y, V_edge_y, H_edge_y, pup_cent_y))
            w = np.concatenate((center_weight, V_edge_weight, H_edge_weight,
                                pup_cent_weight))
            return (ox, oy, w, xc, yc)

        cx, cy, w, _, _ = generate_pseudo_source(radius, additional_psf, density)

        wfs = self._sim.wfs
        pwfs = self._sim.config.p_wfss[nwfs]
        ptel = self._sim.config.p_tel
        pyrsize = pwfs._Nfft
        pixsize = (np.pi * pwfs._qpixsize) / (3600 * 180)
        scale_pos = 2 * np.pi / pyrsize * \
            (pwfs.Lambda * 1e-6 / ptel.diam) / pixsize

        cx = cx.flatten() * scale_pos
        cy = cy.flatten() * scale_pos
        pyr_npts = len(cx)
        pwfs.set_pyr_npts(pyr_npts)
        pwfs.set_pyr_cx(cx)
        pwfs.set_pyr_cy(cy)
        wfs.d_wfs[nwfs].set_pyr_modulation(cx, cy, w, pyr_npts)

    def setPyrMethod(self, pyrmethod):
        CompassSupervisor.setPyrMethod(self, pyrmethod)
        self._sim.rtc.do_centroids(0)  # To be ready for the next getSlopes

    def setGain(self, gain) -> None:
        '''
        Set the scalar gain of feedback controller loop
        '''
        print("canapass")
        if type(gain) in [int, float]:
            self._sim.rtc.d_control[0].set_gain(gain)
        else:
            raise ValueError(
                    "ERROR CANNOT set array gain in canapass (generic + integrator law")


########################## PROTO #############################

    def initModalGain(self, gain, cmatModal, modalBasis, control=0, resetGain=True):
        """
        Given a gain, cmat and btt2v initialise the modal gain mode
        """
        print("TODO: A RECODER !!!!")
        nmode_total = modalBasis.shape[1]
        nactu_total = modalBasis.shape[0]
        nfilt = nmode_total - cmatModal.shape[0]
        ctrl = self._sim.rtc.d_control[control]
        ctrl.set_commandlaw('modal_integrator')
        cmat = np.zeros((nactu_total, cmatModal.shape[1]))
        dec = cmat.shape[0] - cmatModal.shape[0]
        cmat[:-dec, :] += cmatModal  # Fill the full Modal with all non-filtered modes
        modes2V = np.zeros((nactu_total, nactu_total))
        dec2 = modes2V.shape[1] - modalBasis.shape[1]
        modes2V[:, :-dec2] += modalBasis
        mgain = np.ones(len(modes2V)) * gain  # Initialize the gain
        ctrl.set_matE(modes2V)
        ctrl.set_cmat(cmat)
        if resetGain:
            ctrl.set_mgain(mgain)

    def leaveModalGain(self, control=0):
        ctrl = self._sim.rtc.d_control[control]
        ctrl.set_commandlaw('integrator')

    def set_mgain(self, mgain, control=0):
        """
        Wrapper function to let adopt set the modal gain.
        """
        ctrl = self._sim.rtc.d_control[control]
        ctrl.set_mgain(mgain)

    def get_mgain(self, control=0):
        """
        Wrapper function to let adopt get the modal gain.
        """
        ctrl = self._sim.rtc.d_control[control]
        return (np.array(ctrl.d_gain))

    def setModalBasis(self, modalBasis, P):
        """
        Function used to set the modal basis and projector in canapass
        """
        self.modalBasis = modalBasis
        self.P = P

    def getTargetPhase(self, tarnum):
        """
        Returns the target phase
        """
        pup = self.getSpupil()
        ph = self.getTarPhase(tarnum) * pup
        return ph

    def getInfluFunction(self, ndm):
        """
        returns the influence function cube for the given dm

        """
        return self._sim.config.p_dms[ndm]._influ

    def getInfluFunctionIpupilCoords(self, ndm):
        """
        returns the lower left coordinates of the influ function support in the ipupil coord system

        """
        i1 = self._sim.config.p_dm0._i1  # i1 is in the dmshape support coords
        j1 = self._sim.config.p_dm0._j1  # j1 is in the dmshape support coords
        ii1 = i1 + self._sim.config.p_dm0._n1  # in  ipupil coords
        jj1 = j1 + self._sim.config.p_dm0._n1  # in  ipupil coords
        return ii1, jj1

    #def getVolts(self):
    #    return self._sim.rtc.get_voltage(0)
    #
    #def getSlopes(self):
    #     return self._sim.rtc.get_centroids(0)

    def writeDataInFits(self, data, fullpath):
        pfits.writeto(fullpath, data, overwrite=True)

    def recordCB(self, CBcount, subSample=1, tarnum=0, seeAtmos=True,
                 tarPhaseFilePath="", NCPA=False, ncpawfs=None, refSlopes=None,
                 ditchStrehl=True):
        slopesdata = None
        voltsdata = None
        tarPhaseData = None
        aiData = None
        k = 0
        srseList = []
        srleList = []
        gNPCAList = []

        # Resets the target so that the PSF LE is synchro with the data
        # Doesn't reset it if DitchStrehl == False (used for real time gain computation)
        if ditchStrehl:
            for i in range(len(self._sim.config.p_targets)):
                self.resetStrehl(i)

        # Starting CB loop...
        for j in range(CBcount):
            print(j, end="\r")
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

    def recordCB_PyrFocPla_PSFSE(self, CBcount, starting_index, save_path,
                                 image_subSample=20, subSample=1, tarnum=0,
                                 seeAtmos=True, tarPhaseFilePath="", NCPA=False,
                                 ncpawfs=None, refSlopes=None, ditchStrehl=True):
        slopesdata = None
        voltsdata = None
        tarPhaseData = None
        aiData = None
        k = 0
        index_counter = starting_index
        srseList = []
        srleList = []
        gNPCAList = []

        # Resets the target so that the PSF LE is synchro with the data
        # Doesn't reset it if DitchStrehl == False (used for real time gain computation)
        if ditchStrehl:
            for i in range(len(self._sim.config.p_targets)):
                self.resetStrehl(i)

        # Starting CB loop...
        for j in range(CBcount):
            print(j, end="\r")
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
            if (index_counter % image_subSample == 0):
                psfSE = self.getTarImage(0, "se")
                pyrFP = self.getPyrFocalPlane()
                pyrHR = self.getWfsImage()
                wfs = self.getWfsPhase(0)
                pfits.writeto(save_path + 'psfSE_' + str(index_counter).zfill(5), psfSE,
                              overwrite=True)
                pfits.writeto(save_path + 'pyrFP_' + str(index_counter).zfill(5), pyrFP,
                              overwrite=True)
                pfits.writeto(save_path + 'pyrHR_' + str(index_counter).zfill(5), pyrHR,
                              overwrite=True)
                pfits.writeto(save_path + 'wfs_' + str(index_counter).zfill(5), wfs,
                              overwrite=True)
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
            index_counter += 1
        if (tarPhaseFilePath != ""):
            print("Saving tarPhase cube at: ", tarPhaseFilePath)
            pfits.writeto(tarPhaseFilePath, tarPhaseData, overwrite=True)
        psfLE = self.getTarImage(tarnum, "le")
        return slopesdata, voltsdata, aiData, psfLE, srseList, srleList, gNPCAList, index_counter

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
        posx = wao.config.p_wfss[0]._validpuppixx + s2ipup
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
    supervisor = CanapassSupervisor(arguments["<parameters_filename>"], cacao=True)
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
