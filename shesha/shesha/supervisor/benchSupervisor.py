## @package   shesha.supervisor.benchSupervisor
## @brief     Initialization and execution of a Bench supervisor
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   4.4.0
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

import numpy as np

from shesha.constants import CentroiderType, WFSType
from shesha.init.dm_init import dm_init_standalone
from shesha.init.rtc_init import rtc_standalone
from shesha.sutra_wrap import carmaWrap_context

from .aoSupervisor import AoSupervisor

from typing import Callable


class BenchSupervisor(AoSupervisor):

    def __init__(self, configFile: str = None, BRAHMA: bool = False,
                 CACAO: bool = False):
        '''
        Init the COMPASS wih the configFile
        '''
        self.pauseLoop = None
        self.rtc = None
        self.frame = None
        self.BRAHMA = BRAHMA
        self.CACAO = CACAO
        self.iter = 0
        self.slopesIdx = None

        if configFile is not None:
            self.loadConfig(configFile=configFile)

    #     _    _         _                  _
    #    / \  | |__  ___| |_ _ __ __ _  ___| |_
    #   / _ \ | '_ \/ __| __| '__/ _` |/ __| __|
    #  / ___ \| |_) \__ \ |_| | | (_| | (__| |_
    # /_/   \_\_.__/|___/\__|_|  \__,_|\___|\__|
    #
    #  __  __      _   _               _
    # |  \/  | ___| |_| |__   ___   __| |___
    # | |\/| |/ _ \ __| '_ \ / _ \ / _` / __|
    # | |  | |  __/ |_| | | | (_) | (_| \__ \
    # |_|  |_|\___|\__|_| |_|\___/ \__,_|___/

    def singleNext(self) -> None:
        '''
        Move atmos -> getSlope -> applyControl ; One integrator step
        '''
        self.loadNewWfsFrame()
        if (self.pauseLoop is not True):
            self.computeWfsFrame()
            self.setCommand(0, np.array(self.rtc.d_control[0].d_voltage))
        if self.BRAHMA or self.CACAO:
            self.rtc.publish()
        self.iter += 1

    def getTarImage(self, tarID, expoType: str = "se") -> np.ndarray:
        '''
        Get an image from a target
        '''
        raise NotImplementedError("Not implemented")

    def setCommand(self, nctrl: int, command: np.ndarray) -> None:
        ''' TODO
        Immediately sets provided command to DMs - does not affect integrator
        '''
        # Do stuff
        self.dmSetCallback(command)
        # Btw, update the RTC state with the information
        # self.rtc.d_control[nctrl].set_com(command, command.size)

    def getCom(self, nControl: int = 0) -> np.ndarray:
        '''
        Get command from DM, and set it back to nCtrl controller.
        These should be equivalent, unless an external source controls the DM as well
        '''
        # Do something
        command = self.dmGetCallback()
        # Btw, update the RTC state with the information
        # self.rtc.d_control[nControl].set_com(command, command.size)

        return command

    #  ____                  _ _   _        __  __      _   _               _
    # / ___| _ __   ___  ___(_) |_(_) ___  |  \/  | ___| |_| |__   ___   __| |___
    # \___ \| '_ \ / _ \/ __| | __| |/ __| | |\/| |/ _ \ __| '_ \ / _ \ / _` / __|
    #  ___) | |_) |  __/ (__| | |_| | (__  | |  | |  __/ |_| | | | (_) | (_| \__ \
    # |____/| .__/ \___|\___|_|\__|_|\___| |_|  |_|\___|\__|_| |_|\___/ \__,_|___/
    #       |_|

    def __repr__(self):

        s = '--- BenchSupervisor ---\nRTC: ' + repr(self.rtc)
        if hasattr(self, '_cam'):
            s += '\nCAM: ' + repr(self._cam)
        if hasattr(self, '_dm'):
            s += '\nDM: ' + repr(self._dm)
        return s

    def loadNewWfsFrame(self, numWFS: int = 0) -> None:
        '''
            Acquire a new WFS frame and load it.
        '''
        self.frame = self.camCallback()
        if (type(self.frame) is tuple):
            numWFS = len(self.frame)
            for i in range(numWFS):
                self.rtc.d_centro[i].load_img(self.frame[i], self.frame[i].shape[0],
                                              self.frame[i].shape[1])
        else:
            self.rtc.d_centro[numWFS].load_img(self.frame, self.frame.shape[0],
                                               self.frame.shape[1])

    def computeWfsFrame(self):
        '''
            Compute the WFS frame: calibrate, centroid, commands.
        '''
        # for index, centro in enumerate(self.rtc.d_centro):
        for centro in self.rtc.d_centro:
            centro.calibrate_img()
        self.rtc.do_centroids(0)
        self.rtc.do_control(0)
        self.rtc.do_clipping(0)
        self.rtc.comp_voltage(0)

    def setOneActu(self, nctrl: int, ndm: int, nactu: int, ampli: float = 1,
                   reset: bool = True) -> None:
        '''
        Push the selected actuator
        '''
        command = self.getCommand()
        if reset:
            command *= 0
        command[nactu] = ampli
        self.setCommand(nctrl, command)

    def forceContext(self) -> None:
        """
        Active all the GPU devices specified in the parameters file
        Required for using with widgets, due to multithreaded init
        and in case GPU 0 is not used by the simu
        """
        if self.isInit() and self.c is not None:
            current_Id = self.c.activeDevice
            for devIdx in range(len(self.config.p_loop.devices)):
                self.c.set_activeDeviceForce(devIdx)
            self.c.set_activeDevice(current_Id)

    def resetDM(self, nDM: int) -> None:
        '''
        Reset the DM number nDM
        '''
        if hasattr(self, '_dm'):
            self._dm.reset_dm()

    def resetCommand(self, nctrl: int = -1) -> None:
        '''
        Reset the nctrl Controller command buffer, reset all controllers if nctrl  == -1
        '''
        if (nctrl == -1):  # All Dms reset
            for control in self.rtc.d_control:
                control.d_com.reset()
        else:
            self.rtc.d_control[nctrl].d_com.reset()

    def loadConfig(self, configFile: str = None, sim=None) -> None:
        '''
        Init the COMPASS wih the configFile
        '''
        from shesha.util.utilities import load_config_from_file
        load_config_from_file(self, configFile)

    def setCamCallback(self, camCallback: Callable):
        '''
        Set the externally defined function that allows to grab frames
        '''
        self.camCallback = camCallback

    def setDmCallback(self, dmGetCallback: Callable, dmSetCallback: Callable):
        '''
        Set the externally defined function that allows to grab frames
        '''
        self.dmGetCallback = dmGetCallback
        self.dmSetCallback = dmSetCallback

    def isInit(self) -> bool:
        '''
        return the status on COMPASS init
        '''
        return self.is_init

    # TEST J
    def initConfig(self) -> None:
        '''
        Initialize the bench
        '''
        print("->RTC")
        self.wfsNb = len(self.config.p_wfss)
        print("Configuration of", self.wfsNb, "wfs ...")

        if (hasattr(self.config, 'p_loop') and self.config.p_loop.devices.size > 1):
            self.c = carmaWrap_context.get_instance_ngpu(self.config.p_loop.devices.size,
                                                         self.config.p_loop.devices)
        else:
            self.c = carmaWrap_context.get_instance_1gpu(self.config.p_loop.devices[0])
        nact = self.config.p_controllers[0].nactu
        self._nvalid = []
        self._centroiderType = []
        self._delay = []
        self._offset = []
        self._scale = []
        self._gain = []
        self._cMatSize = []
        self._npix = []

        # Get parameters
        for wfs in range(self.wfsNb):

            if self.config.p_wfss[wfs].type == WFSType.SH:
                self._npix.append(self.config.p_wfss[wfs].npix)
                if self.config.p_wfss[wfs]._validsubsx is None or \
                        self.config.p_wfss[wfs]._validsubsy is None:

                    from hraa.tools.doit import makessp
                    roiTab = makessp(self.config.p_wfss[wfs].nxsub, obs=0., rmax=0.98)
                    self.config.p_wfss[wfs]._nvalid = roiTab[0].size
                    self.config.p_wfss[
                            wfs]._validsubsx = roiTab[0] * self.config.p_wfss[wfs].npix
                    self.config.p_wfss[
                            wfs]._validsubsy = roiTab[1] * self.config.p_wfss[wfs].npix
                else:
                    self.config.p_wfss[wfs]._nvalid = self.config.p_wfss[
                            wfs]._validsubsx.size

                self._nvalid.append(
                        np.array([self.config.p_wfss[wfs]._nvalid], dtype=np.int32))
                # print("nvalid : %d" % self._nvalid[wfs])
                self._centroiderType.append(self.config.p_centroiders[wfs].type)
                self._delay.append(self.config.p_controllers[0].delay)  # ???
                self._offset.append((self.config.p_wfss[wfs].npix - 1) / 2)
                self._scale.append(1)
                self._gain.append(1)
                self._cMatSize.append(2 * self._nvalid[wfs][0])

            elif self.config.p_wfss[wfs].type == WFSType.PYRHR or self.config.p_wfss[
                    wfs].type == WFSType.PYRLR:
                self._nvalid.append(
                        np.array([self.config.p_wfss[wfs]._nvalid],
                                 dtype=np.int32))  # Number of valid SUBAPERTURES
                self._centroiderType.append(self.config.p_centroiders[wfs].type)
                self._delay.append(self.config.p_controllers[0].delay)  # ???
                self._offset.append(0)
                self._scale.append(1)
                self._gain.append(1)
                self._cMatSize.append(
                        self.config.p_wfss[wfs].nPupils * self._nvalid[wfs][0])
                self._npix.append(0)
            else:
                raise ValueError('WFS type not supported')

        # Create RTC
        self.rtc = rtc_standalone(self.c, self.wfsNb, self._nvalid, nact,
                                  self._centroiderType, self._delay, self._offset,
                                  self._scale, brahma=self.BRAHMA, cacao=self.CACAO)

        self.slopesIdx = np.cumsum([0] + [wfs.nslopes for wfs in self.rtc.d_centro])

        # Create centroiders
        for wfs in range(self.wfsNb):
            self.rtc.d_centro[wfs].load_validpos(
                    self.config.p_wfss[wfs]._validsubsx,
                    self.config.p_wfss[wfs]._validsubsy,
                    self.config.p_wfss[wfs]._validsubsx.size)
            if self.config.p_centroiders[wfs].type is CentroiderType.BPCOG:
                self.rtc.d_centro[wfs].set_nmax(self.config.p_centroiders[wfs].nmax)
            self.rtc.d_centro[wfs].set_npix(self._npix[wfs])
            # finally
            self.config.p_centroiders[wfs]._nslope = self.rtc.d_centro[wfs].nslopes
            print("wfs ", wfs, " set as ", self._centroiderType[wfs])
        size = sum(self._cMatSize)
        cMat = np.zeros((nact, size), dtype=np.float32)
        print("Size of cMat:", cMat.shape)

        # Initiate RTC
        self.rtc.d_control[0].set_cmat(cMat)
        self.rtc.d_control[0].set_decayFactor(
                np.ones(nact, dtype=np.float32) * (self._gain[0] - 1))
        self.rtc.d_control[0].set_matE(np.identity(nact, dtype=np.float32))
        self.rtc.d_control[0].set_mgain(np.ones(nact, dtype=np.float32) * -self._gain[0])
        self.is_init = True
        print("RTC initialized")

    def adaptiveWindows(self, initConfig=False, numwfs=0):
        '''
        Re-centre the centroiding boxes around the spots, and loads
        the new box coordinates in the slopes computation supervisor
        pipeline.

        Input arguments:
            <initConfig> : True/False(default).
            True resets to the default positions of boxes.
            False does the centring job normally.
        '''
        if initConfig:
            # reset de la configuration initiale
            ijSubap = self.config.p_wfss[numwfs].get_validsub()
            nsubap = ijSubap.shape[1]
            self.rtc.d_centro[numwfs].load_validpos(ijSubap[0], ijSubap[1], nsubap)
        else:
            # acquire slopes first
            nslopes = 10
            s = 0.
            for i in range(nslopes):
                self.loadNewWfsFrame()  # sinon toutes les slopes sont les memes
                self.computeWfsFrame()
                s = s + self.getSlope()[self.slopesIdx[numwfs]:self.slopesIdx[numwfs +
                                                                              1]]
            s /= nslopes
            # get coordinates of valid sub-apertures
            #ijSubap = self.config.p_wfss[numwfs].get_validsub()
            iSubap = np.array(self.rtc.d_centro[numwfs].d_validx)
            jSubap = np.array(self.rtc.d_centro[numwfs].d_validy)
            # get number of subaps
            nsubap = iSubap.shape[0]
            # reshape the array <s> to be conformable with <ijSubap>
            s = np.resize(s, (2, nsubap))
            # re-centre the boxes around the spots
            new_iSubap = (iSubap + s[0, :].round()).astype(int)
            new_jSubap = (jSubap + s[1, :].round()).astype(int)
            # load the new positions of boxes
            self.rtc.d_centro[numwfs].load_validpos(new_iSubap, new_jSubap, nsubap)

    def getCurrentWindowsPos(self, numwfs=0):
        """
        Returns the currently used subapertures positions.

        """
        iSubap = np.array(self.rtc.d_centro[numwfs].d_validx)
        jSubap = np.array(self.rtc.d_centro[numwfs].d_validy)
        return iSubap, jSubap

    def getSlopesIndex(self):
        return self.slopesIdx
