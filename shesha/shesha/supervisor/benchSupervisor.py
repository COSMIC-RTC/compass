## @package   shesha.supervisor.benchSupervisor
## @brief     Initialization and execution of a Bench supervisor
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   4.3.0
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
        self.computeWfsFrame()
        self.setCommand(0, np.array(self.rtc.d_control[0].d_voltage))
        if self.BRAHMA or self.CACAO:
            self.rtc.publish()

    def getTarImage(self, tarID, expoType: str = "se") -> np.ndarray:
        '''
        Get an image from a target
        '''
        raise NotImplementedError("Not implemented")

    def getWfsImage(self, numWFS: int = 0, calPix=False) -> np.ndarray:
        '''
        Get an image from the WFS
        '''
        if (calPix):
            return np.array(self.rtc.d_centro[0].d_img)
        else:
            return np.array(self.rtc.d_centro[0].d_img_raw)

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

    def __init__(self, configFile: str = None, BRAHMA: bool = False,
                 CACAO: bool = False):
        '''
        Init the COMPASS wih the configFile
        '''

        self.rtc = None
        self.frame = None
        self.BRAHMA = BRAHMA
        self.CACAO = CACAO

        if configFile is not None:
            self.loadConfig(configFile=configFile)

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
        self.rtc.d_centro[0].load_img(self.frame, self.frame.shape[0],
                                      self.frame.shape[1])

    def computeWfsFrame(self):
        '''
            Compute the WFS frame: calibrate, centroid, commands.
        '''
        self.rtc.d_centro[0].calibrate_img()
        self.rtc.do_centroids(0)
        # self.rtc.do_control(0)
        # self.rtc.do_clipping(0)
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

    def initConfig(self) -> None:
        '''
        Initialize the bench
        '''
        print("->RTC")
        wfsNb = len(self.config.p_wfss)
        p_wfs = self.config.p_wfss[0]
        if wfsNb > 1:
            raise RuntimeError("multi WFS not supported")

        if (hasattr(self.config, 'p_loop') and self.config.p_loop.devices.size > 1):
            self.c = carmaWrap_context.get_instance_ngpu(self.config.p_loop.devices.size,
                                                         self.config.p_loop.devices)
        else:
            self.c = carmaWrap_context.get_instance_1gpu(self.config.p_loop.devices[0])

        if p_wfs.type == WFSType.SH:
            self.npix = p_wfs.npix
            if p_wfs._validsubsx is None or \
                    p_wfs._validsubsy is None:
                # import rtcData.DataInit as di
                # dataS = di.makeSH(wfsNb=wfsNb, frameSize=self.cam.getWidth(),
                #                   roiSize=p_wfs.nxsub, subSize=self.npix)
                # p_wfs._nvalid = dataS.data["roiTab"].data.shape[1]
                # p_wfs._validsubsx = dataS.data["roiTab"].data[0, :]
                # p_wfs._validsubsy = dataS.data["roiTab"].data[1, :]

                from hraa.tools.doit import makessp
                roiTab = makessp(p_wfs.nxsub, obs=0., rmax=0.98)
                # for pos in self.roiTab: pos *= self.pitch
                p_wfs._nvalid = roiTab[0].size
                p_wfs._validsubsx = roiTab[0] * p_wfs.npix
                p_wfs._validsubsy = roiTab[1] * p_wfs.npix
            else:
                p_wfs._nvalid = p_wfs._validsubsx.size

            nvalid = np.array([p_wfs._nvalid], dtype=np.int32)
            print("nvalid : %d" % nvalid)
            offset = (p_wfs.npix - 1) / 2
            scale = 1
            gain = 1
            nact = self.config.p_dms[0].nact

            self.rtc = rtc_standalone(self.c, wfsNb, nvalid, nact,
                                      self.config.p_centroiders[0].type,
                                      self.config.p_controllers[0].delay, offset, scale,
                                      brahma=self.BRAHMA, cacao=self.CACAO)
            # put pixels in the SH grid coordonates
            self.rtc.d_centro[0].load_validpos(p_wfs._validsubsx, p_wfs._validsubsy,
                                               nvalid)
            if self.config.p_centroiders[0].type is CentroiderType.BPCOG:
                self.rtc.d_centro[0].set_nmax(self.config.p_centroiders[0].nmax)

            self.rtc.d_centro[0].set_npix(self.npix)

            cMat = np.zeros((nact, 2 * nvalid[0]), dtype=np.float32)

        elif p_wfs.type == WFSType.PYRHR or p_wfs.type == WFSType.PYRLR:
            nvalid = np.array([p_wfs._nvalid],
                              dtype=np.int32)  # Number of valid SUBAPERTURES
            nact = sum([p_dm.get_ntotact()
                        for p_dm in self.config.p_dms])  # Number of actu over all DMs
            gain = 1.
            self.rtc = rtc_standalone(self.c, wfsNb, nvalid, nact,
                                      self.config.p_centroiders[0].type,
                                      self.config.p_controllers[0].delay, 0, 1,
                                      brahma=self.BRAHMA, cacao=self.CACAO)

            self.rtc.d_centro[0].load_validpos(p_wfs._validsubsx, p_wfs._validsubsy,
                                               len(p_wfs._validsubsx))

            cMat = np.zeros((nact, p_wfs.nPupils * nvalid[0]), dtype=np.float32)
        else:
            raise ValueError('WFS type not supported')

        self.rtc.d_control[0].set_cmat(cMat)
        self.rtc.d_control[0].set_decayFactor(
                np.ones(nact, dtype=np.float32) * (gain - 1))
        self.rtc.d_control[0].set_matE(np.identity(nact, dtype=np.float32))
        self.rtc.d_control[0].set_mgain(np.ones(nact, dtype=np.float32) * -gain)

        self.is_init = True
