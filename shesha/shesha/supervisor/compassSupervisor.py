from .abstractSupervisor import AbstractSupervisor
import numpy as np

import shesha.constants as scons
from shesha.constants import CONST

from tqdm import trange


class CompassSupervisor(AbstractSupervisor):

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

    def getConfig(self):
        '''
        Returns the configuration in use, in a supervisor specific format ?
        '''
        return self._sim.config

    def enableAtmos(self, enable) -> None:
        ''' TODO
        Set or unset whether atmos is enabled when running loop (see singleNext)
        '''
        self._seeAtmos = enable

    def setCommand(self, command: np.ndarray) -> None:
        '''
        Immediately sets provided command to DMs - does not affect integrator
        '''
        self._sim.dms.set_full_comm((command).astype(np.float32).copy())

    def setPerturbationVoltage(self, command: np.ndarray) -> None:
        '''
        Add this offset value to integrator (will be applied at the end of next iteration)
        '''
        self._sim.rtc.d_control[0].set_perturbcom(command.astype(np.float32).copy())

    def getSlope(self) -> np.ndarray:
        '''
        Immediately gets one slope vector for all WFS at the current state of the system
        '''
        return self.computeSlopes()

    def computeIMatModal(self, M2V: np.ndarray, pushVector: np.ndarray,
                         refOffset: np.ndarray, noise: bool,
                         useAtmos: bool) -> np.ndarray:
        '''
        TODO
        Computes a modal interaction matrix for the given modal matrix
        with given push values (length = nModes)
        around an (optional) offset value
        optionally with noise
        with/without atmos shown to WFS
        '''
        raise NotImplementedError("Not implemented")
        # return np.empty(1)

    def next(self, nbiters, see_atmos=True):
        for _ in trange(nbiters):
            self._sim.next(see_atmos=see_atmos)

    def singleNext(self, moveAtmos: bool=True, showAtmos: bool=True, getPSF: bool=False,
                   getResidual: bool=False) -> None:
        '''
        Move atmos -> getSlope -> applyControl ; One integrator step
        '''
        self._sim.next(see_atmos=showAtmos)  # why not self._seeAtmos?

    def closeLoop(self) -> None:
        '''
        DM receives controller output + pertuVoltage
        '''
        self._sim.rtc.d_control[0].set_openloop(0)  # closeLoop

    def openLoop(self) -> None:
        '''
        Integrator computation goes to /dev/null but pertuVoltage still applied
        '''
        self._sim.rtc.d_control[0].set_openloop(1)  # openLoop

    def setRefSlopes(self, refSlopes: np.ndarray) -> None:
        '''
        Set given ref slopes in controller
        '''
        self._sim.rtc.d_control[0].set_centroids_ref(refSlopes)

    def getRefSlopes(self) -> np.ndarray:
        '''
        Get the currently used reference slopes
        '''
        return np.array(self._sim.rtc.d_control[0].d_centroids_ref)

    def setGain(self, gainMat) -> None:
        '''
        Set the scalar gain of feedback controller loop
        '''
        if type(gainMat) in [int, float]:
            gainMat = np.ones(
                    np.sum(self._sim.config.p_controller0.nactu),
                    dtype=np.float32) * gainMat
        self._sim.rtc.d_control[0].set_mgain(gainMat)

    def setCommandMatrix(self, cMat: np.ndarray) -> None:
        '''
        Set the cmat for the controller to use
        '''
        self._sim.rtc.d_control[0].set_cmat(cMat)

    def setNoise(self, noise, numwfs=0):
        '''
        Set noise value of WFS numwfs
        '''
        self._sim.wfs.d_wfs[numwfs].set_noise(noise)
        print("Noise set to: %d" % noise)

    def setPyrModulation(self, pyrMod: float) -> None:
        '''
        Set pyramid modulation value - in l/D units
        '''
        from shesha.ao.wfs import comp_new_pyr_ampl

        _, _, _, pyr_npts = comp_new_pyr_ampl(0, pyrMod, self._sim.wfs, self._sim.rtc,
                                              self._sim.config.p_wfss,
                                              self._sim.config.p_tel)

        print("PYR modulation set to: %f L/D using %d points" % (pyrMod, pyr_npts))

    def setPyrMethod(self, pyrMethod):
        '''
        Set pyramid compute method
        '''
        self._sim.rtc.d_centro[0].set_pyr_method(pyrMethod)  # Sets the pyr method
        print("PYR method set to: %d" % self._sim.rtc.d_centro[0].pyr_method)

    def getRawWFSImage(self, numWFS: int=0) -> np.ndarray:
        '''
        Get an image from the WFS
        '''
        return np.array(self._sim.wfs.d_wfs[numWFS].d_binimg)

    def getTarImage(self, tarID, expoType: str="se") -> np.ndarray:
        '''
        Get an image from a target
        '''
        self._sim.tar.d_targets[0].comp_image()
        if (expoType == "se"):
            return np.fft.fftshift(np.array(self._sim.tar.d_targets[0].d_image_se))
        elif (expoType == "le"):
            return np.fft.fftshift(np.array(self._sim.tar.d_targets[0].d_image_le))
        else:
            raise ValueError("Unknown exposure type")

    def getIntensities(self) -> np.ndarray:
        '''
        Return sum of intensities in subaps. Size nSubaps, same order as slopes
        '''
        raise NotImplementedError("Not implemented")
        # return np.empty(1)

    def getAllDataLoop(self, nIter: int, slope: bool, command: bool, target: bool,
                       intensity: bool, targetPhase: bool) -> np.ndarray:
        '''
        Returns a sequence of data at continuous loop steps.
        Requires loop to be asynchronously running
        '''
        raise NotImplementedError("Not implemented")
        # return np.empty(1)

    #  ____                  _ _   _        __  __      _   _               _
    # / ___| _ __   ___  ___(_) |_(_) ___  |  \/  | ___| |_| |__   ___   __| |___
    # \___ \| '_ \ / _ \/ __| | __| |/ __| | |\/| |/ _ \ __| '_ \ / _ \ / _` / __|
    #  ___) | |_) |  __/ (__| | |_| | (__  | |  | |  __/ |_| | | | (_) | (_| \__ \
    # |____/| .__/ \___|\___|_|\__|_|\___| |_|  |_|\___|\__|_| |_|\___/ \__,_|___/
    #       |_|

    def __init__(self, configFile: str=None, BRAHMA: bool=False):
        '''
        Init the COMPASS wih the configFile
        '''
        self._sim = None
        self._seeAtmos = False

        if configFile is not None:
            self.loadConfig(configFile, BRAHMA)

    def __repr__(self):
        return str(self._sim)

    def computeSlopes(self):
        for w in self._sim.wfs.d_wfs:
            w.d_gs.comp_image()
        self._sim.rtc.do_centroids(0)
        return np.array(self._sim.rtc.d_control[0].d_centroids)

    def resetDM(self, numdm: int=-1) -> None:
        '''
        Reset the DM number nDM or all DMs if  == -1
        '''
        if (numdm == -1):  # All Dms reset
            for dm in self._sim.dms.d_dms:
                dm.reset_shape()
        else:
            self._sim.dms[numdm].reset_shape()

    def resetSimu(self, noiseList):
        self.resetTurbu()
        self.resetNoise(noiseList)

    def resetTurbu(self):
        ilayer = 0
        for k in range(self._sim.atm.nscreens):
            self._sim.atm.set_seed(k, 1234 + ilayer)
            self._sim.atm.refresh_screen(k)
            ilayer += 1

    def resetNoise(self, noiseList):
        for nwfs in range(len(self._sim.config.p_wfss)):
            self._sim.wfs.d_wfs[nwfs].set_noise(noiseList[nwfs], 1234 + nwfs)

    def resetStrehl(self, nTar: int) -> None:
        '''
        Reset the Strehl Ratio of the target nTar
        '''
        self._sim.tar.d_targets[nTar].reset_strehlmeter()

    def loadConfig(self, configFile: str, BRAMA: bool=False) -> None:
        '''
        Init the COMPASS wih the configFile
        '''
        if self._sim is None:
            if BRAMA:
                from shesha.sim.simulatorBrahma import SimulatorBrahma
                self._sim = SimulatorBrahma(configFile)
            else:
                from shesha.sim.simulator import Simulator
                self._sim = Simulator(configFile)
        else:
            self._sim.clear_init()
            self._sim.load_from_file(configFile)

    def isInit(self) -> bool:
        '''
        return the status on COMPASS init
        '''
        return self._sim.is_init

    def clearInitSim(self) -> None:
        '''
        Clear the initialization of the simulation
        '''
        self._sim.clear_init()

    def forceContext(self) -> None:
        '''
        Clear the initialization of the simulation
        '''
        self._sim.force_context()

    def initConfig(self) -> None:
        '''
        Initialize the simulation
        '''
        self._sim.init_sim()

    def getAtmScreen(self, indx: int) -> np.ndarray:
        '''
        return the selected atmos screen
        '''
        return np.array(self._sim.atm.d_screens[indx].d_screen)

    def getWfsPhase(self, numWFS: int) -> np.ndarray:
        '''
        return the WFS screen of WFS number numWFS
        '''
        return np.array(self._sim.wfs.d_wfs[numWFS].d_gs.d_phase)

    def getDmPhase(self, indx: int) -> np.ndarray:
        '''
        return the selected DM screen
        '''
        return np.array(self._sim.dms.d_dms[indx].d_shape)

    def getTarPhase(self, numTar: int) -> np.ndarray:
        '''
        return the target screen of target number numTar
        '''
        return np.array(self._sim.tar.d_targets[numTar].d_phase)

    def getPyrHRImage(self, numWFS: int=0) -> np.ndarray:
        '''
        Get an HR image from the WFS
        '''
        return np.array(self._sim.wfs.d_wfs[numWFS].d_hrimg)

    def getSlopeGeom(self, numWFS: int) -> np.ndarray:
        '''
        return the slopes geom of WFS number numWFS
        '''
        self._sim.rtc.do_centroids_geom(0)
        slopesGeom = np.array(self._sim.rtc.d_control[0].d_centroids)
        self._sim.rtc.do_centroids(0)
        return slopesGeom

    def getStrehl(self, numTar: int) -> np.ndarray:
        '''
        return the Strehl Ratio of target number numTar
        '''
        src = self._sim.tar.d_targets[numTar]
        src.comp_image()
        src.comp_strehl()
        avgVar = 0
        if (src.phase_var_count > 0):
            avgVar = src.phase_var_avg / src.phase_var_count
        return [src.strehl_se, src.strehl_le, src.phase_var, avgVar]

    def getFrameCounter(self) -> int:
        '''
        return the current frame counter of the loop
        '''
        return self._sim.iter
