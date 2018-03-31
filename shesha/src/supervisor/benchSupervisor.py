from .abstractSupervisor import AbstractSupervisor
import numpy as np

from hraa.devices.camera.fakecam import Fakecam
from shesha_init import rtc_standalone, dm_init_standalone
import rtcData.DataInit as di
from shesha_sim import Simulator, SimulatorBrahma
from sutra_bind.wrap import naga_context


class BenchSupervisor(AbstractSupervisor):

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

    def enableAtmos(self, enable=True) -> None:
        ''' TODO
        Set or unset whether atmos is enabled when running loop (see singleNext)
        '''
        raise NotImplementedError("Not implemented")
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
        self.rtc.set_perturbcom(0, command.astype(np.float32).copy())

    def getSlope(self) -> np.ndarray:
        '''
        Immediately gets one slope vector for all WFS at the current state of the system
        '''
        return self.rtc.get_centroids(0)

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

    def singleNext(self, moveAtmos: bool=True, showAtmos: bool=True, getPSF: bool=False,
                   getResidual: bool=False) -> None:
        '''
        Move atmos -> getSlope -> applyControl ; One integrator step
        '''
        frame = self.cam.getFrame()
        self.rtc.load_rtc_img(0, frame.astype(np.float32))
        #for SH
        self.rtc.fill_rtc_bincube(0, self.npix)
        self.rtc.do_centroids(0)
        self.rtc.do_control(0)
        self.rtc.save_com(0)

    def closeLoop(self) -> None:
        '''
        DM receives controller output + pertuVoltage
        '''
        self.rtc.set_openloop(0, 0)  # closeLoop

    def openLoop(self) -> None:
        '''
        Integrator computation goes to /dev/null but pertuVoltage still applied
        '''
        self.rtc.set_openloop(0, 1)  # openLoop

    def setRefSlopes(self, refSlopes: np.ndarray) -> None:
        '''
        Set given ref slopes in controller
        '''
        self.rtc.set_centroids_ref(0, refSlopes)

    def getRefSlopes(self) -> np.ndarray:
        '''
        Get the currently used reference slopes
        '''
        self.rtc.get_centroids_ref(0)

    def setGain(self, gain: float) -> None:
        '''
        Set the scalar gain of feedback controller loop
        '''
        self.rtc.set_gain(gain)

    def setCommandMatrix(self, cMat: np.ndarray) -> None:
        '''
        Set the cmat for the controller to use
        '''
        self.rtc.set_cmat(cMat)

    def setPyrModulation(self, pyrMod: float) -> None:
        '''
        Set pyramid modulation value - in l/D units
        '''
        raise NotImplementedError("Not implemented")

    def getRawWFSImage(self, numWFS: int=0) -> np.ndarray:
        '''
        Get an image from the WFS
        '''
        return self.cam.getFrame()

    def getTarImage(self, tarID, expoType: str="se") -> np.ndarray:
        '''
        Get an image from a target
        '''
        raise NotImplementedError("Not implemented")

    def getIntensities(self) -> np.ndarray:
        '''
        Return sum of intensities in subaps. Size nSubaps, same order as slopes
        '''
        raise NotImplementedError("Not implemented")

    def getAllDataLoop(self, nIter: int, slope: bool, command: bool, target: bool,
                       intensity: bool, targetPhase: bool) -> np.ndarray:
        '''
        Returns a sequence of data at continuous loop steps.
        Requires loop to be asynchronously running
        '''
        raise NotImplementedError("Not implemented")

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
        self.cam = None
        self.rtc = None
        self.npix = 0

        if configFile is not None:
            self.loadConfig(configFile, BRAHMA)

    def __repr__(self):
        return str(self._sim)

    def resetDM(self, nDM: int) -> None:
        '''
        Reset the DM number nDM
        '''
        self._sim.dms.resetdm(nDM)

    def loadConfig(self, configFile: str, BRAMA: bool=False) -> None:
        '''
        Init the COMPASS wih the configFile
        '''

        if self._sim is None:
            if BRAMA:
                self._sim = SimulatorBrahma(configFile)
            else:
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
        print("->cam")
        self.cam = Fakecam(
                self._sim.config.p_cams[0].camAddr, self._sim.config.p_cams[0].width,
                self._sim.config.p_cams[0].height, self._sim.config.p_cams[0].offset_w,
                self._sim.config.p_cams[0].offset_h,
                self._sim.config.p_cams[0].expo_usec,
                self._sim.config.p_cams[0].framerate)

        dataS = di.makeSH(wfsNb=1, frameSize=self._sim.config.p_cams[0].width,
                          roiSize=self._sim.config.p_wfss[0].nxsub,
                          subSize=self._sim.config.p_wfss[0].npix)
        dataL = di.makeSCAO(dataS, cmdNb=4096, gain=0)

        self.context = naga_context(devices=np.array([0], dtype=np.int32))
        nvalid = np.array([dataS.data["roiTab"].data.shape[1]], dtype=np.int32)
        offset = dataS.param["subSize"] // 2 + 0.5
        self.rtc = rtc_standalone(self.context, dataS.param["wfsNb"], nvalid,
                                  dataL.param["cmdNb"], b"cog", 1, offset * 0,
                                  dataS.param["modFac"])
        self.npix = dataS.param["subSize"]
        cMat = dataL.data["cmdMat"].data
        s = cMat.shape
        self.rtc.set_cmat(0, cMat.reshape(s[0], s[1] * s[2] * s[3]))
        self.rtc.set_decayFactor(0,
                                 np.ones(dataL.param["cmdNb"], dtype=np.float32) *
                                 -(1 - dataL.param["gain"]))
        self.rtc.set_matE(0, np.identity(dataL.param["cmdNb"], dtype=np.float32))
        self.rtc.set_mgain(
                0,
                np.ones(dataL.param["cmdNb"], dtype=np.float32) * dataL.param["gain"] *
                (-1))
        xvalid = dataS.data["roiTab"].data[1, :] / self.npix
        yvalid = dataS.data["roiTab"].data[0, :] / self.npix
        self.rtc.load_rtc_validpos(0, xvalid.astype(np.int32), yvalid.astype(np.int32))
        self._sim.is_init = True

    def getWfsPhase(self, numWFS: int) -> np.ndarray:
        '''
        return the WFS screen of WFS number numWFS
        '''
        raise NotImplementedError("Not implemented")
        return self._sim.atm.get_screen(numWFS)

    def getDmPhase(self, dm_type: str, alt: int) -> np.ndarray:
        '''
        return the DM screen of type dm_type conjugatide at the altitude alt
        '''
        return self._sim.dms.get_dm(dm_type, alt)

    def getTarPhase(self, numTar: int) -> np.ndarray:
        '''
        return the target screen of target number numTar
        '''
        raise NotImplementedError("Not implemented")
        return self._sim.tar.get_phase(numTar)

    def getFrameCounter(self) -> int:
        '''
        return the current frame counter of the loop
        '''
        return self._sim.iter
