from .abstractSupervisor import AbstractSupervisor
import numpy as np

import shesha_ao as ao
import shesha_sim
import shesha_constants as scons
from shesha_constants import CONST

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
        return self.sim.config

    def enableAtmos(self, enable=True) -> None:
        ''' TODO
        Set or unset whether atmos is enabled when running loop (see singleNext) 
        '''
        self.seeAtmos = enable

    def setCommand(self, command: np.ndarray) -> None:
        ''' 
        Immediately sets provided command to DMs - does not affect integrator 
        '''
        self.sim.dms.set_full_comm((command).astype(np.float32).copy())

    def setPerturbationVoltage(self, command: np.ndarray) -> None:
        ''' 
        Add this offset value to integrator (will be applied at the end of next iteration)
        '''
        self.sim.rtc.setPertuVoltages(0, command.astype(np.float32).copy())

    def getSlope(self) -> np.ndarray:
        ''' 
        Immediately gets one slope vector for all WFS at the current state of the system 
        '''
        return self.sim.rtc.get_centroids(0)

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
        return np.empty(1)

    def singleNext(self, moveAtmos: bool=True, showAtmos: bool=True, getPSF: bool=False,
                   getResidual: bool=False) -> None:
        ''' 
        Move atmos -> getSlope -> applyControl ; One integrator step 
        '''
        self.sim.next(see_atmos=showAtmos)  # why not self.seeAtmos?

    def closeLoop(self) -> None:
        ''' 
        DM receives controller output + pertuVoltage 
        '''
        self.sim.rtc.set_openloop(0, 0)  # closeLoop

    def openLoop(self) -> None:
        ''' 
        Integrator computation goes to /dev/null but pertuVoltage still applied 
        '''
        self.sim.rtc.set_openloop(0, 1)  # openLoop

    def setRefSlopes(self, refSlopes: np.ndarray) -> None:
        ''' 
        Set given ref slopes in controller 
        '''
        self.sim.rtc.set_centroids_ref(0, refSlopes)

    def getRefSlopes(self) -> np.ndarray:
        ''' 
        Get the currently used reference slopes 
        '''
        self.sim.rtc.get_centroids_ref(0)

    def setGain(self, gain: float) -> None:
        ''' 
        Set the scalar gain of feedback controller loop 
        '''
        self.sim.rtc.set_gain(0, gain)

    def setCommandMatrix(self, cMat: np.ndarray) -> None:
        ''' 
        Set the cmat for the controller to use 
        '''
        self.sim.rtc.set_cmat(0, cMat)

    def setPyrModulation(self, pyrMod: float) -> None:
        ''' 
        Set pyramid modulation value - in l/D units 
        '''
        self.sim.rtc.set_pyr_ampl(0, pyrMod, self.sim.config.p_wfss,
                                  self.sim.config.p_tel)

    def getRawWFSImage(self, numWFS: int=0) -> np.ndarray:
        ''' 
        Get an image from the WFS 
        '''
        if self.sim.p_config.p_wfss[numWFS] == scons.WFSType.PYRHR:
            return self.sim.wfs.get_pyrimg(numWFS)
        elif self.sim.p_config.p_wfss[numWFS] == scons.WFSType.SH:
            return self.sim.wfs.get_binimg(numWFS)
        else:
            raise "WFSType not handled"

    def getTarImage(self, tarID, expoType:str="se") -> np.ndarray:
        ''' 
        Get an image from a target 
        '''
        return self.sim.tar.get_image(tarID, bytes(expoType, "utf-8"))

    def getIntensities(self) -> np.ndarray:
        ''' 
        Return sum of intensities in subaps. Size nSubaps, same order as slopes 
        '''
        return np.empty(1)

    def getAllDataLoop(self, nIter: int, slope: bool, command: bool, target: bool,
                       intensity: bool, targetPhase: bool) -> np.ndarray:
        '''
        Returns a sequence of data at continuous loop steps.
        Requires loop to be asynchronously running
        '''
        return np.empty(1)

    #  ____                  _ _   _        __  __      _   _               _     
    # / ___| _ __   ___  ___(_) |_(_) ___  |  \/  | ___| |_| |__   ___   __| |___ 
    # \___ \| '_ \ / _ \/ __| | __| |/ __| | |\/| |/ _ \ __| '_ \ / _ \ / _` / __|
    #  ___) | |_) |  __/ (__| | |_| | (__  | |  | |  __/ |_| | | | (_) | (_| \__ \
    # |____/| .__/ \___|\___|_|\__|_|\___| |_|  |_|\___|\__|_| |_|\___/ \__,_|___/
    #       |_|                                                                   

    def __init__(self, configFile:str=None, BRAHMA:bool=False):
        '''
        Init the COMPASS wih the configFile
        '''
        self.sim = None
        self.seeAtmos = False

        if configFile is not None:
            self.loadConfig(configFile, BRAHMA)

    def __repr__(self):
        print(self.sim)

    def resetDM(self, nDM: int) -> None:
        '''
        Reset the DM number nDM
        '''
        self.sim.dms.resetdm(nDM)

    def resetStrehl(self, nTar: int) -> None:
        '''
        Reset the Strehl Ratio of the target nTar
        '''
        self.sim.tar.reset_strehl(nTar)

    def loadConfig(self, configFile:str, BRAMA:bool=False) -> None:
        '''
        Init the COMPASS wih the configFile
        '''
        if self.sim is None:       
            if BRAMA:
                self.sim = shesha_sim.SimulatorBrahma(configFile)
            else:
                self.sim = shesha_sim.Simulator(configFile)
        else:
            self.sim.clear_init()
            self.sim.load_from_file(configFile)

    def isInit(self) -> bool:
        '''
        return the status on COMPASS init
        '''
        return self.sim is None

    def clearInitSim(self) -> None:
        '''
        Clear the initialization of the simulation
        '''
        self.sim.clear_init()

    def forceContext(self) -> None:
        '''
        Clear the initialization of the simulation
        '''
        self.sim.force_context()

    def initConfig(self) -> None:
        '''
        Initialize the simulation
        '''
        self.sim.init_sim()

    def getAtmScreen(self, alt:int) -> np.ndarray:
        '''
        return the atmos screen at the altitude alt
        '''
        return self.sim.atm.get_screen(alt)

    def getWfsPhase(self, numWFS:int) -> np.ndarray:        
        '''
        return the WFS screen of WFS number numWFS
        '''
        return self.sim.atm.get_screen(numWFS)

    def getDmPhase(self, dm_type:str, alt:int) -> np.ndarray:        
        '''
        return the DM screen of type dm_type conjugatide at the altitude alt
        '''
        return self.sim.dms.get_dm(dm_type, alt)

    def getTarPhase(self, numTar:int) -> np.ndarray:        
        '''
        return the target screen of target number numTar
        '''
        return self.sim.tar.get_phase(numTar)

    def getPyrHRImage(self, numWFS: int=0) -> np.ndarray:
        ''' 
        Get an HR image from the WFS 
        '''
        if self.sim.p_config.p_wfss[numWFS] == scons.WFSType.PYRHR:
            return self.sim.wfs.get_pyrimghr(numWFS)
        else:
            raise "WFSType not handled"

    def getSlopeGeom(self, numWFS:int) -> np.ndarray:        
        '''
        return the slopes geom of WFS number numWFS
        '''
        self.sim.wfs.slopes_geom(index, numTar)

        return self.sim.rtc.get_slopes(numTar)

    def getStrehl(self, numTar:int) -> np.ndarray: 
        '''
        return the Strehl Ratio of target number numTar
        '''
        self.sim.tar.comp_image(numTar)
        return self.sim.tar.get_strehl(numTar)

    def getFrameCounter(self) -> int:
        '''
        return the current frame counter of the loop
        '''
        return self.sim.iter