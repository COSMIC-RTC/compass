""" @package shesha.supervisor.compassSupervisor
Widget to simulate a closed loop

Usage:
  compassSupervisor.py [<parameters_filename>]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
"""
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

    def setOneActu(self, ndm: int, nactu: int, ampli: float = 1) -> None:
        '''
        Push the selected actuator
        '''
        self._sim.dms.d_dms[ndm].comp_oneactu(nactu, ampli)

    def setDmShapeFrom(self, command: np.ndarray) -> None:
        '''
        Immediately sets provided command to DMs - does not affect integrator
        '''
        self._sim.dms.set_full_com(command)

    def setCommand(self, nctrl: int, command: np.ndarray) -> None:
        '''
        Set the RTC command vector
        '''
        self._sim.rtc.d_control[nctrl].set_com(command, command.size)

    def setPerturbationVoltage(self, nControl: int, name: str,
                               command: np.ndarray) -> None:
        '''
        Add this offset value to integrator (will be applied at the end of next iteration)
        '''
        if len(command.shape) == 1:
            self._sim.rtc.d_control[nControl].set_perturb_voltage(name, command, 1)
        elif len(command.shape) == 2:
            self._sim.rtc.d_control[nControl].set_perturb_voltage(
                    name, command, command.shape[0])
        else:
            raise AttributeError("command should be a 1D or 2D array")

    def resetPerturbationVoltage(self, nControl: int) -> None:
        '''
        Reset the perturbation voltage of the nControl controller
        (i.e. will remove ALL perturbation voltages.)
	If you want to reset just one, see the function removePerturbationVoltage().
        '''
        self._sim.rtc.d_control[nControl].reset_perturb_voltage()

    def removePerturbationVoltage(self, nControl: int, name: str) -> None:
        '''
        Remove the perturbation voltage called <name>, from the
 	controller number <nControl>.
	If you want to remove all of them, see function resetPerturbationVoltage().

        '''
        self._sim.rtc.d_control[nControl].remove_perturb_voltage(name)

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

    def singleNext(self, moveAtmos: bool = True, showAtmos: bool = True,
                   getPSF: bool = False, getResidual: bool = False) -> None:
        '''
        Move atmos -> getSlope -> applyControl ; One integrator step
        '''
        self._sim.next(see_atmos=showAtmos)  # why not self._seeAtmos?

    def closeLoop(self, ncontrol: int = 0) -> None:
        """
        DM receives controller output + pertuVoltage
        Parameters:
            ncontrol: (int, optionnal): controller index. Default is 0
        """
        self._sim.rtc.d_control[ncontrol].set_openloop(0)  # closeLoop

    def openLoop(self, rst=True, ncontrol: int = 0) -> None:
        '''
        Integrator computation goes to /dev/null but pertuVoltage still applied
        '''
        self._sim.rtc.d_control[ncontrol].set_openloop(1, rst)  # openLoop

    def setRefSlopes(self, refSlopes: np.ndarray) -> None:
        '''
        Set given ref slopes in controller
        '''
        self._sim.rtc.set_centroids_ref(refSlopes)

    def getRefSlopes(self) -> np.ndarray:
        '''
        Get the currently used reference slopes
        '''
        refSlopes = np.empty(0)
        for centro in self._sim.rtc.d_centro:
            refSlopes = np.append(refSlopes, np.array(centro.d_centroids_ref))
        return refSlopes

    def setGain(self, gainMat, ncontrol = 0) -> None:
        '''
        Set the scalar gain of feedback controller loop
        '''
        if type(gainMat) in [int, float]:
            gainMat = np.ones(
                    np.sum(self._sim.config.p_controller0.nactu),
                    dtype=np.float32) * gainMat
        self._sim.rtc.d_control[ncontrol].set_mgain(gainMat)

    def setCommandMatrix(self, cMat: np.ndarray , nControl: int = 0) -> None:
        '''
        Set the cmat for the controller to use
        '''
        self._sim.rtc.d_control[nControl].set_cmat(cMat)

    def setNoise(self, noise, numwfs=0, seed=1234):
        '''
        Set noise value of WFS numwfs
        '''
        self._sim.wfs.d_wfs[numwfs].set_noise(noise, int(seed + numwfs))
        print("Noise set to: %f on WFS %d" % (noise, numwfs))

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
        print("PYR method set to " + self._sim.rtc.d_centro[0].pyr_method)

    def getPyrMethod(self):
        return self._sim.rtc.d_centro[0].pyr_method

    def setGSmag(self, mag, numwfs=0):
        numwfs = int(numwfs)
        sim = self._sim
        wfs = sim.wfs.d_wfs[numwfs]
        if (sim.config.p_wfs0.type == "pyrhr"):
            r = wfs.comp_nphot(sim.config.p_loop.ittime,
                               sim.config.p_wfss[numwfs].optthroughput,
                               sim.config.p_tel.diam, sim.config.p_tel.cobs,
                               sim.config.p_wfss[numwfs].zerop, mag)
        else:
            r = wfs.comp_nphot(sim.config.p_loop.ittime,
                               sim.config.p_wfss[numwfs].optthroughput,
                               sim.config.p_tel.diam, sim.config.p_wfss[numwfs].nxsub,
                               sim.config.p_wfss[numwfs].zerop, mag)
        if (r == 0):
            print("GS magnitude is now %f on WFS %d" % (mag, numwfs))

    def getTarImage(self, tarID, expoType: str = "se") -> np.ndarray:
        '''
        Get an image from a target
        '''
        if (expoType == "se"):
            return np.fft.fftshift(np.array(self._sim.tar.d_targets[tarID].d_image_se))
        elif (expoType == "le"):
            return np.fft.fftshift(np.array(self._sim.tar.d_targets[tarID].d_image_le)
                                   ) / self._sim.tar.d_targets[tarID].strehl_counter
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

    def __init__(self, configFile: str = None, cacao: bool = False,
                 use_DB: bool = False):
        '''
        Init the COMPASS supervisor

        Parameters
        ------------
        configFile: (str): (optionnal) Path to the parameter file
        cacao: (bool): (optionnal) Flag to enable cacao
        use_DB: (bool): (optionnal) Flag to enable database
        '''
        self._sim = None
        self._seeAtmos = False
        self.config = None
        self.cacao = cacao
        self.use_DB = use_DB

        if configFile is not None:
            self.loadConfig(configFile=configFile)

    def __repr__(self):
        return object.__repr__(self) + str(self._sim)

    def loop(self, n: int = 1, monitoring_freq: int = 100, **kwargs):
        """
        Perform the AO loop for n iterations

        :parameters:
            n: (int): (optional) Number of iteration that will be done
            monitoring_freq: (int): (optional) Monitoring frequency [frames]
        """
        self._sim.loop(n, monitoring_freq=monitoring_freq, **kwargs)

    def forceContext(self) -> None:
        '''
        Clear the initialization of the simulation
        '''
        self._sim.force_context()

    def computeSlopes(self, ncontrol = 0):
        for w in self._sim.wfs.d_wfs:
            w.d_gs.comp_image()
        self._sim.rtc.do_centroids(ncontrol)
        return np.array(self._sim.rtc.d_control[ncontrol].d_centroids)

    def resetDM(self, numdm: int = -1) -> None:
        '''
        Reset the DM number nDM or all DMs if  == -1
        '''
        if (numdm == -1):  # All Dms reset
            for dm in self._sim.dms.d_dms:
                dm.reset_shape()
        else:
            self._sim.dms.d_dms[numdm].reset_shape()

    def resetCommand(self, nctrl: int = -1) -> None:
        '''
        Reset the nctrl Controller command buffer, reset all controllers if nctrl  == -1
        '''
        if (nctrl == -1):  # All Dms reset
            for control in self._sim.rtc.d_control:
                control.d_com.reset()
        else:
            self._sim.rtc.d_control[nctrl].d_com.reset()

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

    def resetTarPhase(self, nTar: int) -> None:
        '''
        Reset the phase screen of the target nTar
        '''
        self._sim.tar.d_targets[nTar].d_phase.reset()

    def loadConfig(self, configFile: str = None, sim=None) -> None:
        '''
        Init the COMPASS simulator wih the configFile
        '''
        if self._sim is None:
            if sim is None:
                if self.cacao:
                    from shesha.sim.simulatorCacao import SimulatorCacao as Simulator
                else:
                    from shesha.sim.simulator import Simulator
                self._sim = Simulator(filepath=configFile, use_DB=self.use_DB)
            else:
                self._sim = sim
        else:
            self._sim.clear_init()
            self._sim.load_from_file(configFile)
        self.config = self._sim.config

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

    def initConfig(self) -> None:
        '''
        Initialize the simulation
        '''
        self._sim.init_sim()
        self.enableAtmos(True)

    def getNcpaWfs(self, wfsnum):
        return np.array(self._sim.wfs.d_wfs[wfsnum].d_gs.d_ncpa_phase)

    def getNcpaTar(self, tarnum):
        return np.array(self._sim.tar.d_targets[tarnum].d_ncpa_phase)

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

    def getDmShape(self, indx: int) -> np.ndarray:
        '''
        return the selected DM screen
        '''
        return np.array(self._sim.dms.d_dms[indx].d_shape)

    def getTarPhase(self, numTar: int) -> np.ndarray:
        '''
        return the target screen of target number numTar
        '''
        return np.array(self._sim.tar.d_targets[numTar].d_phase)

    def getPyrHRImage(self, numWFS: int = 0) -> np.ndarray:
        '''
        Get an HR image from the WFS
        '''
        return np.array(self._sim.wfs.d_wfs[numWFS].d_hrimg)

    def getWfsImage(self, numWFS: int = 0) -> np.ndarray:
        '''
        Get an image from the WFS
        '''
        return np.array(self._sim.wfs.d_wfs[numWFS].d_binimg)

    def getSlopeGeom(self, numWFS: int, ncontrol: int = 0) -> np.ndarray:
        '''
        return the slopes geom of WFS number numWFS
        '''
        self._sim.rtc.do_centroids_geom(ncontrol)
        slopesGeom = np.array(self._sim.rtc.d_control[ncontrol].d_centroids)
        self._sim.rtc.do_centroids(ncontrol)
        return slopesGeom

    def getStrehl(self, numTar: int, do_fit: bool = True) -> np.ndarray:
        '''
        return the Strehl Ratio of target number numTar
        '''
        src = self._sim.tar.d_targets[numTar]
        src.comp_strehl(do_fit)
        avgVar = 0
        if (src.phase_var_count > 0):
            avgVar = src.phase_var_avg / src.phase_var_count
        return [src.strehl_se, src.strehl_le, src.phase_var, avgVar]

    def getFrameCounter(self) -> int:
        '''
        return the current frame counter of the loop
        '''
        return self._sim.iter

    def getIFsparse(self, nControl: int):
        '''
        Return the IF of DM as a sparse matrix
        '''
        return self._sim.rtc.d_control[nControl].d_IFsparse.get_csr()

    def getIFtt(self, nControl: int):
        '''
        Return the IF of a TT DM as a sparse matrix
        '''
        return np.array(self._sim.rtc.d_control[nControl].d_TT)

    def getIFdm(self, nDM: int):
        '''
        Return the IF of a Petal DM made with M4
        '''

        from shesha.ao import basis
        if_sparse = basis.compute_DMbasis(self._sim.dms.d_dms[nDM], self._sim.config.p_dms[nDM], self._sim.config.p_geom)

        return if_sparse

    def getCentroids(self, nControl: int):
        '''
        Return the centroids of the nControl controller
        '''
        return np.array(self._sim.rtc.d_control[nControl].d_centroids)

    def getCom(self, nControl: int):
        '''
        Get command from nControl controller
        '''
        return np.array(self._sim.rtc.d_control[nControl].d_com)

    def getErr(self, nControl: int):
        '''
        Get command increment from nControl controller
        '''
        return np.array(self._sim.rtc.d_control[nControl].d_err)

    def getVoltage(self, nControl: int):
        '''
        Get voltages from nControl controller
        '''
        return np.array(self._sim.rtc.d_control[nControl].d_voltage)

    def setIntegratorLaw(self, ncontrol = 0):
        self._sim.rtc.d_control[ncontrol].set_commandlaw("integrator")

    def setDecayFactor(self, decay, ncontrol = 0):
        self._sim.rtc.d_control[ncontrol].set_decayFactor(decay)

    def setEMatrix(self, eMat, ncontrol = 0):
        self._sim.rtc.d_control[ncontrol].set_matE(eMat)

    def doRefslopes(self):
        print("Doing refslopes...")
        self._sim.rtc.do_centroids_ref(0)
        print("refslopes done")

    def resetRefslopes(self):
        for centro in self._sim.rtc.d_centro:
            centro.d_centroids_ref.reset()

    def setNcpaWfs(self, ncpa, wfsnum):
        self._sim.wfs.d_wfs[wfsnum].d_gs.set_ncpa(ncpa)

    def setNcpaTar(self, ncpa, tarnum):
        self._sim.tar.d_targets[tarnum].set_ncpa(ncpa)

    def setWfsPhase(self, numwfs, phase):
        self._sim.wfs.d_wfs[numwfs].d_gs.set_phase(phase)

    def setMpupil(self, mpupil, numwfs=0):
        oldmpup = self.getMpupil()
        dimx = oldmpup.shape[0]
        dimy = oldmpup.shape[1]
        if ((mpupil.shape[0] != dimx) or (mpupil.shape[1] != dimy)):
            print("Error mpupil shape on wfs %d must be: (%d,%d)" % (numwfs, dimx, dimy))
        else:
            self._sim.wfs.d_wfs[numwfs].set_pupil(mpupil.copy())

    def getIpupil(self):
        return self._sim.config.p_geom._ipupil

    def getSpupil(self):
        return self._sim.config.p_geom._spupil

    def getMpupil(self):
        return self._sim.config.p_geom._mpupil

    def getTarAmplipup(self, tarnum):
        return self._sim.config.tar.get_amplipup(tarnum)

    def getImat(self, nControl: int):
        """
        Return the interaction matrix of the controller

        Parameters
        ------------
        nControl: (int): controller index
        """
        return np.array(self._sim.rtc.d_control[nControl].d_imat)

    def getCmat(self, nControl: int):
        """
        Return the command matrix of the controller

        Parameters
        ------------
        nControl: (int): controller index
        """
        return np.array(self._sim.rtc.d_control[nControl].d_cmat)

    def setCentroThresh(self, nCentro: int, thresh: float):
        """
        Set the threshold value of a thresholded COG

        Parameters
        ------------
        nCentro: (int): centroider index
        thresh: (float): new threshold value
        """
        self._sim.rtc.d_centro[nCentro].set_threshold(thresh)

    def getPyrFocalPlane(self, nwfs: int = 0):
        """
        No arguments
        Returns the psf in the focal plane of the pyramid.
        """
        return np.fft.fftshift(np.array(self._sim.wfs.d_wfs[nwfs].d_pyrfocalplane))

    def reset(self,tar=-1,rst=True):
        self.resetTurbu()
        if(tar<0):
            for tar in range(self._sim.tar.ntargets):
                self.resetStrehl(tar)
        else:
            self.resetStrehl(tar)
        self.resetDM()
        self.openLoop(rst=rst)
        self.closeLoop()
