"""Widget to simulate a closed loop

Usage:
  widget_fab.py [<parameters_filename>] [--expert]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
"""

import os, sys
import numpy as np
import time
from collections import OrderedDict

import pyqtgraph as pg
sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/AOlib")
sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/src/shesha_util")
from tools import plsh, plpyr
from tqdm import tqdm
import threading
import astropy.io.fits as pfits
from PyQt5 import QtGui, QtWidgets
from PyQt5.uic import loadUiType
from PyQt5.QtCore import QThread, QObject, QTimer, pyqtSignal
from threading import Thread
from supervisor.abstractSupervisor import AbstractSupervisor
from functools import partial
from subprocess import Popen, PIPE

import shesha_ao as ao
import shesha_sim
import shesha_constants as scons
from shesha_constants import CONST

import compassConfigToFile as cf

from typing import Any, Dict, Tuple, Callable, List
"""
low levels debugs:
gdb --args python -i widget_ao.py
"""
from docopt import docopt

import Pyro4
"""
IMPORTANT PYRO V4: To re-enable the np.array serialisation
add in the .bashrc (before launching Pyro4 NS):
export PYRO_SERIALIZERS_ACCEPTED=serpent,json,marshal,pickle,dill
export PYRO_LOGFILE=pyro.log
export PYRO_LOGLEVEL=DEBUG

Add on the python Server Side: Pyro4.config.SERIALIZERS_ACCEPTED = set(['pickle','json', 'marshal', 'serpent'])
Add on the python client Side: Pyro4.config.SERIALIZER='pickle'
"""
Pyro4.config.SERIALIZERS_ACCEPTED = set(['pickle', 'json', 'marshal', 'serpent'])

sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/data/par/")
WindowTemplate, TemplateBaseClass = loadUiType(
        os.environ["SHESHA_ROOT"] + "/widgets/widget_ao.ui")  # type: type, type
from widget_ao import widgetAOWindow


class MergeMetaClass(type(AbstractSupervisor), type(widgetAOWindow)):
    pass


@Pyro4.expose
class widgetAOWindowPyro(AbstractSupervisor, widgetAOWindow, metaclass=MergeMetaClass):

    def __init__(self, configFile: Any=None, BRAMA: bool=False,
                 expert: bool=False) -> None:
        widgetAOWindow.__init__(self, configFile, BRAMA)
        #Pyro.core.ObjBase.__init__(self)

        self.CB = {}
        self.wpyr = None
        #############################################################
        #                 CONNECTED BUTTONS                         #
        #############################################################
        # Default path for config files

        self.ui.actionShow_Pyramid_Tools.toggled.connect(self.showPyrTools)
        self.ph2modes = None
        self.KL2V = None
        self.P = None
        self.currentBuffer = 1
        self.wpyrNbBuffer = 1
        #############################################################
        #                       METHODS                             #
        #############################################################

    """
         ____  _   _ ____  _____ ______     _____ ____   ___  ____
        / ___|| | | |  _ \| ____|  _ \ \   / /_ _/ ___| / _ \|  _ \
        \___ \| | | | |_) |  _| | |_) \ \ / / | |\___ \| | | | |_) |
         ___) | |_| |  __/| |___|  _ < \ V /  | | ___) | |_| |  _ <
        |____/ \___/|_|   |_____|_| \_\ \_/  |___|____/ \___/|_| \_\

            _    ____ ____ _____ ____      _    ____ _____
           / \  | __ ) ___|_   _|  _ \    / \  / ___|_   _|
          / _ \ |  _ \___ \ | | | |_) |  / _ \| |     | |
         / ___ \| |_) |__) || | |  _ <  / ___ \ |___  | |
        /_/   \_\____/____/ |_| |_| \_\/_/   \_\____| |_|

         __  __ _____ _____ _   _  ___  ____  ____
        |  \/  | ____|_   _| | | |/ _ \|  _ \/ ___|
        | |\/| |  _|   | | | |_| | | | | | | \___ \
        | |  | | |___  | | |  _  | |_| | |_| |___) |
        |_|  |_|_____| |_| |_| |_|\___/|____/|____/

    """

    def getConfig(self, path):
        ''' Returns the configuration in use, in a supervisor specific format '''
        return self.writeConfigOnFile(path)

    def loadConfig(self) -> None:
        ''' Load the configuration for the supervisor'''
        widgetAOWindow.loadConfig(self)
        print("switching to a generic controller")
        self.sim.config.p_controllers[0].type = scons.ControllerType.GENERIC

    def initConfig(self):
        ''' Init the configuration for the supervisor'''
        self.InitConfig()

    def setCommand(self, commandsVector):
        ''' Immediately sets provided command to DMs - does not affect integrator '''
        self.sim.dms.set_full_comm(commandsVector)

    def setPerturbationVoltage(self, actusVolts):
        ''' Add this offset value to integrator (will be applied at the end of next iteration)'''
        self.sim.rtc.set_perturbcom(0, actusVolts.astype(np.float32).copy())

    def getSlope(self):
        ''' Immediately gets one slope vector for all WFS at the current state of the system '''
        return self.computeSlopes()

    def singleNext(self, moveAtmos: bool=True, showAtmos: bool=True, getPSF: bool=False,
                   getResidual: bool=False):
        ''' Move atmos -> getSlope -> applyControl ; One integrator step '''
        self.loopOnce()

    def closeLoop(self):
        ''' DM receives controller output + (optional) pertuVoltage '''
        self.sim.rtc.set_openloop(0, 0)

    def openLoop(self):
        ''' Integrator computation goes to /dev/null but pertuVoltage still applied (if any)'''
        self.sim.rtc.set_openloop(0, 1)

    def setRefSlopes(self, refSlopes):
        ''' Set given ref slopes in controller '''
        self.loadRefSlopes(refSlopes)

    def getRefSlopes(self):
        ''' Get the currently used reference slopes '''
        raise ValueError("ERROR getRefSlopes not Implemented")

    def setGain(self, gainMat):
        ''' Set the gain of feedback controller loop '''
        if ((type(gainMat) is float) or (type(gainMat) is int)):
            gainMat = np.ones(
                    np.sum(self.sim.config.p_controller0.nactu),
                    dtype=np.float32) * gainMat
        self.sim.rtc.set_mgain(0, gainMat)

    def setCommandMatrix(self, cMat):
        ''' Set the cmat for the controller to use '''
        return self.sim.rtc.set_cmat(0, cMat)

    def getTarImage(self, tarID):
        ''' Get an image from a target '''
        return self.getTarImage(tarID)

    def getIntensities(self):
        ''' Return sum of intensities in subaps. Size nSubaps, same order as slopes '''
        raise ValueError("ERROR getIntensities not Implemented")

    def getAllDataLoop(self, nIter: int, slope: bool, command: bool, target: bool,
                       intensity: bool, targetPhase: bool) -> np.ndarray:
        '''
        Returns a sequence of data at continuous loop steps.
        Requires loop to be asynchronously running
        '''
        raise ValueError("ERROR getAllDataLoop not Implemented")

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

    def InitConfig(self):
        try:
            ps = PyroServer(self)
            ps.start()
        except:
            print("Warning: Error while starting Pyro server")
        widgetAOWindow.InitConfig(self)

    def InitConfigFinished(self) -> None:
        widgetAOWindow.InitConfigFinished(self)

    def loopOnce(self) -> None:
        widgetAOWindow.loopOnce(self)
        start = time.time()
        refreshDisplayTime = 1. / self.ui.wao_frameRate.value()

        if (self.ui.actionShow_Pyramid_Tools.isChecked()):  # PYR only
            self.wpyr.Fe = 1 / self.sim.config.p_loop.ittime  # needs Fe for PSD...
            if (self.wpyr.CBNumber == 1):
                self.ai = self.computeModalResiduals()
                self.setPyrToolsParams(self.ai)
            else:
                if (self.currentBuffer == 1):  # First iter of the CB
                    aiVect = self.computeModalResiduals()
                    self.ai = aiVect[np.newaxis, :]
                    self.currentBuffer += 1  # Keep going

                else:  # Keep filling the CB
                    aiVect = self.computeModalResiduals()
                    self.ai = np.concatenate((self.ai, aiVect[np.newaxis, :]))
                    if (self.currentBuffer < self.wpyr.CBNumber):
                        self.currentBuffer += 1  # Keep going
                    else:
                        self.currentBuffer = 1  # reset buffer
                        self.setPyrToolsParams(self.ai)  # display

    def next(self, nbIters, see_atmos=True):
        ''' Move atmos -> getSlope -> applyControl ; One integrator step '''
        for i in tqdm(range(nbIters)):
            self.sim.next(see_atmos=see_atmos)

    def resetDm(self, numdm=-1):
        if (numdm == -1):  # All Dms reset
            for numdm in range(len(self.sim.config.p_dms)):
                self.sim.dms.resetdm(self.sim.config.p_dms[numdm].type,
                                     self.sim.config.p_dms[numdm].alt)
        else:
            self.sim.dms.resetdm(self.sim.config.p_dms[numdm].type,
                                 self.sim.config.p_dms[numdm].alt)

    def resetSimu(self, noiseList):
        self.resetTurbu()
        time.sleep(1)
        self.resetNoise(noiseList)

    def resetTurbu(self):
        ilayer = 0
        for layerAlt in self.sim.atm.list_alt().tolist():
            self.sim.atm.set_seed(layerAlt, 1234 + ilayer)
            self.sim.atm.refresh_screen(layerAlt)
            ilayer += 1

    def resetNoise(self, noiseList):
        for nwfs in range(len(self.sim.config.p_wfss)):
            self.sim.wfs.set_noise(nwfs, noiseList[nwfs], 1234 + nwfs)

    def initPyrTools(self):
        ADOPTPATH = os.getenv("ADOPTPATH")
        sys.path.append(ADOPTPATH + "/widgets")
        from pyrStats import widget_pyrStats
        print("OK Pyramid Tools Widget initialized")
        self.wpyr = widget_pyrStats()
        self.wpyrNbBuffer = self.wpyr.CBNumber
        self.wpyr.show()

    def setPyrToolsParams(self, ai):
        self.wpyr.pup = self.getSpupil()
        self.wpyr.phase = self.getTargetPhase(0)
        self.wpyr.updateResiduals(ai)
        if (self.ph2modes is None):
            print('computing phase 2 Modes basis')
            self.computePh2Modes()
        self.wpyr.ph2modes = self.ph2modes

    def showPyrTools(self):
        if (self.wpyr is None):
            try:
                print("Lauching pyramid widget...")
                self.initPyrTools()
                print("Done")
            except:
                raise ValueError("ERROR: ADOPT  not found. Cannot launch Pyramid tools")
        else:
            if (self.ui.actionShow_Pyramid_Tools.isChecked()):
                self.wpyr.show()
            else:
                self.wpyr.hide()

    def computePh2Modes(self):
        self.aoLoopClicked(False)
        self.ui.wao_run.setChecked(False)
        oldnoise = self.sim.config.p_wfs0.noise
        self.setNoise(-1)

        if (self.modalBasis is None):
            self.modalBasis, self.P = self.getModes2VBasis("Btt")
        nbmode = self.modalBasis.shape[1]
        pup = self.sim.config.p_geom._spupil
        ph = self.sim.tar.get_phase(0)
        ph2KL = np.zeros((nbmode, ph.shape[0], ph.shape[1]))
        S = np.sum(pup)
        for i in range(nbmode):
            self.sim.tar.reset_phase(0)
            self.sim.dms.set_full_comm((self.modalBasis[:, i]).astype(np.float32).copy())
            self.sim.next(see_atmos=False)
            ph = self.sim.tar.get_phase(0) * pup
            # Normalisation pour les unites rms en microns !!!
            norm = np.sqrt(np.sum((ph)**2) / S)
            ph2KL[i] = ph / norm
            self.printInPlace("mode #%d/%d" % (i, nbmode))
        self.ph2modes = ph2KL
        self.setNoise(oldnoise)
        self.aoLoopClicked(True)
        self.ui.wao_run.setChecked(True)
        return ph2KL

    def getModes2VBasis(self, ModalBasisType):
        if (ModalBasisType == "KL2V"):
            print("Computing KL2V basis...")
            self.modalBasis, _ = self.returnkl2V()
            return self.modalBasis, 0
        elif (ModalBasisType == "Btt"):
            print("Computing Btt basis...")
            self.modalBasis, self.P = self.compute_Btt2()
            return self.modalBasis, self.P

    def returnkl2V(self):
        """
        KL2V = ao.compute_KL2V(self.sim.config.p_controllers[
                               0], self.sim.dms, self.sim.config.p_dms, self.sim.config.p_geom, self.sim.config.p_atmos, self.sim.config.p_tel)
        """
        if (self.KL2V is None):
            print("Computing KL2V...")
            KL2V = ao.compute_KL2V(self.sim.config.p_controllers[0], self.sim.dms,
                                   self.sim.config.p_dms, self.sim.config.p_geom,
                                   self.sim.config.p_atmos, self.sim.config.p_tel)
            print("KL2V Done!")
            self.KL2V = KL2V
            return KL2V, 0
        else:
            return self.KL2V, 0

    def compute_Btt2(self):
        IF = self.sim.rtc.get_IFsparse(1).T
        N = IF.shape[0]
        n = IF.shape[1]
        #T = IF[:,-2:].copy()
        T = self.sim.rtc.get_IFtt(1)
        #IF = IF[:,:n-2]
        n = IF.shape[1]

        delta = IF.T.dot(IF).toarray() / N

        # Tip-tilt + piston
        Tp = np.ones((T.shape[0], T.shape[1] + 1))
        Tp[:, :2] = T.copy()  #.toarray()
        deltaT = IF.T.dot(Tp) / N
        # Tip tilt projection on the pzt dm
        tau = np.linalg.inv(delta).dot(deltaT)

        # Famille generatrice sans tip tilt
        G = np.identity(n)
        tdt = tau.T.dot(delta).dot(tau)
        subTT = tau.dot(np.linalg.inv(tdt)).dot(tau.T).dot(delta)
        G -= subTT

        # Base orthonormee sans TT
        gdg = G.T.dot(delta).dot(G)
        U, s, V = np.linalg.svd(gdg)
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
        #IF = rtc.get_IFsparse(1).T
        delta[:-2, :-2] = IF.T.dot(IF).toarray() / N
        delta[-2:, -2:] = T.T.dot(T) / N
        P = Btt.T.dot(delta)

        return Btt.astype(np.float32), P.astype(np.float32)

    def getAi(self):
        return self.wpyr.ai

    def updateSRSE(self, SRSE):
        self.ui.wao_strehlSE.setText(SRSE)

    def updateSRLE(self, SRLE):
        self.ui.wao_strehlLE.setText(SRLE)

    def updateCurrentLoopFrequency(self, freq):
        self.ui.wao_currentFreq.setValue(freq)

    def doImatModal(self, ampliVec, KL2V, Nslopes, noise=False, nmodesMax=0,
                    withTurbu=False, pushPull=False):
        iMatKL = np.zeros((KL2V.shape[1], Nslopes))
        self.aoLoopClicked(False)
        self.ui.wao_run.setChecked(False)
        time.sleep(1)
        st = time.time()
        #currentVolts = self.sim.rtc.get_voltage(0).copy()[None,:]

        if (nmodesMax):
            KLMax = nmodesMax
        else:
            KLMax = KL2V.shape[1]
        for kl in range(KLMax):
            v = ampliVec[kl] * KL2V[:, kl:kl + 1].T.copy()
            if ((pushPull is True) or
                (withTurbu is True)):  # with turbulence/aberrations => push/pull
                self.sim.rtc.set_perturbcom(
                        0, v)  # Adding Perturbation voltage on current iteration
                devpos = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                self.sim.rtc.set_perturbcom(0, -v)
                devmin = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                iMatKL[kl, :] = (devpos - devmin) / (2. * ampliVec[kl])
                #imat[:-2, :] /= pushDMMic
                #if(nmodesMax == 0):# i.e we measured all modes including TT
                #imat[-2:, :] /= pushTTArcsec
            else:  # No turbulence => push only
                self.sim.rtc.set_openloop(0, 1)  # openLoop
                self.sim.rtc.set_perturbcom(0, v)
                iMatKL[kl, :] = self.applyVoltGetSlopes(noise=noise) / ampliVec[kl]
            print("Doing KL interaction matrix on mode: #%d\r" % kl, end=' ')
            os.sys.stdout.flush()

        print("Modal interaction matrix done in %3.0f seconds" % (time.time() - st))
        self.aoLoopClicked(True)
        self.ui.wao_run.setChecked(True)

        return iMatKL

    def doImatPhase(self, cubePhase, Nslopes, noise=False, nmodesMax=0, withTurbu=False,
                    pushPull=False, wfsnum=0):
        iMatPhase = np.zeros((cubePhase.shape[0], Nslopes))
        self.aoLoopClicked(False)
        self.ui.wao_run.setChecked(False)
        time.sleep(1)
        st = time.time()

        for nphase in range(cubePhase.shape[0]):
            if ((pushPull is True) or
                (withTurbu is True)):  # with turbulence/aberrations => push/pull
                self.setNcpaWfs(cubePhase[nphase, :, :], wfsnum=wfsnum)
                devpos = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                self.setNcpaWfs(-cubePhase[nphase, :, :], wfsnum=wfsnum)
                devmin = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                iMatPhase[nphase, :] = (devpos - devmin) / 2
            else:  # No turbulence => push only
                self.sim.rtc.set_openloop(0, 1)  # openLoop
                self.setNcpaWfs(cubePhase[nphase, :, :], wfsnum=wfsnum)
                iMatPhase[nphase, :] = self.applyVoltGetSlopes(noise=noise)
            print("Doing Phase interaction matrix on mode: #%d\r" % nphase, end=' ')
            os.sys.stdout.flush()
        self.setNcpaWfs(cubePhase[nphase, :, :] * 0.,
                        wfsnum=wfsnum)  # Remove the Phase on WFS
        _ = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
        print("Modal interaction matrix done in %3.0f seconds" % (time.time() - st))
        self.aoLoopClicked(True)
        self.ui.wao_run.setChecked(True)

        return iMatPhase

    def computeSlopes(self):
        for w in range(len(self.sim.config.p_wfss)):
            self.sim.wfs.comp_img(w)
        self.sim.rtc.do_centroids(0)
        return self.sim.rtc.get_centroids(0)

    def applyVoltGetSlopes(self, noise=False, turbu=False, reset=1):
        self.sim.rtc.apply_control(0, self.sim.dms)
        for w in range(len(self.sim.config.p_wfss)):
            if (turbu):
                self.sim.wfs.raytrace(w, b"all", self.sim.tel, self.sim.atm,
                                      self.sim.dms, rst=reset, ncpa=1)
            else:
                self.sim.wfs.raytrace(w, b"dm", self.sim.tel, self.sim.atm, self.sim.dms,
                                      rst=reset, ncpa=1)
            self.sim.wfs.comp_img(w, noise=noise)
        self.sim.rtc.do_centroids(0)
        return self.sim.rtc.get_centroids(0)

    def computeModalResiduals(self):
        self.sim.rtc.do_control_geo(1, self.sim.dms, self.sim.tar, 0)
        v = self.sim.rtc.get_com(
                1
        )  # We compute here the residual phase on the DM modes. Gives the Equivalent volts to apply/
        if (self.P is None):
            self.modalBasis, self.P = self.getModes2VBasis("Btt")
        ai = self.P.dot(v) * 1000.  # np rms units
        return ai

    def writeConfigOnFile(self,
                          filepath=os.environ["SHESHA_ROOT"] + "/widgets/canapass.conf"):
        aodict = OrderedDict()

        aodict.update({"Fe": 1 / self.sim.config.p_loop.ittime})
        aodict.update({"teldiam": self.sim.config.p_tel.diam})
        aodict.update({"telobs": self.sim.config.p_tel.cobs})

        # WFS
        aodict.update({"nbWfs": len(self.sim.config.p_wfss)})
        aodict.update({"nbTargets": int(self.sim.config.p_target.ntargets)})
        aodict.update({"nbCam": aodict["nbWfs"]})
        aodict.update({"nbOffaxis": 0})
        aodict.update({"nbNgsWFS": 1})
        aodict.update({"nbLgsWFS": 0})
        aodict.update({"nbFigSensor": 0})
        aodict.update({"nbSkyWfs": aodict["nbWfs"]})
        aodict.update({"nbOffNgs": 0})

        # DMS
        aodict.update({"nbDms": len(self.sim.config.p_dms)})
        aodict.update({"Nactu": sum(self.sim.config.p_controllers[0].nactu)})

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
                    self.sim.config.p_wfss[i]._isvalid))  # Valid subap array
            new_hduwfsl[i].header["DATATYPE"] = "valid_wfs%d" % i

            xytab = np.zeros((self.sim.config.p_wfss[i]._validsubsx.shape[0],
                              self.sim.config.p_wfss[i]._validsubsy.shape[0]))
            xytab[0] = self.sim.config.p_wfss[i]._validsubsx
            xytab[1] = self.sim.config.p_wfss[i]._validsubsy

            new_hduwfsSubapXY.append(
                    pfits.ImageHDU(xytab))  # Valid subap array inXx Y on the detector
            new_hduwfsSubapXY[i].header["DATATYPE"] = "validXY_wfs%d" % i

            pixsize.append(self.sim.config.p_wfss[i].pixsize)
            NslopesList.append(self.sim.config.p_wfss[i]._nvalid * 2)  # slopes per wfs
            NsubapList.append(self.sim.config.p_wfss[i]._nvalid)  # subap per wfs
            listWfsType.append(self.sim.config.p_wfss[i].type)
            xPosList.append(self.sim.config.p_wfss[i].xpos)
            yPosList.append(self.sim.config.p_wfss[i].ypos)
            fstopsize.append(self.sim.config.p_wfss[i].fssize)
            fstoptype.append(self.sim.config.p_wfss[i].fstop)
            nxsubList.append(self.sim.config.p_wfss[i].nxsub)
            nysubList.append(self.sim.config.p_wfss[i].nxsub)
            lambdaList.append(self.sim.config.p_wfss[i].Lambda)
            dms_seen.append(list(self.sim.config.p_wfss[i].dms_seen))
            noise.append(self.sim.config.p_wfss[i].noise)

            if (self.sim.config.p_wfss[i].type == b"pyrhr"):
                pyrModulationList.append(self.sim.config.p_wfss[i].pyr_ampl)
                pyr_npts.append(self.sim.config.p_wfss[i].pyr_npts)
                pyr_pupsep.append(self.sim.config.p_wfss[i].pyr_pup_sep)
                npixPerSub.append(1)
            else:
                pyrModulationList.append(0)
                pyr_npts.append(0)
                pyr_pupsep.append(0)
                npixPerSub.append(self.sim.config.p_wfss[i].npix)
        confname = filepath.split("/")[-1].split('.conf')[0]
        new_hduwfsl.writeto(filepath.split(".conf")[0] + '_wfsConfig.fits', clobber=True)
        new_hduwfsSubapXY.writeto(
                filepath.split(".conf")[0] + '_wfsValidXYConfig.fits', clobber=True)
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
            listDmsType.append(self.sim.config.p_dms[j].type)
            NactuX.append(self.sim.config.p_dms[j].nact)
            unitPerVolt.append(self.sim.config.p_dms[j].unitpervolt)
            push4imat.append(self.sim.config.p_dms[j].push4imat)
            coupling.append(self.sim.config.p_dms[j].coupling)
            tmp = []
            if (self.sim.config.p_dms[j].type != 'tt'):
                tmpdata = np.zeros((2, len(self.sim.config.p_dm0._i1)))
                tmpdata[0, :] = self.sim.config.p_dm0._j1
                tmpdata[1, :] = self.sim.config.p_dm0._i1
                new_hdudmsl.append(pfits.ImageHDU(tmpdata))  # Valid subap array
                new_hdudmsl[j].header["DATATYPE"] = "valid_dm%d" % j
            #for k in range(aodict["nbWfs"]):
            #    tmp.append(self.sim.computeDMrange(j, k))

            push4iMatArcSec.append(tmp)
        new_hdudmsl.writeto(filepath.split(".conf")[0] + '_dmsConfig.fits', clobber=True)
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
            listTargetsLambda.append(self.sim.config.p_target.Lambda[k])
            listTargetsXpos.append(self.sim.config.p_target.xpos[k])
            listTargetsYpos.append(self.sim.config.p_target.ypos[k])
            listTargetsMag.append(self.sim.config.p_target.mag[k])
            listTargetsDmsSeen.append(self.sim.config.p_target.dms_seen[k])

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

    def setIntegratorLaw(self):
        self.sim.rtc.set_commandlaw(0, b"integrator")

    def setDecayFactor(self, decay):
        self.sim.rtc.set_decayFactor(0, decay.astype(np.float32).copy())

    def setEMatrix(self, eMat):
        self.sim.rtc.set_matE(0, eMat.astype(np.float32).copy())

    def doRefslopes(self):
        self.sim.rtc.do_centroids_ref(0)
        print("refslopes done")

    def loadRefSlopes(self, ref):
        self.sim.rtc.set_centroids_ref(0, ref)

    def resetRefslopes(self):
        self.sim.rtc.set_centroids_ref(0, self.getSlopes() * 0.)

    def setPyrModulation(self, pyrmod):
        self.sim.rtc.set_pyr_ampl(0, pyrmod, self.sim.config.p_wfss,
                                  self.sim.config.p_tel)
        print("PYR modulation set to: %f L/D" % pyrmod)
        self.sim.rtc.do_centroids(0)  # To be ready for the next getSlopes

    def setPyrMethod(self, pyrmethod):
        self.sim.rtc.set_pyr_method(0, pyrmethod,
                                    self.sim.config.p_centroiders)  # Sets the pyr method
        print("PYR method set to: %d" % self.sim.rtc.get_pyr_method(0))
        self.sim.rtc.do_centroids(0)  # To be ready for the next getSlopes

    def setNoise(self, noise, numwfs=0):
        self.sim.wfs.set_noise(numwfs, noise)
        print("Noise set to: %d" % noise)

    def setNcpaWfs(self, ncpa, wfsnum):
        self.sim.wfs.set_ncpa_phase(wfsnum, ncpa.astype(np.float32).copy())

    def setNcpaTar(self, ncpa, tarnum):
        self.sim.tar.set_ncpa_phase(tarnum, ncpa.astype(np.float32).copy())

    def set_phaseWFS(self, numwfs, phase):
        pph = phase.astype(np.float32)
        self.sim.wfs.set_phase(0, pph)
        _ = self.computeSlopes()

    def getSlopesGeom(self, nb):
        self.sim.rtc.do_centroids_geom(0)
        slopesGeom = self.sim.rtc.get_centroids(0)
        self.sim.rtc.do_centroids(0)
        return slopesGeom

    def getIpupil(self):
        return self.sim.config.p_geom._ipupil

    def getSpupil(self):
        return self.sim.config.p_geom._spupil

    def getTarImage(self, tarnum, tartype):
        return self.sim.tar.get_image(tarnum, tartype)

    def getMpupil(self):
        return self.sim.config.p_geom._mpupil

    def getAmplipup(self, tarnum):
        return self.sim.config.tar.get_amplipup(tarnum)

    def getPhase(self, tarnum):
        return self.sim.tar.get_phase(tarnum)

    def getWFSPhase(self, wfsnum):
        return self.sim.wfs.get_phase(wfsnum)

    def getTargetPhase(self, tarnum):
        pup = self.getSpupil()
        ph = self.sim.tar.get_phase(tarnum) * pup
        return ph

    def getNcpaWfs(self, wfsnum):
        return self.sim.wfs.get_ncpa_phase(wfsnum)

    def getNcpaTar(self, tarnum):
        return self.sim.tar.get_ncpa_phase(tarnum)

    #def getVolts(self):
    #    return self.sim.rtc.get_voltage(0)
    #
    #def getSlopes(self):
    #     return self.sim.rtc.get_centroids(0)

    def recordCB(self, CBcount, tarnum=0, restart=False, seeAtmos=True):
        self.aoLoopClicked(False)
        self.ui.wao_run.setChecked(False)
        time.sleep(1)

        slopesdata = None
        voltsdata = None
        tarPhaseData = None
        aiData = None
        # Resets the target so that the PSF LE is synchro with the data
        for i in range(self.sim.config.p_target.ntargets):
            self.sim.tar.reset_strehl(i)

        # Starting CB loop...
        for j in tqdm(range(CBcount)):
            self.sim.next(see_atmos=seeAtmos)
            for t in range(self.sim.config.p_target.ntargets):
                self.sim.tar.comp_image(t)

            aiVector = self.computeModalResiduals()
            if (aiData is None):
                aiData = np.zeros((len(aiVector), CBcount))
            aiData[:, j] = aiVector

            slopesVector = self.sim.rtc.get_centroids(0)
            if (slopesdata is None):
                slopesdata = np.zeros((len(slopesVector), CBcount))
            slopesdata[:, j] = slopesVector

            voltsVector = self.sim.rtc.get_com(0)
            if (voltsdata is None):
                voltsdata = np.zeros((len(voltsVector), CBcount))
            voltsdata[:, j] = voltsVector

            tarPhaseArray = self.sim.tar.get_phase(
                    tarnum) * self.sim.config.p_geom._spupil
            if (tarPhaseData is None):
                tarPhaseData = np.zeros((*tarPhaseArray.shape, CBcount))
            tarPhaseData[:, :, j] = tarPhaseArray

        psfLE = self.sim.tar.get_image(tarnum, b"le")
        if (restart):
            self.aoLoopClicked(True)
            self.ui.wao_run.setChecked(True)
        return slopesdata, voltsdata, tarPhaseData, aiData, psfLE

        #wao.sim.config.p_geom._ipupil
        """
        plt.matshow(wao.sim.config.p_geom._ipupil, fignum=1)

        # DM positions in iPupil:
        dmposx = wao.sim.config.p_dm0._xpos
        dmposy = wao.sim.config.p_dm0._ypos
        plt.scatter(dmposx, dmposy, color="blue")




        #WFS position in ipupil
        ipup = wao.sim.config.p_geom._ipupil
        spup = wao.sim.config.p_geom._spupil
        s2ipup = (ipup.shape[0] - spup.shape[0]) / 2.
        posx = wao.sim.config.p_wfss[0]._istart + s2ipup
        posx = posx *  wao.sim.config.p_wfss[0]._isvalid
        posx = posx[np.where(posx > 0)] - ipup.shape[0] / 2 - 1
        posy = wao.sim.config.p_wfss[0]._jstart + s2ipup
        posy = posy * wao.sim.config.p_wfss[0]._isvalid
        posy = posy.T[np.where(posy > 0)] - ipup.shape[0] / 2 - 1

        #center of ssp position in ipupil
        demissp = (posx[1]-posx[0])/2
        sspx = posx+ipup.shape[0]/2+demissp
        sspy = posy+ipup.shape[0]/2+demissp
        plt.scatter(sspx, sspy, color="red")
        imat = wao.sim.rtc.get_imat(0)

        influDM = wao.sim.dms.get_influ(b"pzt", 0)
        influTT = wao.sim.dms.get_influ(b"tt", 0)

        """


try:

    class PyroServer(Thread):
        """
        Pyro object Thread

        """

        def __init__(self, obj):
            Thread.__init__(self)
            self.setDaemon(1)
            self.ready = False
            self.object = obj
            print("initThread")

        def run(self):
            print("Starting Pyro Server...")
            daemon = Pyro4.Daemon()
            ns = Pyro4.locateNS()
            self.ready = True
            try:
                p = Popen("whoami", shell=True, stdout=PIPE, stderr=PIPE)
                out, err = p.communicate()
                if (err != b''):
                    print(err)
                    raise ValueError("ERROR CANNOT RECOGNIZE USER")
                else:
                    user = out.split(b"\n")[0].decode("utf-8")
                    print("User is " + user)
                    ns.remove("waoconfig_" + user)
            except:
                # ns.deleteGroup(':GS')
                # ns.createGroup(":GS")
                pass
            # print self.object.getVar()
            uri = daemon.register(self.object)
            ns.register("waoconfig_" + user, uri)
            print("starting daemon")
            daemon.requestLoop()
            print("daemon started")
except:
    print("Error while initializing Pyro object")
    pass

if __name__ == '__main__':
    arguments = docopt(__doc__)
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('cleanlooks')
    wao = widgetAOWindowPyro(arguments["<parameters_filename>"], BRAMA=True,
                             expert=arguments["--expert"])
    wao.show()
