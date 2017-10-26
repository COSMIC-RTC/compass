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

import pyqtgraph as pg
sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/AOlib")
sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/src/shesha_util")
from tools import plsh, plpyr

import threading

from PyQt5 import QtGui, QtWidgets
from PyQt5.uic import loadUiType
from PyQt5.QtCore import QThread, QObject, QTimer, pyqtSignal
from threading import Thread

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
from aoCalib import adoptCalib_class

sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/data/par/")
WindowTemplate, TemplateBaseClass = loadUiType(
        os.environ["SHESHA_ROOT"] + "/widgets/widget_ao.ui")  # type: type, type
from widget_ao import widgetAOWindow


@Pyro4.expose
class widgetAOWindowPyro(widgetAOWindow):

    def __init__(self, configFile: Any=None, BRAMA: bool=False,
                 expert: bool=False) -> None:
        widgetAOWindow.__init__(self, configFile, BRAMA)
        #Pyro.core.ObjBase.__init__(self)

        self.aoCalib = None

        #############################################################
        #                 CONNECTED BUTTONS                         #
        #############################################################
        # Default path for config files

        self.ui.actionShow_Pyramid_Tools.toggled.connect(self.showPyrTools)
        self.ph2modes = None
        self.KL2V = None
        self.P = None

        #############################################################
        #                       METHODS                             #
        #############################################################

    def loadConfig(self) -> None:
        super().loadConfig()
        print("switching to a generic controller")
        self.sim.config.p_controllers[0].type = scons.ControllerType.GENERIC

    def initPyrTools(self):
        ADOPTPATH = os.getenv("ADOPTPATH")
        sys.path.append(ADOPTPATH + "/widgets")
        from pyrStats import widget_pyrStats
        print("OK Pyramid Tools Widget initialized")
        self.wpyr = widget_pyrStats()
        self.wpyr.show()

    def setPyrToolsParams(self, ai):
        #if(self.ph2KL is None):
        #self.computePh2KL()
        self.wpyr.pup = self.getSpupil()
        self.wpyr.phase = self.getTargetPhase(0)
        self.wpyr.updateResiduals(ai)
        self.wpyr.ph2modes = self.ph2modes

    def showPyrTools(self):
        if (self.wpyr is None):
            try:
                self.initPyrTools()
            except:
                raise ValueError("ERROR: ADOPT  not found. Cannot launch Pyramid tools")
        else:
            if (self.ui.actionShow_Pyramid_Tools.isChecked()):
                self.wpyr.show()
            else:
                self.wpyr.hide()

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

    def getSlopesGeom(self, nb):
        self.sim.rtc.do_centroids_geom(0)
        slopesGeom = self.sim.rtc.get_centroids(0)
        self.sim.rtc.do_centroids(0)
        return slopesGeom

    def getConfig(self, path):
        return self.aoCalib.getConfig(self, path)
        #return cf.returnConfigfromWao(self, filepath=path)

    def getIpupil(self):
        return self.sim.config.p_geom.get_ipupil()

    def getSpupil(self):
        return self.sim.config.p_geom.get_spupil()

    def getImage2(self, tarnum, tartype):
        return self.sim.tar.get_image(tarnum, tartype)

    def getMpupil(self):
        return self.sim.config.p_geom.get_mpupil()

    def getAmplipup(self, tarnum):
        return self.sim.config.tar.get_amplipup(tarnum)

    def getPhase(self, tarnum):
        return self.sim.tar.get_phase(tarnum)

    def getWFSPhase(self, wfsnum):
        return self.sim.wfs.get_phase(wfsnum)

    def computePh2KL(self):
        self.aoLoopClicked(False)
        self.ui.wao_run.setChecked(False)
        oldnoise = self.sim.config.p_wfs0.noise
        self.setNoise(-1)

        if (self.KL2V is None):
            self.KL2V, _ = self.returnkl2V()
        nbmode = self.KL2V.shape[1]
        pup = self.getSpupil()
        ph = self.sim.tar.get_phase(0)
        ph2KL = np.zeros((nbmode, ph.shape[0], ph.shape[1]))
        S = np.sum(pup)
        for i in range(nbmode):
            self.sim.tar.reset_phase(0)
            self.sim.dms.set_full_comm((self.KL2V[:, i]).astype(np.float32).copy())
            self.sim.tar.dmtrace(0, self.sim.dms)
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

    def getTargetPhase(self, tarnum):
        pup = self.getSpupil()
        ph = self.sim.tar.get_phase(tarnum) * pup
        return ph

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

    def getModes2VBasis(self, ModalBasisType):
        if (ModalBasisType == "KL2V"):
            print("Computing KL2V basis...")
            return self.returnkl2V()
        elif (ModalBasisType == "Btt"):
            print("Computing Btt basis...")
            return self.compute_Btt2()

    def computeModalResiduals(self):
        self.sim.rtc.do_control_geo(1, self.sim.dms, self.sim.tar, 0)
        #self.sim.rtc.do_control_geo_on(1, self.sim.dms,self.sim.tar, 0)
        v = self.sim.rtc.getCom(1)
        if (self.P is None):
            self.modalBasis, self.P = self.getModes2VBasis("Btt")
        ai = self.P.dot(v) * 1000.  # np rms units
        return ai

    def returnTest(self, dataType):
        if (dataType == "float"):
            return 1.0
        elif (dataType == "list"):
            return [0, 1, "5"]
        elif (dataType == "array"):
            return np.ones((10, 10))
        elif (dataType == "arraylist"):
            return list(np.ones((10, 10)))

        else:
            print("oups")

    def returnkl2V(self):
        """
        KL2V = ao.compute_KL2V(self.sim.config.p_controllers[
                               0], self.sim.dms, self.sim.config.p_dms, self.sim.config.p_geom, self.sim.config.p_atmos, self.sim.config.p_tel)
        """
        if (self.KL2V is None):
            print("1")
            KL2V = ao.compute_KL2V(self.sim.config.p_controllers[0], self.sim.dms,
                                   self.sim.config.p_dms, self.sim.config.p_geom,
                                   self.sim.config.p_atmos, self.sim.config.p_tel)
            print("ici")
            self.KL2V = KL2V
            return KL2V, 0
        else:
            return self.KL2V, 0

    def _getNcpaWfs(self, wfsnum):
        return self.sim.wfs.get_ncpa_phase(wfsnum)

    def _getNcpaTar(self, tarnum):
        return self.sim.tar.get_ncpa_phase(tarnum)

    def _setNcpaWfs(self, ncpa, wfsnum):
        self.sim.wfs.set_ncpa_phase(wfsnum, ncpa.astype(np.float32).copy())

    def _setNcpaTar(self, ncpa, tarnum):
        self.sim.tar.set_ncpa_phase(tarnum, ncpa.astype(np.float32).copy())

    def doRefslopes(self):
        self.sim.rtc.do_centroids_ref(0)
        print("refslopes done")

    def loadRefSlopes(self, ref):
        self.sim.rtc.set_centroids_ref(0, ref)

    def resetRefslopes(self):
        self.sim.rtc.set_centroids_ref(0, self.getSlopes() * 0.)

    def setCommandMatrix(self, cMat):
        return self.sim.rtc.set_cmat(0, cMat)

    def setPertuVoltages(self, actusVolts):
        self.sim.rtc.setPertuVoltages(0, actusVolts.astype(np.float32).copy())

    def getVolts(self):
        return self.sim.rtc.get_voltage(0)

    def getSlopes(self):
        return self.sim.rtc.get_centroids(0)

    def setIntegratorLaw(self):
        self.sim.rtc.set_commandlaw(0, "integrator")

    def setGain(self, gain):
        self.sim.rtc.set_gain(0, gain)

    def setDecayFactor(self, decay):
        self.sim.rtc.set_decayFactor(0, decay.astype(np.float32).copy())

    def setEMatrix(self, eMat):
        self.sim.rtc.set_matE(0, eMat.astype(np.float32).copy())

    def closeLoop(self):
        self.sim.rtc.set_openloop(0, 0)

    def openLoop(self):
        self.sim.rtc.set_openloop(0, 1)

    def updateSRSE(self, SRSE):
        self.ui.wao_strehlSE.setText(SRSE)

    def updateSRLE(self, SRLE):
        self.ui.wao_strehlLE.setText(SRLE)

    def set_phaseWFS(self, numwfs, phase):
        pph = phase.astype(np.float32)
        self.sim.wfs.set_phase(0, pph)
        _ = self.computeSlopes()

    def updateCurrentLoopFrequency(self, freq):
        self.ui.wao_currentFreq.setValue(freq)

    def InitConfig(self):
        widgetAOWindow.InitConfig(self)
        try:
            ps = PyroServer(wao)
            ps.start()
        except:
            print("Warning: Error while starting Pyro server")

    def InitConfigFinished(self) -> None:
        widgetAOWindow.InitConfigFinished(self)
        self.aoCalib = adoptCalib_class(self.sim.config, self.sim.wfs, self.sim.tel,
                                        self.sim.atm, self.sim.dms, self.sim.tar,
                                        self.sim.rtc, ao)

    def measureIMatKL(self, ampliVec, KL2V, Nslopes, noise=False, nmodesMax=0,
                      withTurbu=False, pushPull=False):
        iMatKL = np.zeros((KL2V.shape[1], Nslopes))
        self.aoLoopClicked(False)
        self.ui.wao_run.setChecked(False)
        time.sleep(1)
        st = time.time()
        currentVolts = self.sim.rtc.get_voltage(0) * 0.

        if (nmodesMax):
            KLMax = nmodesMax
        else:
            KLMax = KL2V.shape[1]
        for kl in range(KLMax):
            v = ampliVec[kl] * KL2V[:, kl:kl + 1].T.copy()
            if ((pushPull is True) or
                (withTurbu is True)):  # with turbulence/aberrations => push/pull
                self.sim.rtc.set_perturbcom(0, v + currentVolts)
                devpos = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                self.sim.rtc.set_perturbcom(0, -v + currentVolts)
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

    def hello(self):
        return "Hello!"

    def computeSlopes(self):
        for w in range(len(self.sim.config.p_wfss)):
            self.sim.wfs.comp_img(w)
        self.sim.rtc.do_centroids(0)
        return self.sim.rtc.get_centroids(0)

    def applyVoltGetSlopes(self, noise=False, turbu=False):
        self.sim.rtc.apply_control(0, self.sim.dms)
        for w in range(len(self.sim.config.p_wfss)):
            if (turbu):
                self.sim.wfs.raytrace(w, b"all", self.sim.tel, self.sim.atm,
                                      self.sim.dms, rst=1, ncpa=1)
            else:
                self.sim.wfs.raytrace(w, b"dm", self.sim.tel, self.sim.atm, self.sim.dms,
                                      rst=1, ncpa=1)
            self.sim.wfs.comp_img(w, noise=noise)
        self.sim.rtc.do_centroids(0)
        return self.sim.rtc.get_centroids(0)

    def computeModalResiduals(self):
        self.sim.rtc.do_control_geo(1, self.sim.dms, self.sim.tar, 0)
        #self.sim.rtc.do_control_geo_on(1, self.sim.dms,self.sim.tar, 0)
        v = self.sim.rtc.getCom(1)
        if (self.P is None):
            self.modalBasis, self.P = self.getModes2VBasis("Btt")
        ai = self.P.dot(v) * 1000.  # np rms units
        return ai

    def loopOnce(self) -> None:
        widgetAOWindow.loopOnce(self)
        start = time.time()
        self.sim.next(see_atmos=self.see_atmos)

        refreshDisplayTime = 1. / self.ui.wao_frameRate.value()

        if (time.time() - self.refreshTime > refreshDisplayTime):
            if (self.ui.actionShow_Pyramid_Tools.isChecked()):
                ai = self.computeModalResiduals()
                self.setPyrToolsParams(ai)


try:

    class PyroServer(Thread):
        """
        Main Geometry Server Thread

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
                ns.unregister("waoconfig")
            except:
                # ns.deleteGroup(':GS')
                # ns.createGroup(":GS")
                pass
            # print self.object.getVar()
            uri = daemon.register(self.object)
            ns.register("waoconfig", uri)
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
