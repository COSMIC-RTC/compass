"""
widget_canapass.

import cProfile
import pstats as ps
"""

import sys
import os
import numpy as np
import naga as ch
import shesha as ao
import time
import matplotlib.pyplot as plt
import pyqtgraph as pg
import tools
import hdf5_utils as h5u
import threading
from PyQt4.uic import loadUiType
from PyQt4 import QtGui
from PyQt4.Qt import QThread, QObject
from PyQt4.QtCore import QTimer, SIGNAL
from functools import partial
from pandas import HDFStore, read_hdf
from threading import Thread
try:
    import qdarkstyle
    darkStyle=True
except:
    darkStyle=False
import astropy.io.fits as pfits

sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/widgets/")
sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/data/par/")
sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/lib/")
import compassConfigToFile as cf
from aoCalib import adoptCalib_class

WindowTemplate, TemplateBaseClass = loadUiType(
    os.environ["SHESHA_ROOT"] + "/widgets/widget_canapass.ui")

import compassConfigToFile as cf
plt.rcParams['image.cmap'] = 'viridis'

"""
low levels debugs:
gdb --args python -i widget_canapass.py

"""

PYROVERSION = int(os.getenv("PYROVERSION"))
if PYROVERSION == None:
    raise Exception('PYROVERSION environment variable is unset!')

if(PYROVERSION == 4):
    try:
        import Pyro4
        USE_PYRO = True
    except:
        USE_PYRO = False
        raise Exception("ERROR could not import Pyro4 (is it installed?)!!!!")

if(PYROVERSION == 3):
    try:
        import Pyro.core
        USE_PYRO = True

    except:
        USE_PYRO = False
        raise Exception("ERROR could not import Pyro (v3) (is it installed?)!!!!")


@Pyro4.expose
class widgetAOWindow(TemplateBaseClass):
    def __init__(self):
        TemplateBaseClass.__init__(self)

        self.SRLE = []
        self.SRSE = []
        self.numiter = []

        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        self.pyroVersion = PYROVERSION
        #############################################################
        #                   ATTRIBUTES                              #
        #############################################################
        self.c = None
        self.atm = None
        self.tel = None
        self.wfs = None
        self.rtc = None
        self.tar = None
        self.dms = None

        self.config = None
        # self.displayLock = threading.Lock()
        self.loopLock = threading.Lock()
        self.iter = 0
        self.loaded = False
        self.stop = False

        self.refreshTime = time.time()

        self.loop = None
        self.assistant = None
        self.selector_init = None

        self.brama_rtc_flag = 1
        self.brama_tar_flag = 1
        self.see_atmos = 0
        self.aoCalib = None
        #############################################################
        #             PYQTGRAPH WINDOW INIT                         #
        #############################################################

        self.img = pg.ImageItem(border='w')  # create image area
        self.img.setTransform(QtGui.QTransform(
            0, 1, 1, 0, 0, 0))  # flip X and Y
        # self.p1 = self.ui.wao_pgwindow.addPlot() # create pyqtgraph plot area
        self.p1 = self.ui.wao_pgwindow.addViewBox()
        self.p1.setAspectLocked(True)
        self.p1.addItem(self.img)  # Put image in plot area

        self.hist = pg.HistogramLUTItem()  # Create an histogram
        self.hist.setImageItem(self.img)  # Compute histogram from img
        self.ui.wao_pgwindow.addItem(self.hist)
        self.hist.autoHistogramRange()  # init levels
        self.hist.setMaximumWidth(100)

        #############################################################
        #                  CONNECTED BUTTONS                        #
        #############################################################
        # Default path for config files
        self.defaultParPath = os.environ[
            "SHESHA_ROOT"] + "/data/par/MICADO/"
        self.ui.wao_loadConfig.clicked.connect(self.addConfigFromFile)
        self.ui.wao_init.clicked.connect(self.InitConfig)
        self.ui.wao_run.setCheckable(True)
        self.ui.wao_run.clicked[bool].connect(self.aoLoopClicked)
        self.ui.wao_next.clicked.connect(self.mainLoop)
        self.imgType = str(self.ui.wao_selectScreen.currentText())
        self.ui.wao_unzoom.clicked.connect(self.p1.autoRange)
        self.ui.wao_resetSR.clicked.connect(self.resetSR)
        self.ui.wao_selectScreen.currentIndexChanged.connect(
            partial(self.updateNumberSelector, textType=None))
        self.ui.wao_selectNumber.currentIndexChanged.connect(
            self.setNumberSelection)
        self.ui.wao_rtcWindowMPL.hide()
        self.ui.wao_frameRate.setValue(2)
        #self.ui.wao_PSFlogscale.clicked.connect(self.updateDisplay)
        self.ui.wao_PSFlogscale.toggled.connect(self.updateDisplay)
        self.ui.wao_atmosphere.setCheckable(True)
        self.ui.wao_atmosphere.clicked[bool].connect(self.set_atmos)
        self.dispStatsInTerminal = False
        self.ui.wao_clearSR.clicked.connect(self.clearSR)
        self.ui.wao_Display.setCheckState(True)
        self.ui.wao_Display.stateChanged.connect(self.updateDisplay)
        #self.ui.StatsInTerminal.stateChanged.connect(self.updateStatsInTerminal)
        self.ui.StatsInTerminal.toggled.connect(self.updateStatsInTerminal)
        self.ui.actionQuit.toggled.connect(self.quitGUI)

        self.SRcircleAtmos = {}
        self.SRcircleWFS = {}
        self.SRcircleDM = {}
        self.SRcircleTarget = {}

        self.ui.wao_loadConfig.setDisabled(False)
        self.ui.wao_init.setDisabled(True)
        self.ui.wao_run.setDisabled(True)
        self.ui.wao_next.setDisabled(True)
        self.ui.wao_unzoom.setDisabled(True)
        self.ui.wao_resetSR.setDisabled(True)

        # self.addConfigFromFile(
        #     filepath=os.environ["SHESHA_ROOT"] + "/data/par/canapass.py")
        # self.InitConfig()

    def writefits(self, data, path):
        pfits.writeto(path, data, clobber=True)

    def getConfig(self, path):
        return self.aoCalib.getConfig(self, path)
        #return cf.returnConfigfromWao(self, filepath=path)
    def returnkl2V(self, path):
        print("computing KL2V...")
        KL2V = self._returnkl2V()
        print("KL2V done")
        print(path)
        if(self.pyroVersion == 4):
            print("using V4 ...")
            self.writefits(KL2V, path)
        elif(self.pyroVersion == 3):
            print("using V3 ...")
            return KL2V
        else:
            raise Exception("ERROR pyro version not set")


    def _returnkl2V(self):

        #KL2V = ao.compute_KL2V(self.config.p_controllers[0], self.dms, self.config.p_dms, self.config.p_geom, self.config.p_atmos, self.config.p_tel)
        KL2V = np.ones(256)
        time.sleep(10)
        return KL2V

    def doRefslopes(self):
        self.rtc.do_centroids_ref(0)
        print("refslopes done")

    def setCommandMatrix(self, cMat):
        return self.rtc.set_cmat(0, cMat)

    def setPertuVoltages(self, actusVolts):
        # self.dms.set_full_comm(actusVolts.astype(np.float32).copy())
        self.rtc.setPertuVoltages(0, actusVolts.astype(np.float32).copy())

    def getVolts(self):
        return self.rtc.getVoltage(0)

    def getSlopes(self):
        return self.rtc.getCentroids(0)

    def setIntegratorLaw(self):
        self.rtc.set_commandlaw(0, "integrator")

    def setGain(self, gain):
        self.rtc.set_gain(0, gain)

    def setDecayFactor(self, decay):
        self.rtc.set_decayFactor(0, decay.astype(np.float32).copy())

    def setEMatrix(self, eMat):
        self.rtc.set_matE(0, eMat.astype(np.float32).copy())

    def closeLoop(self):
        self.rtc.set_openloop(0, 0)

    def openLoop(self):
        self.rtc.set_openloop(0, 1)

    def updateStatsInTerminal(self):
        if(self.ui.StatsInTerminal.isChecked()):
            self.dispStatsInTerminal = True
        else:
            self.dispStatsInTerminal = False

    def updateSRSE(self, SRSE):
        self.ui.wao_strehlSE.setText(SRSE)

    def updateSRLE(self, SRLE):
        self.ui.wao_strehlLE.setText(SRLE)

    def updateCurrentLoopFrequency(self, freq):
        self.ui.wao_currentFreq.setValue(freq)

    def set_atmos(self, atmos):
        self.see_atmos = atmos

    def resetSR(self):
        for t in range(self.config.p_target.ntargets):
            self.tar.reset_strehl(t)

    def quitGUI(self):
        reply = QtGui.QMessageBox.question(self, 'Message',
                                           "Are you sure to quit?", QtGui.QMessageBox.Yes |
                                           QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            self.stop = True
            if self.loop is not None:
                self.loop.join()
            # super(widgetAOWindow, self).closeEvent(event)
            quit()
            # sys.exit()
        else:
            print("Exit aborted")

    def closeEvent(self, event):
        reply = QtGui.QMessageBox.question(self, 'Message',
                                           "Are you sure to quit?", QtGui.QMessageBox.Yes |
                                           QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
            self.stop = True
            if self.loop is not None:
                self.loop.join()
            # super(widgetAOWindow, self).closeEvent(event)
            quit()
            # sys.exit()
        else:
            event.ignore()

        #############################################################
        #                        METHODS                            #
        #############################################################

    def addConfigFromFile(self, filepath=False):
        if filepath is False:
            filepath = str(QtGui.QFileDialog(directory=self.defaultParPath).getOpenFileName(
                self, "Select parameter file", "", "parameters file (*.py);;hdf5 file (*.h5);;all files (*)"))
        self.loaded = False
        filename = filepath.split('/')[-1]
        if(filepath.split('.')[-1] == "py"):
            pathfile = filepath.split(filename)[0]
            # if (pathfile not in sys.path):
            sys.path.insert(0, pathfile)

            if self.config is not None:
                print("Removing previous config")
                self.config = None
                config = None

            print("loading ", filename.split(".py")[0])
            exec("import %s as config" % filename.split(".py")[0])

            sys.path.remove(pathfile)
        else:
            raise ValueError("Parameter file extension must be .py or .h5")

        print("switching to generic Controller")
        config.p_controller0.set_type("generic")

        self.config = config
        self.ui.wao_selectScreen.clear()
        if(self.config.p_wfss[0].type_wfs == "pyrhr"):
            self.selector_init = ["Phase - Atmos", "Phase - WFS", "Pyrimg - LR",
                                  "Pyrimg - HR", "Centroids - WFS", "Slopes - WFS",
                                  "Phase - Target", "Phase - DM",
                                  "PSF LE", "PSF SE"]
        else:
            self.selector_init = ["Phase - Atmos", "Phase - WFS", "Spots - WFS",
                                  "Centroids - WFS", "Slopes - WFS",
                                  "Phase - Target", "Phase - DM",
                                  "PSF LE", "PSF SE"]
        self.ui.wao_selectScreen.addItems(self.selector_init)
        self.ui.wao_selectScreen.setCurrentIndex(0)
        self.updateNumberSelector(textType=self.imgType)

        self.ui.wao_init.setDisabled(False)
        self.ui.wao_run.setDisabled(True)
        self.ui.wao_next.setDisabled(True)
        self.ui.wao_unzoom.setDisabled(True)
        self.ui.wao_resetSR.setDisabled(True)

    def aoLoopClicked(self, pressed):
        if(pressed):
            self.c.set_activeDeviceForce(0, 1)
            self.stop = False
            self.refreshTime = time.time()
            self.run()
            self.ui.wao_Display.setCheckState(True)
            # self.loop = threading.Thread(target=self.run)
            # self.loop.start()
        else:
            self.stop = True
            self.ui.wao_Display.setCheckState(False)
            # self.loop.join()
            # self.loop = None

    def setNumberSelection(self):
        if(self.ui.wao_selectNumber.currentIndex() > -1):
            self.numberSelected = self.ui.wao_selectNumber.currentIndex()
        else:
            self.numberSelected = 0
        self.updateDisplay()

    def updateNumberSelector(self, textType=None):
        if(textType is None):
            textType = str(self.ui.wao_selectScreen.currentText())
        self.imgType = textType
        self.ui.wao_selectNumber.clear()
        if(textType == "Phase - Atmos"):
            n = self.config.p_atmos.nscreens
        elif(textType == "Phase - WFS" or textType == "Spots - WFS" or textType == "Centroids - WFS" or textType == "Slopes - WFS" or textType == "Pyrimg - HR" or textType == "Pyrimg - LR"):
            n = len(self.config.p_wfss)
        elif(textType == "Phase - Target" or textType == "PSF LE" or textType == "PSF SE"):
            n = self.config.p_target.ntargets
        elif(textType == "Phase - DM"):
            n = len(self.config.p_dms)
        else:
            n = 0
        self.ui.wao_selectNumber.addItems([str(i) for i in range(n)])
        self.updateDisplay()

    def InitConfig(self):
        self.loaded = False
        self.ui.wao_loadConfig.setDisabled(True)
        self.ui.wao_init.setDisabled(True)
        thread = WorkerThread(self, self.InitConfigThread)
        QObject.connect(thread, SIGNAL(
            "jobFinished( PyQt_PyObject )"), self.InitConfigFinished)
        thread.start()

    def InitConfigThread(self):
        if(hasattr(self, "atm")):
            del self.atm
        if(hasattr(self, "tel")):
            del self.tel
        if(hasattr(self, "wfs")):
            del self.wfs
        if(hasattr(self, "rtc")):
            del self.rtc
        if(hasattr(self, "tar")):
            del self.tar
        if(hasattr(self, "dms")):
            del self.dms

        self.iter = 0

        gpudevice = self.ui.wao_deviceNumber.value()  # using GUI value
        gpudevice = self.config.p_loop.devices
        # gpudevice = "ALL"  # using all GPU avalaible
        # gpudevice = np.arange(4, dtype=np.int32)  # using 4 GPUs: 0-3
        gpudevice = np.array([4, 5, 6, 7], dtype=np.int32)  # using 4 GPUs: 4-7
        self.ui.wao_deviceNumber.setDisabled(True)

        print("-> using GPU", gpudevice)

        if not self.c:
            if type(gpudevice) is np.ndarray:
                self.c = ch.naga_context(devices=gpudevice)
            elif type(gpudevice) is int:
                self.c = ch.naga_context(gpudevice)
            else:
                self.c = ch.naga_context()

        self.wfs, self.tel = ao.wfs_init(self.config.p_wfss, self.config.p_atmos, self.config.p_tel,
                                         self.config.p_geom, self.config.p_target, self.config.p_loop,
                                         self.config.p_dms)

        self.atm = ao.atmos_init(self.c, self.config.p_atmos, self.config.p_tel,
                                 self.config.p_geom, self.config.p_loop,
                                 self.config.p_wfss, self.wfs, self.config.p_target, clean=1, load={})

        self.dms = ao.dm_init(
            self.config.p_dms, self.config.p_wfss, self.wfs, self.config.p_geom, self.config.p_tel)

        self.tar = ao.target_init(self.c, self.tel, self.config.p_target, self.config.p_atmos,
                                  self.config.p_geom, self.config.p_tel, self.config.p_dms, brama=self.brama_tar_flag)

        self.rtc = ao.rtc_init(self.tel, self.wfs, self.config.p_wfss, self.dms, self.config.p_dms,
                               self.config.p_geom, self.config.p_rtc, self.config.p_atmos,
                               self.atm, self.config.p_tel, self.config.p_loop, do_refslp=False, clean=1, simul_name="",
                               load={}, brama=self.brama_rtc_flag, g_tar=self.tar)
        self.rtc.set_openloop(0, 1)
        self.loaded = True
        self.aoCalib = adoptCalib_class(self.config, self.wfs, self.tel, self.atm, self.dms, self.tar, self.rtc, ao)

    def InitConfigFinished(self):
        self.ui.wao_loadConfig.setDisabled(False)

        self.currentViewSelected = None
        self.SRCrossX = None
        self.SRCrossY = None

        # remove previous pupil materialisation
#        vb = self.p1.getViewBox()
#        for it in vb.items():
#            if type(it) is pg.ScatterPlotItem:
#                self.p1.removeItem(it)
        for i in self.SRcircleAtmos:
            self.p1.removeItem(self.SRcircleAtmos[i])
        for i in self.SRcircleWFS:
            self.p1.removeItem(self.SRcircleWFS[i])
        for i in self.SRcircleDM:
            self.p1.removeItem(self.SRcircleDM[i])
        for i in self.SRcircleTarget:
            self.p1.removeItem(self.SRcircleTarget[i])

        for i in range(len(self.config.p_atmos.alt)):
            data = self.atm.get_screen(self.config.p_atmos.alt[i])
            cx, cy = self.circleCoords(
                self.config.p_geom.pupdiam / 2, 1000, data.shape[0], data.shape[1])
            self.SRcircleAtmos[i] = pg.ScatterPlotItem(cx, cy, pen='r', size=1)
            self.p1.addItem(self.SRcircleAtmos[i])
            self.SRcircleAtmos[i].setPoints(cx, cy)
            self.SRcircleAtmos[i].hide()

        for i in range(len(self.config.p_wfss)):
            data = self.wfs.get_phase(i)
            cx, cy = self.circleCoords(
                self.config.p_geom.pupdiam / 2, 1000, data.shape[0], data.shape[1])
            self.SRcircleWFS[i] = pg.ScatterPlotItem(cx, cy, pen='r', size=1)
            self.p1.addItem(self.SRcircleWFS[i])
            self.SRcircleWFS[i].setPoints(cx, cy)
            self.SRcircleWFS[i].hide()

        for i in range(len(self.config.p_dms)):
            dm_type = self.config.p_dms[i].type_dm
            alt = self.config.p_dms[i].alt
            data = self.dms.get_dm(dm_type, alt)
            cx, cy = self.circleCoords(
                self.config.p_geom.pupdiam / 2, 1000, data.shape[0], data.shape[1])
            self.SRcircleDM[i] = pg.ScatterPlotItem(cx, cy, pen='r', size=1)
            self.p1.addItem(self.SRcircleDM[i])
            self.SRcircleDM[i].setPoints(cx, cy)
            self.SRcircleDM[i].hide()

        for i in range(self.config.p_target.ntargets):
            data = self.tar.get_phase(i)
            cx, cy = self.circleCoords(
                self.config.p_geom.pupdiam / 2, 1000, data.shape[0], data.shape[1])
            self.SRcircleTarget[i] = pg.ScatterPlotItem(
                cx, cy, pen='r', size=1)
            self.p1.addItem(self.SRcircleTarget[i])
            self.SRcircleTarget[i].setPoints(cx, cy)
            self.SRcircleTarget[i].show()

        print("====================")
        print("init done")
        print("====================")
        print("objects initialized on GPU:")
        print("--------------------------------------------------------")
        print(self.atm)
        print(self.wfs)
        print(self.dms)
        print(self.tar)
        print(self.rtc)

        self.updateDisplay()
        self.p1.autoRange()

        self.ui.wao_init.setDisabled(True)
        self.ui.wao_run.setDisabled(False)
        self.ui.wao_next.setDisabled(False)
        self.ui.wao_unzoom.setDisabled(False)
        self.ui.wao_resetSR.setDisabled(False)

    def circleCoords(self, ampli, npts, datashape0, datashape1):
        # ampli = self.config.p_geom.pupdiam/2
        # npts = 100
        cx = ampli * np.sin((np.arange(npts) + 1) * 2. *
                            np.pi / npts) + datashape0 / 2
        cy = ampli * np.cos((np.arange(npts) + 1) * 2. *
                            np.pi / npts) + datashape1 / 2
        return cx, cy

    def resetDM(self):
        if(self.dms):
            ndm = self.ui.wao_selectDM.currentIndex()
            if(ndm > -1):
                self.dms.resetdm(
                    str(self.ui.wao_dmTypeSelector.currentText()), self.ui.wao_dmAlt.value())
                self.updateDisplay()
                print("DM " + str(ndm) + " reset")
            else:
                print("Invalid DM : please select a DM to reset")
        else:
            print("There is not any dm to reset")

    def computeDMrange(self, numdm, numwfs, push4imat=None, unitpervolt=None):
        i = numdm
        if(push4imat is None or push4imat == 0):
            push4imat = self.config.p_dms[i].push4imat
        if(unitpervolt is None or unitpervolt == 0):
            unitpervolt = self.config.p_dms[i].unitpervolt

        actuPushInMicrons = push4imat * unitpervolt
        coupling = self.config.p_dms[i].coupling
        a = coupling * actuPushInMicrons
        b = 0
        c = actuPushInMicrons
        d = coupling * actuPushInMicrons
        if(self.config.p_dms[i].type_dm is not "tt"):
            dist = self.config.p_tel.diam
        else:
            dist = self.config.p_tel.diam / self.config.p_wfss[numwfs].nxsub
        Delta = (1e-6 * (((c + d) / 2) - ((a + b) / 2)))
        actuPushInArcsecs = 206265. * Delta / dist
        return actuPushInArcsecs

    def setupDisp(self, fig="pg"):
        if fig == "pg":
            widToShow = self.ui.wao_pgwindow
            widToHide = self.ui.wao_rtcWindowMPL
        elif fig == "MPL":
            widToShow = self.ui.wao_rtcWindowMPL
            widToHide = self.ui.wao_pgwindow
        else:
            return

        if(not widToShow.isVisible()):
            widToShow.show()
            widToHide.hide()

    def measureIMatKL(self, ampliVec, KL2V, Nslopes, noise=False, nmodesMax=0, withTurbu=False):
        iMatKL = np.zeros((KL2V.shape[1], Nslopes))
        self.aoLoopClicked(False)
        self.ui.wao_run.setChecked(False)
        time.sleep(1)
        st = time.time()
        currentVolts = wao.rtc.getVoltage(0)*0.

        if(nmodesMax):
            KLMax = nmodesMax
        else:
            KLMax = KL2V.shape[1]
        for kl in range(KLMax):
            v = ampliVec[kl] * KL2V[:, kl:kl + 1].T.copy()
            if(withTurbu): # with turbulence/aberrations => push/pull
                self.rtc.set_perturbcom(0, v + currentVolts)
                devpos = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                self.rtc.set_perturbcom(0, -v + currentVolts)
                devmin = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                iMatKL[kl, :] = (devpos - devmin) / (2. * ampliVec[kl])
                #imat[:-2, :] /= pushDMMic
                #if(nmodesMax == 0):# i.e we measured all modes including TT
                #imat[-2:, :] /= pushTTArcsec
            else: # No turbulence => push only
                self.rtc.set_openloop(0, 1)  # openLoop
                self.rtc.set_perturbcom(0, v)
                iMatKL[kl, :] = self.applyVoltGetSlopes(noise=noise) / ampliVec[kl]
            print("Doing KL interaction matrix on mode: #%d\r" % kl, end=' ')
            os.sys.stdout.flush()

        print("Modal interaction matrix done in %3.0f seconds" % (time.time() - st))
        self.aoLoopClicked(True)
        self.ui.wao_run.setChecked(True)

        return iMatKL


    def applyVoltGetSlopes(self, noise=False, turbu=False):
        self.rtc.applycontrol(0, self.dms)
        for w in range(len(self.config.p_wfss)):
            if(turbu):
                self.wfs.sensors_trace(w, "all", self.tel, self.atm, self.dms, rst=1)
            else:
                self.wfs.sensors_trace(w, "dm", self.tel, self.atm, self.dms, rst=1)
            self.wfs.sensors_compimg(w, noise=noise)
        self.rtc.docentroids(0)
        return self.rtc.getcentroids(0)


    def clearSR(self):
        self.SRLE = [self.SRLE[-1]]
        self.SRSE = [self.SRSE[-1]]
        self.numiter = [self.numiter[-1]]

    def updateSRDisplay(self, SRLE, SRSE, numiter):
        self.SRLE.append(SRLE)
        self.SRSE.append(SRSE)
        self.numiter.append(numiter)
        if(len(self.SRSE) > 100):  # Clipping last 100 points...
            self.SRLE = self.SRLE[-100:]
            self.SRSE = self.SRSE[-100:]
            self.numiter = self.numiter[-100:]
        self.ui.wao_SRPlotWindow.canvas.axes.clear()
        self.ui.wao_SRPlotWindow.canvas.axes.yaxis.set_label("SR")
        self.ui.wao_SRPlotWindow.canvas.axes.xaxis.set_label("num iter")
        self.ui.wao_SRPlotWindow.canvas.axes.plot(
            self.numiter, self.SRSE, linestyle="--", color="red", marker="o", label="SR SE")
        self.ui.wao_SRPlotWindow.canvas.axes.plot(
            self.numiter, self.SRLE, linestyle="--", color="blue", marker="o", label="SR LE")
        # self.ui.wao_SRPlotWindow.canvas.axes.grid()
        self.ui.wao_SRPlotWindow.canvas.draw()

    def updateDisplay(self):
        if (not self.loaded) or (not self.ui.wao_Display.isChecked()):
            # print " widget not fully initialized"
            return

        data = None
        if not self.loopLock.acquire(False):
            # print "Loop locked"
            return
        else:
            try:
                if(self.SRCrossX and (self.imgType in ["Phase - Target", "Phase - DM", "Phase - Atmos", "Phase - WFS", "Spots - WFS", "Centroids - WFS", "Slopes - WFS"])):
                    self.SRCrossX.hide()
                    self.SRCrossY.hide()

                # if(self.SRcircle and (self.imgType in ["Spots - WFS",
                # "Centroids - WFS", "Slopes - WFS","PSF SE","PSF LE"])):
                for i in range(len(self.config.p_atmos.alt)):
                    self.SRcircleAtmos[i].hide()
                for i in range(len(self.config.p_wfss)):
                    self.SRcircleWFS[i].hide()
                for i in range(len(self.config.p_dms)):
                    self.SRcircleDM[i].hide()
                for i in range(self.config.p_target.ntargets):
                    self.SRcircleTarget[i].hide()

                if(self.atm):
                    if(self.imgType == "Phase - Atmos"):
                        self.setupDisp("pg")
                        data = self.atm.get_screen(
                            self.config.p_atmos.alt[self.numberSelected])
                        if(self.imgType != self.currentViewSelected):
                            self.p1.setRange(
                                xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                        self.currentViewSelected = self.imgType
                        self.SRcircleAtmos[self.numberSelected].show()

                if(self.wfs):
                    if(self.imgType == "Phase - WFS"):
                        self.setupDisp("pg")
                        data = self.wfs.get_phase(self.numberSelected)
                        if(self.imgType != self.currentViewSelected):
                            self.p1.setRange(
                                xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                        self.currentViewSelected = self.imgType
                        self.SRcircleWFS[self.numberSelected].show()

                    if(self.imgType == "Spots - WFS"):
                        self.setupDisp("pg")
                        if(self.config.p_wfss[self.numberSelected].type_wfs == "sh"):
                            data = self.wfs.get_binimg(self.numberSelected)
                        elif(self.config.p_wfss[self.numberSelected].type_wfs == "pyr"):
                            data = self.wfs.get_pyrimg(self.numberSelected)
                        if(self.imgType != self.currentViewSelected):
                            self.p1.setRange(
                                xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                        self.currentViewSelected = self.imgType

                    if(self.imgType == "Pyrimg - LR"):
                        self.setupDisp("pg")
                        if(self.config.p_wfss[self.numberSelected].type_wfs == "pyrhr"):
                            data = self.wfs.get_pyrimg(self.numberSelected)
                        if(self.imgType != self.currentViewSelected):
                            self.p1.setRange(
                                xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                        self.currentViewSelected = self.imgType

                    if(self.imgType == "Pyrimg - HR"):
                        self.setupDisp("pg")
                        if(self.config.p_wfss[self.numberSelected].type_wfs == "pyrhr"):
                            data = self.wfs.get_pyrimghr(
                                self.numberSelected)
                        if(self.imgType != self.currentViewSelected):
                            self.p1.setRange(
                                xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                        self.currentViewSelected = self.imgType

                    if(self.imgType == "Centroids - WFS"):
                        self.setupDisp("MPL")
                        self.ui.wao_rtcWindowMPL.canvas.axes.clear()
                        # retrieving centroids
                        centroids = self.rtc.getCentroids(0)
                        nvalid = [
                            2 * o._nvalid for o in self.config.p_wfss]
                        ind = np.sum(nvalid[:self.numberSelected])
                        if(self.config.p_wfss[self.numberSelected].type_wfs == "pyrhr"):
                            tools.plpyr(centroids[ind:ind + nvalid[self.numberSelected]], self.config.p_wfs0._isvalid)
                        else:
                            x, y, vx, vy = tools.plsh(centroids[ind:ind + nvalid[self.numberSelected]], self.config.p_wfss[
                                                  self.numberSelected].nxsub, self.config.p_tel.cobs, returnquiver=True)  # Preparing mesh and vector for display
                        self.ui.wao_rtcWindowMPL.canvas.axes.quiver(
                            x, y, vx, vy, pivot='mid')
                        self.ui.wao_rtcWindowMPL.canvas.draw()
                        self.currentViewSelected = self.imgType

                        return
                    if(self.imgType == "Slopes - WFS"):
                        self.setupDisp("MPL")
                        self.ui.wao_rtcWindowMPL.canvas.axes.clear()
                        self.wfs.slopes_geom(self.numberSelected, 0)
                        slopes = self.wfs.get_slopes(self.numberSelected)
                        x, y, vx, vy = tools.plsh(slopes, self.config.p_wfss[
                                                  self.numberSelected].nxsub, self.config.p_tel.cobs, returnquiver=True)  # Preparing mesh and vector for display
                        self.ui.wao_rtcWindowMPL.canvas.axes.quiver(
                            x, y, vx, vy, pivot='mid')
                        self.ui.wao_rtcWindowMPL.canvas.draw()
                        self.currentViewSelected = self.imgType

                        return

                if(self.dms):
                    if(self.imgType == "Phase - DM"):
                        self.setupDisp("pg")
                        dm_type = self.config.p_dms[
                            self.numberSelected].type_dm
                        alt = self.config.p_dms[self.numberSelected].alt
                        data = self.dms.get_dm(dm_type, alt)

                        if(self.imgType != self.currentViewSelected):
                            self.p1.setRange(
                                xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                        self.currentViewSelected = self.imgType
                        self.SRcircleDM[self.numberSelected].show()
                if(self.tar):
                    if(self.imgType == "Phase - Target"):
                        self.setupDisp("pg")
                        data = self.tar.get_phase(self.numberSelected)
                        if(self.imgType != self.currentViewSelected):
                            self.p1.setRange(
                                xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                        self.currentViewSelected = self.imgType
                        self.SRcircleTarget[self.numberSelected].show()

                    if(self.imgType == "PSF SE"):
                        self.setupDisp("pg")
                        data = self.tar.get_image(
                            self.numberSelected, "se")
                        if(self.ui.wao_PSFlogscale.isChecked()):
                            if np.any(data <= 0):
                                print((
                                    "\nzero founds, log display disabled\n", RuntimeWarning))
                                self.ui.wao_PSFlogscale.setCheckState(False)
                            else:
                                data = np.log10(data)

                        if (not self.SRCrossX):
                            Delta = 5
                            self.SRCrossX = pg.PlotCurveItem(np.array([data.shape[0] / 2 + 0.5 - Delta,
                                                                       data.shape[0] / 2 + 0.5 + Delta]),
                                                             np.array([data.shape[1] / 2 + 0.5,
                                                                       data.shape[1] / 2 + 0.5]),
                                                             pen='r')
                            self.SRCrossY = pg.PlotCurveItem(np.array([data.shape[0] / 2 + 0.5,
                                                                       data.shape[0] / 2 + 0.5]),
                                                             np.array([data.shape[1] / 2 + 0.5 - Delta,
                                                                       data.shape[1] / 2 + 0.5 + Delta]),
                                                             pen='r')
                            # Put image in plot area
                            self.p1.addItem(self.SRCrossX)
                            # Put image in plot area
                            self.p1.addItem(self.SRCrossY)

                        if(self.imgType != self.currentViewSelected):
                            zoom = 50
                            self.SRCrossX.show()
                            self.SRCrossY.show()
                            self.p1.setRange(xRange=(data.shape[0] / 2 + 0.5 - zoom,
                                                     data.shape[0] / 2 + 0.5 + zoom),
                                             yRange=(data.shape[1] / 2 + 0.5 - zoom,
                                                     data.shape[1] / 2 + 0.5 + zoom),)
                        self.currentViewSelected = self.imgType

                    if(self.imgType == "PSF LE"):
                        self.setupDisp("pg")
                        data = self.tar.get_image(
                            self.numberSelected, "le")
                        if(self.ui.wao_PSFlogscale.isChecked()):
                            data = np.log10(data)
                        if (not self.SRCrossX):
                            Delta = 5
                            self.SRCrossX = pg.PlotCurveItem(np.array([data.shape[0] / 2 + 0.5 - Delta,
                                                                       data.shape[0] / 2 + 0.5 + Delta]),
                                                             np.array([data.shape[1] / 2 + 0.5,
                                                                       data.shape[1] / 2 + 0.5]),
                                                             pen='r')
                            self.SRCrossY = pg.PlotCurveItem(np.array([data.shape[0] / 2 + 0.5,
                                                                       data.shape[0] / 2 + 0.5]),
                                                             np.array([data.shape[1] / 2 + 0.5 - Delta,
                                                                       data.shape[1] / 2 + 0.5 + Delta]),
                                                             pen='r')

                            # Put image in plot area
                            self.p1.addItem(self.SRCrossX)
                            # Put image in plot area
                            self.p1.addItem(self.SRCrossY)
                        if(self.imgType != self.currentViewSelected):
                            zoom = 50
                            self.p1.setRange(xRange=(data.shape[0] / 2 + 0.5 - zoom,
                                                     data.shape[0] / 2 + 0.5 + zoom),
                                             yRange=(data.shape[1] / 2 + 0.5 - zoom,
                                                     data.shape[1] / 2 + 0.5 + zoom))
                            self.SRCrossX.show()
                            self.SRCrossY.show()

                        self.currentViewSelected = self.imgType

                if(data is not None):
                    autoscale = self.ui.wao_autoscale.isChecked()
                    if(autoscale):
                        # inits levels
                        self.hist.setLevels(data.min(), data.max())
                    self.img.setImage(data, autoLevels=autoscale)
                    # self.p1.autoRange()
            finally:
                self.loopLock.release()

            refreshDisplayTime = 1000. / self.ui.wao_frameRate.value()
            if(self.ui.wao_Display.isChecked()):
                # Update GUI plots
                QTimer.singleShot(refreshDisplayTime, self.updateDisplay)

    def mainLoop(self):
        if not self.loopLock.acquire(False):
            # print " Loop locked"
            return
        else:
            try:
                start = time.time()
                self.atm.move_atmos()
                if(self.config.p_controllers[0].type_control == "geo"):
                    for t in range(self.config.p_target.ntargets):
                        if wao.see_atmos:
                            self.tar.atmos_trace(t, self.atm, self.tel)
                        else:
                            self.tar.reset_phase(t)
                        self.rtc.docontrol_geo(0, self.dms, self.tar, 0)
                        self.rtc.applycontrol(0, self.dms)
                        self.tar.dmtrace(0, self.dms)
                else:
                    for t in range(self.config.p_target.ntargets):
                        if wao.see_atmos:
                            self.tar.atmos_trace(t, self.atm, self.tel)
                        else:
                            self.tar.reset_phase(t)
                        self.tar.dmtrace(t, self.dms)
                    for w in range(len(self.config.p_wfss)):
                        if wao.see_atmos:
                            self.wfs.sensors_trace(
                                w, "atmos", self.tel, self.atm, self.dms)
                        else:
                            self.wfs.reset_phase(w)
                        if not self.config.p_wfss[w].openloop:
                            self.wfs.sensors_trace(
                                w, "dm", self.tel, self.atm, self.dms)
                        self.wfs.sensors_compimg(w)

                    self.rtc.docentroids(0)
                    self.rtc.docontrol(0)
                    self.rtc.doclipping(0, -1e5, 1e5)
                    self.rtc.applycontrol(0, self.dms)

                if(wao.brama_rtc_flag):
                    self.rtc.publish()
                if(wao.brama_tar_flag):
                    self.tar.publish()

                refreshDisplayTime = 1. / self.ui.wao_frameRate.value()
                if(time.time() - self.refreshTime > refreshDisplayTime):
                    signal_le = ""
                    signal_se = ""
                    for t in range(self.config.p_target.ntargets):
                        SR = self.tar.get_strehl(t)
                        if(t==self.numberSelected): # Plot on the wfs selected
                            self.updateSRDisplay(SR[1], SR[0], self.iter)

                        signal_se += "%1.2f   " % SR[0]
                        signal_le += "%1.2f   " % SR[1]

                    loopTime = time.time() - start
                    currentFreq = 1 / loopTime
                    refreshFreq = 1 / (time.time() - self.refreshTime)

                    self.ui.wao_strehlSE.setText(signal_se)
                    self.ui.wao_strehlLE.setText(signal_le)
                    self.ui.wao_currentFreq.setValue(1 / loopTime)
                    if(self.dispStatsInTerminal):
                        self.printInPlace("iter #%d SR: (L.E, S.E.)= %s, %srunning at %4.1fHz (real %4.1fHz)" % (
                            self.iter, signal_le, signal_se, refreshFreq, currentFreq))
                    self.refreshTime = time.time()
                self.iter += 1
            finally:
                self.loopLock.release()

    def printInPlace(self, text):
        # This seems to trigger the GUI and keep it responsive
        print("\r" + text, end=' ')
        sys.stdout.flush()
        # sys.stdout.write(text)

    def run(self):
        # print "Loop started"
        self.mainLoop()
        if not self.stop:
            QTimer.singleShot(0, self.run)

        # print "Loop stopped"


class WorkerThread( QThread ):
    def __init__( self, parentThread, parentLoop):
        QThread.__init__( self, parentThread)
        self.loop = parentLoop
    def run( self ):
        self.running = True
        self.loop()
        success = True
        self.emit( SIGNAL( "jobFinished( PyQt_PyObject )" ), success )
    def stop( self ):
        self.running = False
        pass
    def cleanUp( self):
        pass


try:
    class widgetAOWindowPyro(widgetAOWindow):
        def __init__(self):
            widgetAOWindow.__init__(self)
            #Pyro.core.ObjBase.__init__(self)

        def InitConfig(self):
            widgetAOWindow.InitConfig(self)
            try:
                ps = PyroServer(wao)
                ps.start()
            except:
                print("Warning: Error while starting Pyro server")

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
    app = QtGui.QApplication(sys.argv)
    # app.setStyleSheet(qdarkstyle.load_stylesheet(pyside=False))
    app.setStyle('cleanlooks')
    if USE_PYRO:
        wao = widgetAOWindowPyro()
    else:
        wao = widgetAOWindow()

    wao.show()
    """
    locator = Pyro.naming.NameServerLocator()
    ns = locator.getNS()
    Pyro.core.initServer()
    daemon = Pyro.core.Daemon()
    daemon.useNameServer(ns)  # use current ns
    daemon.connect(test, ":waoCanapass")
    daemon.requestLoop()
    """

    """
    daemon = Pyro4.Daemon()                # make a Pyro daemon
    ns = Pyro4.locateNS()                  # find the name server
    uri = daemon.register(wao)   # register the greeting maker as a Pyro object
    ns.register("example.greeting", uri)   # register the object with a name in the name server

        print("Ready.")
        daemon.requestLoop()
    except:
        from warnings import warn
        warn("pyro4 not found", RuntimeWarning)
    """
    # app.connect(wao.ui._quit,QtCore.SIGNAL("clicked()"),app,QtCore.SLOT("quit()"))
