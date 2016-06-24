"""
widget_ao.

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
import glob
import tools
import hdf5_utils as h5u
import threading
from PyQt4.uic import loadUiType
from PyQt4 import QtGui
from functools import partial
import subprocess

sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/data/par/")

WindowTemplate, TemplateBaseClass = loadUiType(os.environ["SHESHA_ROOT"] +
                                               "/widgets/widget_ao.ui")

plt.ion()

"""
low levels debugs:
gdb --args python -i widget_ao.py

"""


class widgetAOWindow(TemplateBaseClass):

    def __init__(self):
        TemplateBaseClass.__init__(self)

        self.ui = WindowTemplate()
        self.ui.setupUi(self)

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
        self.displayLock = threading.Lock()
        self.loopLock = threading.Lock()
        self.iter = 0
        self.loaded = False
        self.stop = False
        self.nbiters = 1000
        self.startTime = 0
        self.loop = None
        self.assistant = None
        self.selector_init = None


        #############################################################
        #                 PYQTGRAPH WINDOW INIT                     #
        #############################################################

        self.img = pg.ImageItem(border='w')  # create image area
        # self.p1 = self.ui.wao_pgwindow.addPlot()  # create pyqtgraph plot
        # area
        self.p1 = self.ui.wao_pgwindow.addViewBox()
        self.p1.setAspectLocked(True)
        self.p1.addItem(self.img)  # Put image in plot area

        self.hist = pg.HistogramLUTItem()  # Create an histogram
        self.hist.setImageItem(self.img)  # Compute histogram from img
        self.ui.wao_pgwindow.addItem(self.hist)
        self.hist.autoHistogramRange()  # init levels
        self.hist.setMaximumWidth(100)

        #############################################################
        #                 CONNECTED BUTTONS                         #
        #############################################################
        # Default path for config files
        self.defaultParPath = os.environ[
            "SHESHA_ROOT"] + "/data/par/par4bench/"
        self.ui.wao_loadConfig.clicked.connect(self.loadConfig)
        self.loadDefaultConfig()
        self.ui.wao_init.clicked.connect(self.InitConfig)
        self.ui.wao_run.setCheckable(True)
        self.ui.wao_run.clicked[bool].connect(self.aoLoopClicked)
        self.ui.wao_openLoop.setCheckable(True)
        self.ui.wao_openLoop.clicked[bool].connect(self.aoLoopOpen)
        self.ui.wao_next.clicked.connect(self.mainLoop)
        self.imgType = str(self.ui.wao_selectScreen.currentText())
        self.ui.wao_configFromFile.clicked.connect(self.addConfigFromFile)
        self.ui.wao_unzoom.clicked.connect(self.p1.autoRange)
        self.ui.wao_selectScreen.currentIndexChanged.connect(
            partial(self.updateNumberSelector, textType=None))
        self.ui.wao_selectNumber.currentIndexChanged.connect(
            self.setNumberSelection)
        self.ui.wao_wfs_plotSelector.currentIndexChanged.connect(
            self.updatePlotWfs)
        self.ui.wao_selectAtmosLayer.currentIndexChanged.connect(
            self.setLayerSelection)
        self.ui.wao_selectWfs.currentIndexChanged.connect(self.setWfsSelection)
        self.ui.wao_selectDM.currentIndexChanged.connect(self.setDmSelection)
        self.ui.wao_selectCentro.currentIndexChanged.connect(
            self.setCentroSelection)
        self.ui.wao_selectTarget.currentIndexChanged.connect(
            self.setTargetSelection)
        self.ui.wao_setAtmos.clicked.connect(self.setAtmosParams)
        self.ui.wao_setWfs.clicked.connect(self.setWfsParams)
        self.ui.wao_setDM.clicked.connect(self.setDmParams)
        self.ui.wao_setControl.clicked.connect(self.setRtcParams)
        self.ui.wao_setCentro.clicked.connect(self.setRtcParams)
        self.ui.wao_setTelescope.clicked.connect(self.setTelescopeParams)
        self.ui.wao_resetDM.clicked.connect(self.resetDM)
        self.ui.wao_update_gain.clicked.connect(self.updateGain)
        self.ui.wao_update_pyr_ampl.clicked.connect(self.updatePyrAmpl)
        self.ui.wao_selectRtcMatrix.currentIndexChanged.connect(
            self.displayRtcMatrix)
        self.ui.wao_rtcWindowMPL.hide()
        self.ui.wao_commandBtt.clicked.connect(self.BttCommand)
        self.ui.wao_Display.clicked.connect(self.updateFrameRate)
        self.ui.wao_frameRate.valueChanged.connect(self.updateFrameRate)
        self.RTDisplay = self.ui.wao_Display.isChecked()
        self.RTDFreq = self.ui.wao_frameRate.value()
        self.ui.wao_PSFlogscale.clicked.connect(self.updateDisplay)
        self.ui.wao_resetSR.clicked.connect(self.resetSR)
        self.ui.wao_actionHelp_Contents.triggered.connect(
            self.on_help_triggered)

        self.ui.wao_dmUnitPerVolt.valueChanged.connect(self.updateDMrangeGUI)
        self.ui.wao_dmpush4iMat.valueChanged.connect(self.updateDMrangeGUI)
        self.ui.wao_pyr_ampl.valueChanged.connect(self.updateAmpliCompGUI)
        self.ui.wao_dmActuPushArcSecNumWFS.currentIndexChanged.connect(self.updateDMrange)
        self.ui.wao_nbiters.valueChanged.connect(self.updateNBiters)
        self.SRcircleAtmos = {}
        self.SRcircleWFS = {}
        self.SRcircleDM = {}
        self.SRcircleTarget = {}
        self.ui.splitter.setSizes([2000, 10])

    def updateNBiters(self):
        self.nbiters = int(self.ui.wao_nbiters.value())

    def on_help_triggered(self, i=None):
        if i is None:
            return

        if not self.assistant or \
           not self.assistant.poll() is None:

            helpcoll = os.environ["COMPASS_ROOT"] + "/doc/COMPASS.qhc"
            cmd = "assistant -enableRemoteControl -collectionFile %s" % helpcoll
            self.assistant = subprocess.Popen(cmd,
                                              shell=True,
                                              stdin=subprocess.PIPE)
#         self.assistant.stdin.write("SetSource qthelp://org.sphinx.compassshesha.r763/doc/index.html\n")

    def resetSR(self):
        tarnum = self.ui.wao_resetSR_tarNum.value()
        print "reset SR on target %d" % tarnum
        self.tar.reset_strehl(tarnum)

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
        #                       METHODS                             #
        #############################################################
    def updateGain(self):
        if(self.rtc):
            self.rtc.set_gain(0,float(self.ui.wao_controlGain.value()))
            print "Loop gain updated on GPU"

    def updateAmpliCompGUI(self):
        diffract = self.config.p_wfss[0].Lambda*1e-6 / self.config.p_tel.diam*206265.
        self.ui.wao_pyr_ampl_arcsec.setValue(self.ui.wao_pyr_ampl.value()*diffract)

    def updateAmpliComp(self):
        diffract = self.config.p_wfss[0].Lambda*1e-6 / self.config.p_tel.diam*206265.
        self.ui.wao_pyr_ampl_arcsec.setValue(self.config.p_wfss[0].pyr_ampl*diffract)

    def updatePyrAmpl(self):
        if(self.rtc):
            self.rtc.set_pyr_ampl(0,self.ui.wao_pyr_ampl.value(),self.config.p_wfss,self.config.p_tel)
            print "Pyramid modulation updated on GPU"
            self.updatePlotWfs()

    def updateFrameRate(self):
        if(self.ui.wao_Display.isChecked()):
            self.RTDisplay = True
        else:
            self.RTDisplay = False
        self.RTDFreq = self.ui.wao_frameRate.value()

    def updateDMrangeGUI(self):
        push4imat=self.ui.wao_dmpush4iMat.value()
        unitpervolt=self.ui.wao_dmUnitPerVolt.value()
        self.updateDMrange(push4imat=push4imat, unitpervolt=unitpervolt)

    def updateDMrange(self, push4imat=None, unitpervolt=None):
        numdm = str(self.ui.wao_selectDM.currentText())
        numwfs = str(self.ui.wao_dmActuPushArcSecNumWFS.currentText())
        if((numdm is not "") and (numwfs is not "") and (push4imat!=0) and (unitpervolt!=0)):
            arcsecDMrange = self.computeDMrange(int(numdm) ,int(numwfs), push4imat=push4imat, unitpervolt=unitpervolt)
            self.ui.wao_dmActPushArcsec.setValue(arcsecDMrange)

    def updateTelescopePanel(self):
        self.ui.wao_zenithAngle.setValue(self.config.p_geom.zenithangle)
        self.ui.wao_diamTel.setValue(self.config.p_tel.diam)
        self.ui.wao_cobs.setValue(self.config.p_tel.cobs)

    def updateDmPanel(self):
        ndm = self.ui.wao_selectDM.currentIndex()
        if(ndm < 0):
            ndm = 0
        self.ui.wao_dmActuPushArcSecNumWFS.clear()
        self.ui.wao_dmActuPushArcSecNumWFS.addItems([str(i) for i in range(len(self.config.p_wfss))])
        self.ui.wao_numberofDMs.setText(str(len(self.config.p_dms)))
        self.ui.wao_dmTypeSelector.setCurrentIndex(
            self.ui.wao_dmTypeSelector.findText(self.config.p_dms[ndm].type_dm))
        self.ui.wao_dmAlt.setValue(self.config.p_dms[ndm].alt)
        self.ui.wao_dmNactu.setValue(self.config.p_dms[ndm].nact)
        self.ui.wao_dmUnitPerVolt.setValue(self.config.p_dms[ndm].unitpervolt)
        self.ui.wao_dmpush4iMat.setValue(self.config.p_dms[ndm].push4imat)
        self.ui.wao_dmCoupling.setValue(self.config.p_dms[ndm].coupling)
        self.ui.wao_dmThresh.setValue(self.config.p_dms[ndm].thresh)
        self.updateDMrange()

    def updateWfsPanel(self):
        nwfs = self.ui.wao_selectWfs.currentIndex()
        if(nwfs < 0):
            nwfs = 0
        self.ui.wao_numberofWfs.setText(str(len(self.config.p_wfss)))
        self.ui.wao_wfsType.setText(str(self.config.p_wfss[nwfs].type_wfs))
        self.ui.wao_wfsNxsub.setValue(self.config.p_wfss[nwfs].nxsub)
        self.ui.wao_wfsNpix.setValue(self.config.p_wfss[nwfs].npix)
        self.ui.wao_wfsPixSize.setValue(self.config.p_wfss[nwfs].pixsize)
        self.ui.wao_wfsXpos.setValue(self.config.p_wfss[nwfs].xpos)
        self.ui.wao_wfsYpos.setValue(self.config.p_wfss[nwfs].ypos)
        self.ui.wao_wfsFracsub.setValue(self.config.p_wfss[nwfs].fracsub)
        self.ui.wao_wfsLambda.setValue(self.config.p_wfss[nwfs].Lambda)
        self.ui.wao_wfsMagnitude.setValue(self.config.p_wfss[nwfs].gsmag)
        self.ui.wao_wfsZp.setValue(np.log10(self.config.p_wfss[nwfs].zerop))
        self.ui.wao_wfsThrough.setValue(self.config.p_wfss[nwfs].optthroughput)
        self.ui.wao_wfsNoise.setValue(self.config.p_wfss[nwfs].noise)
        self.ui.wao_pyr_ampl.setValue(self.config.p_wfss[nwfs].pyr_ampl)

        # LGS panel
        if(self.config.p_wfss[nwfs].gsalt > 0):
            self.ui.wao_wfsIsLGS.setChecked(True)
            self.ui.wao_wfsGsAlt.setValue(self.config.p_wfss[nwfs].gsalt)
            self.ui.wao_wfsLLTx.setValue(self.config.p_wfss[nwfs].lltx)
            self.ui.wao_wfsLLTy.setValue(self.config.p_wfss[nwfs].llty)
            self.ui.wao_wfsLGSpower.setValue(
                self.config.p_wfss[nwfs].laserpower)
            self.ui.wao_wfsReturnPerWatt.setValue(
                self.config.p_wfss[nwfs].lgsreturnperwatt)
            self.ui.wao_wfsBeamSize.setValue(self.config.p_wfss[nwfs].beamsize)
            self.ui.wao_selectLGSProfile.setCurrentIndex(
                self.ui.wao_selectLGSProfile.findText(self.config.p_wfss[nwfs].proftype))

        else:
            self.ui.wao_wfsIsLGS.setChecked(False)

        if(self.config.p_wfss[nwfs].type_wfs == "pyrhr" or self.config.p_wfss[nwfs].type_wfs == "pyr"):
            self.ui.wao_wfs_plotSelector.setCurrentIndex(3)
        self.updatePlotWfs()


    def updateAtmosPanel(self):
        nscreen = self.ui.wao_selectAtmosLayer.currentIndex()
        if(nscreen < 0):
            nscreen = 0
        self.ui.wao_r0.setValue(self.config.p_atmos.r0)
        self.ui.wao_atmosNlayers.setValue(self.config.p_atmos.nscreens)
        self.ui.wao_atmosAlt.setValue(self.config.p_atmos.alt[nscreen])
        self.ui.wao_atmosFrac.setValue(self.config.p_atmos.frac[nscreen])
        self.ui.wao_atmosL0.setValue(self.config.p_atmos.L0[nscreen])
        self.ui.wao_windSpeed.setValue(self.config.p_atmos.windspeed[nscreen])
        self.ui.wao_windDirection.setValue(
            self.config.p_atmos.winddir[nscreen])
        if(self.config.p_atmos.dim_screens is not None):
            self.ui.wao_atmosDimScreen.setText(
                str(self.config.p_atmos.dim_screens[nscreen]))
        self.ui.wao_atmosWindow.canvas.axes.cla()
        width = (self.config.p_atmos.alt.max() / 20. + 0.1) / 1000.
        self.ui.wao_atmosWindow.canvas.axes.barh(
            self.config.p_atmos.alt / 1000. - width / 2., self.config.p_atmos.frac, width, color="blue")
        self.ui.wao_atmosWindow.canvas.draw()

    def updateRtcPanel(self):
        # Centroider panel
        ncentro = self.ui.wao_selectCentro.currentIndex()
        if(ncentro < 0):
            ncentro = 0
        self.ui.wao_centroTypeSelector.setCurrentIndex(
            self.ui.wao_centroTypeSelector.findText(self.config.p_centroiders[ncentro].type_centro))
        self.ui.wao_centroThresh.setValue(
            self.config.p_centroiders[ncentro].thresh)
        self.ui.wao_centroNbrightest.setValue(
            self.config.p_centroiders[ncentro].nmax)
        self.ui.wao_centroThresh.setValue(
            self.config.p_centroiders[ncentro].thresh)
        if(self.config.p_centroiders[ncentro].type_fct):
            self.ui.wao_centroFunctionSelector.setCurrentIndex(
                self.ui.wao_centroFunctionSelector.findText(self.config.p_centroiders[ncentro].type_fct))
        self.ui.wao_centroWidth.setValue(
            self.config.p_centroiders[ncentro].width)

        # Controller panel
        type_contro = self.config.p_controllers[0].type_control
        if(type_contro == "ls" and self.config.p_controllers[0].modopti == 0):
            self.ui.wao_controlTypeSelector.setCurrentIndex(0)
        elif(type_contro == "mv"):
            self.ui.wao_controlTypeSelector.setCurrentIndex(1)
        elif(type_contro == "geo"):
            self.ui.wao_controlTypeSelector.setCurrentIndex(2)
        elif(type_contro == "ls" and self.config.p_controllers[0].modopti):
            self.ui.wao_controlTypeSelector.setCurrentIndex(3)
        elif(type_contro == "cured"):
            self.ui.wao_controlTypeSelector.setCurrentIndex(4)
        else:
            print "pffff...."

        self.ui.wao_controlCond.setValue(self.config.p_controllers[0].maxcond)
        self.ui.wao_controlDelay.setValue(self.config.p_controllers[0].delay)
        self.ui.wao_controlGain.setValue(self.config.p_controllers[0].gain)
        self.ui.wao_controlTTcond.setValue(
            self.config.p_controllers[0].maxcond)  # TODO : TTcond

    def updateTargetPanel(self):
        ntarget = self.ui.wao_selectTarget.currentIndex()
        if(ntarget < 0):
            ntarget = 0
        self.ui.wao_numberofTargets.setText(str(self.config.p_target.ntargets))
        self.ui.wao_targetMag.setValue(self.config.p_target.mag[ntarget])
        self.ui.wao_targetXpos.setValue(self.config.p_target.xpos[ntarget])
        self.ui.wao_targetYpos.setValue(self.config.p_target.ypos[ntarget])
        self.ui.wao_targetLambda.setValue(self.config.p_target.Lambda[ntarget])

        self.ui.wao_targetWindow.canvas.axes.cla()
        xmax = np.max(np.abs(self.config.p_target.xpos))
        ymax = np.max(np.abs(self.config.p_target.ypos))
        if(self.config.p_wfss):
            self.ui.wao_targetWindow.canvas.axes.plot([w.xpos for w in self.config.p_wfss], [
                                                      w.ypos for w in self.config.p_wfss], 'o', color="green")
            xmax = np.max(
                [xmax, np.max(np.abs([w.xpos for w in self.config.p_wfss]))])
            ymax = np.max(
                [ymax, np.max(np.abs([w.ypos for w in self.config.p_wfss]))])
        self.ui.wao_targetWindow.canvas.axes.plot(
            self.config.p_target.xpos, self.config.p_target.ypos, '*', color="red")
        self.ui.wao_targetWindow.canvas.axes.set_xlim(-xmax - 10, xmax + 10)
        self.ui.wao_targetWindow.canvas.axes.set_ylim(-ymax - 10, ymax + 10)
        self.ui.wao_targetWindow.canvas.axes.grid()
        self.ui.wao_targetWindow.canvas.draw()

    def updatePanels(self):
        self.updateTelescopePanel()
        self.updateLayerSelection()
        self.updateAtmosPanel()
        self.updateWfsSelection()
        self.updateWfsPanel()
        self.updateDmSelection()
        self.updateDmPanel()
        self.updateCentroSelection()
        self.updateRtcPanel()
        self.updateTargetSelection()
        self.updateTargetPanel()

    def setTelescopeParams(self):
        self.config.p_tel.set_diam(self.ui.wao_diamTel.value())
        self.config.p_tel.set_cobs(self.ui.wao_cobs.value())
        self.config.p_geom.set_zenithangle(self.ui.wao_zenithAngle.value())
        print "New telescope parameters set"

    def setAtmosParams(self):
        nscreen = self.ui.wao_selectAtmosLayer.currentIndex()
        if(nscreen < 0):
            nscreen = 0
        self.config.p_atmos.alt[nscreen] = self.ui.wao_atmosAlt.value()
        self.config.p_atmos.frac[nscreen] = self.ui.wao_atmosFrac.value()
        self.config.p_atmos.L0[nscreen] = self.ui.wao_atmosL0.value()
        self.config.p_atmos.windspeed[nscreen] = self.ui.wao_windSpeed.value()
        self.config.p_atmos.winddir[
            nscreen] = self.ui.wao_windDirection.value()
        print "New atmos parameters set"

    def setRtcParams(self):
        # Centroider params
        ncentro = self.ui.wao_selectCentro.currentIndex()
        if(ncentro < 0):
            ncentro = 0
        self.config.p_centroiders[ncentro].set_type(
            str(self.ui.wao_centroTypeSelector.currentText()))
        self.config.p_centroiders[ncentro].set_thresh(
            self.ui.wao_centroThresh.value())
        self.config.p_centroiders[ncentro].set_nmax(
            self.ui.wao_centroNbrightest.value())
        self.config.p_centroiders[ncentro].set_thresh(
            self.ui.wao_centroThresh.value())
        self.config.p_centroiders[ncentro].set_type_fct(
            str(self.ui.wao_centroFunctionSelector.currentText()))
        self.config.p_centroiders[ncentro].set_width(
            self.ui.wao_centroWidth.value())

        # Controller panel
        type_contro = str(self.ui.wao_controlTypeSelector.currentText())
        if(type_contro == "LS"):
            self.config.p_controllers[0].set_type("ls")
        elif(type_contro == "MV"):
            self.config.p_controllers[0].set_type("mv")
        elif(type_contro == "PROJ"):
            self.config.p_controllers[0].set_type("geo")
        elif(type_contro == "OptiMods"):
            self.config.p_controllers[0].set_type("ls")
            self.config.p_controllers[0].set_modopti(1)

        self.config.p_controllers[0].set_maxcond(
            self.ui.wao_controlCond.value())
        self.config.p_controllers[0].set_delay(
            self.ui.wao_controlDelay.value())
        self.config.p_controllers[0].set_gain(self.ui.wao_controlGain.value())
        # self.config.p_controllers[0].set_TTcond(self.ui.wao_controlTTcond.value())
        # # TODO : TTcond
        print "New rtc parameters set"

    def setWfsParams(self):
        nwfs = self.ui.wao_selectWfs.currentIndex()
        if(nwfs < 0):
            nwfs = 0
        self.config.p_wfss[nwfs].set_nxsub(self.ui.wao_wfsNxsub.value())
        self.config.p_wfss[nwfs].set_npix(self.ui.wao_wfsNpix.value())
        self.config.p_wfss[nwfs].set_pixsize(self.ui.wao_wfsPixSize.value())
        self.config.p_wfss[nwfs].set_xpos(self.ui.wao_wfsXpos.value())
        self.config.p_wfss[nwfs].set_ypos(self.ui.wao_wfsYpos.value())
        self.config.p_wfss[nwfs].set_fracsub(self.ui.wao_wfsFracsub.value())
        self.config.p_wfss[nwfs].set_Lambda(self.ui.wao_wfsLambda.value())
        self.config.p_wfss[nwfs].set_gsmag(self.ui.wao_wfsMagnitude.value())
        # TODO: find a way to correctly set zerop (limited by the maximum value
        # allowed by the double spin box)
        self.config.p_wfss[nwfs].set_zerop(10 ** (self.ui.wao_wfsZp.value()))
        self.config.p_wfss[nwfs].set_optthroughput(
            self.ui.wao_wfsThrough.value())
        self.config.p_wfss[nwfs].set_noise(self.ui.wao_wfsNoise.value())

        # LGS params
        if(self.ui.wao_wfsIsLGS.isChecked()):
            self.config.p_wfss[nwfs].set_gsalt(self.ui.wao_wfsGsAlt.value())
            self.config.p_wfss[nwfs].set_lltx(self.ui.wao_wfsLLTx.value())
            self.config.p_wfss[nwfs].set_llty(self.ui.wao_wfsLLTy.value())
            self.config.p_wfss[nwfs].set_laserpower(
                self.ui.wao_wfsLGSpower.value())
            self.config.p_wfss[nwfs].set_lgsreturnperwatt(
                self.ui.wao_wfsReturnPerWatt.value())
            self.config.p_wfss[nwfs].set_beamsize(
                self.ui.wao_wfsBeamSize.value())
            self.config.p_wfss[nwfs].set_proftype(
                str(self.ui.wao_selectLGSProfile.currentText()))
        print "New wfs parameters set"

    def setDmParams(self):
        ndm = self.ui.wao_selectDM.currentIndex()
        if(ndm < 0):
            ndm = 0
        self.config.p_dms[ndm].set_type(
            str(self.ui.wao_dmTypeSelector.currentText()))
        self.config.p_dms[ndm].set_alt(self.ui.wao_dmAlt.value())
        self.config.p_dms[ndm].set_nact(self.ui.wao_dmNactu.value())
        self.config.p_dms[ndm].set_unitpervolt(
            self.ui.wao_dmUnitPerVolt.value())
        self.config.p_dms[ndm].set_coupling(self.ui.wao_dmCoupling.value())
        self.config.p_dms[ndm].set_thresh(self.ui.wao_dmThresh.value())
        print "New DM parameters set"

    def updateLayerSelection(self):
        self.ui.wao_selectAtmosLayer.clear()
        self.ui.wao_selectAtmosLayer.addItems(
            [str(i) for i in range(self.config.p_atmos.nscreens)])

    def updateTargetSelection(self):
        self.ui.wao_selectTarget.clear()
        self.ui.wao_selectTarget.addItems(
            [str(i) for i in range(self.config.p_target.ntargets)])

    def updateWfsSelection(self):
        self.ui.wao_selectWfs.clear()
        self.ui.wao_selectWfs.addItems(
            [str(i) for i in range(len(self.config.p_wfss))])

    def updateDmSelection(self):
        self.ui.wao_selectDM.clear()
        self.ui.wao_selectDM.addItems(
            [str(i) for i in range(len(self.config.p_dms))])

    def updateCentroSelection(self):
        self.ui.wao_selectCentro.clear()
        self.ui.wao_selectCentro.addItems(
            [str(i) for i in range(len(self.config.p_centroiders))])

    def setCentroSelection(self):
        self.updateRtcPanel()

    def setLayerSelection(self):
        self.updateAtmosPanel()

    def setTargetSelection(self):
        self.updateTargetPanel()

    def setWfsSelection(self):
        self.updateWfsPanel()

    def setDmSelection(self):
        self.updateDmPanel()

    def addConfigFromFile(self):
        filepath = str(QtGui.QFileDialog(directory=self.defaultParPath).getOpenFileName(
            self, "Select parameter file", "", "parameters file (*.py);;hdf5 file (*.h5);;all files (*)"))
        self.configpath = filepath
        filename = filepath.split('/')[-1]
        if(filepath.split('.')[-1] == "py"):
            self.ui.wao_selectConfig.addItem(filename, 0)
            pathfile = filepath.split(filename)[0]
            if (pathfile not in sys.path):
                sys.path.insert(0, pathfile)
            exec("import %s as config" % filename.split(".py")[0])
            sys.path.remove(pathfile)
        elif(filepath.split('.')[-1] == "h5"):
            sys.path.insert(0, self.defaultParPath)
            import scao_16x16_8pix as config
            sys.path.remove(self.defaultParPath)
            h5u.configFromH5(filepath, config)
        else:
            print "Parameter file extension must be .py or .h5"
            return
        self.config = config
        self.ui.wao_selectConfig.clear()
        self.ui.wao_selectConfig.addItem(filename)
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
        self.updatePanels()
        self.ui.wao_init.setDisabled(False)

    def aoLoopClicked(self, pressed):
        if(pressed):
            self.loop = threading.Thread(target=self.run)
            self.loop.start()
        else:
            self.stop = True
            self.loop.join()
            self.loop = None

    def aoLoopOpen(self, pressed):
        if(pressed):
            self.rtc.set_openloop(0, 1)
        else:
            self.rtc.set_openloop(0, 0)

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
        if(textType == "Phase - WFS" or textType == "Spots - WFS" or textType == "Centroids - WFS" or textType == "Slopes - WFS" or textType == "Pyrimg - HR" or textType == "Pyrimg - LR"):
            n = len(self.config.p_wfss)
        if(textType == "Phase - Target" or textType == "PSF LE" or textType == "PSF SE"):
            n = self.config.p_target.ntargets
        if(textType == "Phase - DM"):
            n = len(self.config.p_dms)
        self.ui.wao_selectNumber.addItems([str(i) for i in range(n)])
        self.updateDisplay()

    def loadConfig(self):
        configfile = str(self.ui.wao_selectConfig.currentText())
        self.configpath = self.defaultParPath + configfile
        sys.path.insert(0, self.defaultParPath)
        exec("import %s as config" % configfile.split(".py")[0])
        self.config = config
        self.loaded = False
        self.ui.wao_selectScreen.clear()
        if(self.config.p_wfss[0].type_wfs == "pyrhr"):
            self.selector_init = ["Phase - Atmos", "Phase - WFS", "Pyrimg - LR",
                                  "Pyrimg - HR", "Centroids - WFS", "Slopes - WFS",
                                  "Phase - Target", "Phase - DM",
                                  "PSF LE", "PSF SE"]
            self.ui.ui_modradiusPanel.show()
            self.ui.ui_modradiusPanelarcesec.show()
            self.ui.wao_pyr_ampl.show()
            self.ui.wao_pyr_ampl_arcsec.show()
            self.ui.wao_update_pyr_ampl.show()
        else:
            self.selector_init = ["Phase - Atmos", "Phase - WFS", "Spots - WFS",
                                  "Centroids - WFS", "Slopes - WFS",
                                  "Phase - Target", "Phase - DM",
                                  "PSF LE", "PSF SE"]
            self.ui.ui_modradiusPanel.hide()
            self.ui.ui_modradiusPanelarcesec.hide()
            self.ui.wao_pyr_ampl.hide()
            self.ui.wao_pyr_ampl_arcsec.hide()
            self.ui.wao_update_pyr_ampl.hide()
        self.ui.wao_selectScreen.addItems(self.selector_init)
        self.ui.wao_selectScreen.setCurrentIndex(0)
        self.updateNumberSelector(textType=self.imgType)
        self.updatePanels()
        self.ui.wao_init.setDisabled(False)

    def loadDefaultConfig(self):
        parlist = sorted(glob.glob(self.defaultParPath + "*.py"))
        self.ui.wao_selectConfig.clear()
        self.ui.wao_selectConfig.addItems(
            [parlist[i].split('/')[-1] for i in range(len(parlist))])

    def InitConfig(self):
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
        self.currentViewSelected = None
        self.SRCrossX = None
        self.SRCrossY = None

        # set simulation name
        if(hasattr(self.config, "simul_name")):
            if(self.config.simul_name is None):
                simul_name = ""
            else:
                simul_name = self.config.simul_name
        else:
            simul_name = ""
        matricesToLoad = {}
        if(simul_name == "" or not self.ui.wao_useDatabase.isChecked()):
            clean = 1
        else:
            clean = 0
            param_dict = h5u.params_dictionary(self.config)
            matricesToLoad = h5u.checkMatricesDataBase(
                os.environ["SHESHA_ROOT"] + "/data/", self.config, param_dict)

        gpudevice = self.ui.wao_deviceNumber.value()
        self.ui.wao_deviceNumber.setDisabled(True)
        print "-> using GPU%d" % gpudevice
#        self.c = ch.naga_context()
#        self.c.set_activeDevice(device)
        if not self.c:
            #self.c = ch.naga_context(gpudevice)
            self.c=ch.naga_context(devices=np.array([4,5,6,7], dtype=np.int32))

        self.wfs, self.tel = ao.wfs_init(self.config.p_wfss, self.config.p_atmos, self.config.p_tel,
                                         self.config.p_geom, self.config.p_target, self.config.p_loop,
                                         self.config.p_dms)

        self.atm = ao.atmos_init(self.c, self.config.p_atmos, self.config.p_tel,
                                 self.config.p_geom, self.config.p_loop,
                                 self.config.p_wfss, self.config.p_target,
                                 clean=clean, load=matricesToLoad)
        self.ui.wao_atmosDimScreen.setText(
            str(self.config.p_atmos.dim_screens[0]))

        self.dms = ao.dm_init(
            self.config.p_dms, self.config.p_wfss, self.wfs, self.config.p_geom, self.config.p_tel)

        self.tar = ao.target_init(self.c, self.tel, self.config.p_target, self.config.p_atmos,
                                  self.config.p_geom, self.config.p_tel, self.config.p_wfss,
                                  self.wfs, self.config.p_dms)

        self.rtc = ao.rtc_init(self.tel, self.wfs, self.config.p_wfss, self.dms, self.config.p_dms,
                               self.config.p_geom, self.config.p_rtc, self.config.p_atmos,
                               self.atm, self.config.p_tel, self.config.p_loop, clean=clean, simul_name=simul_name, load=matricesToLoad)

        if(not clean):
            h5u.validDataBase(
                os.environ["SHESHA_ROOT"] + "/data/", matricesToLoad)

        for i in range(len(self.config.p_atmos.alt)):
            data = self.atm.get_screen(self.config.p_atmos.alt[i])
            cx, cy = self.circleCoords(self.config.p_geom.pupdiam/2, 1000, data.shape[0], data.shape[1])
            self.SRcircleAtmos[i] = pg.ScatterPlotItem(cx,cy, pen='r',size=1)
            self.p1.addItem(self.SRcircleAtmos[i])
            self.SRcircleAtmos[i].setPoints(cx, cy)
            self.SRcircleAtmos[i].hide()

        for i in range(len(self.config.p_wfss)):
            data = self.wfs.get_phase(i)
            cx, cy = self.circleCoords(self.config.p_geom.pupdiam/2, 1000, data.shape[0], data.shape[1])
            self.SRcircleWFS[i] = pg.ScatterPlotItem(cx,cy, pen='r',size=1)
            self.p1.addItem(self.SRcircleWFS[i])
            self.SRcircleWFS[i].setPoints(cx, cy)
            self.SRcircleWFS[i].hide()

        for i in range(len(self.config.p_dms)):
            dm_type = self.config.p_dms[i].type_dm
            alt = self.config.p_dms[i].alt
            data = self.dms.get_dm(dm_type, alt)
            cx, cy = self.circleCoords(self.config.p_geom.pupdiam/2, 1000, data.shape[0], data.shape[1])
            self.SRcircleDM[i] = pg.ScatterPlotItem(cx,cy, pen='r',size=1)
            self.p1.addItem(self.SRcircleDM[i])
            self.SRcircleDM[i].setPoints(cx, cy)
            self.SRcircleDM[i].hide()

        for i in range(self.config.p_target.ntargets):
            data = self.tar.get_phase(i)
            cx, cy = self.circleCoords(self.config.p_geom.pupdiam/2, 1000, data.shape[0], data.shape[1])
            self.SRcircleTarget[i] = pg.ScatterPlotItem(cx,cy, pen='r',size=1)
            self.p1.addItem(self.SRcircleTarget[i])
            self.SRcircleTarget[i].setPoints(cx, cy)
            self.SRcircleTarget[i].show()

        print "===================="
        print "init done"
        print "===================="
        print "objects initialized on GPU:"
        print "--------------------------------------------------------"
        print self.atm
        print self.wfs
        print self.dms
        print self.tar
        print self.rtc
        self.loaded = True
        self.updateDisplay()
        self.displayRtcMatrix()
        self.updatePlotWfs()
        self.p1.autoRange()

        self.ui.wao_PSFlogscale.clicked.connect(self.updateDisplay)
        self.ui.wao_run.setDisabled(False)
        self.ui.wao_next.setDisabled(False)
        self.ui.wao_openLoop.setDisabled(False)
        self.ui.wao_unzoom.setDisabled(False)
        self.ui.wao_resetSR.setDisabled(False)

    def resetDM(self):
        if(self.dms):
            ndm = self.ui.wao_selectDM.currentIndex()
            if(ndm > -1):
                self.dms.resetdm(
                    str(self.ui.wao_dmTypeSelector.currentText()), self.ui.wao_dmAlt.value())
                self.updateDisplay()
                print "DM "+str(ndm)+" reset"
            else:
                print "Invalid DM : please select a DM to reset"
        else:
            print "There is not any dm to reset"

    def BttCommand(self):
        if(self.rtc):
            nfilt = int(self.ui.wao_filterBtt.value())
            ao.command_on_Btt(self.rtc, nfilt)
            print "Loop is commanded from Btt basis now"

    def updatePlotWfs(self):
        RASC = 180./np.pi * 3600.
        typeText = str(self.ui.wao_wfs_plotSelector.currentText())
        n = self.ui.wao_selectWfs.currentIndex()
        self.ui.wao_wfsWindow.canvas.axes.clear()
        ax = self.ui.wao_wfsWindow.canvas.axes
        if(self.config.p_wfss[n].type_wfs == "pyrhr" and typeText == "Pyramid mod. pts" and self.loaded):
            scale_fact = 2*np.pi/self.config.p_wfss[n]._Nfft * \
                         self.config.p_wfss[n].Lambda*1e-6 / self.config.p_tel.diam / self.config.p_wfss[n]._qpixsize * RASC
            cx = self.config.p_wfss[n]._pyr_cx/scale_fact
            cy = self.config.p_wfss[n]._pyr_cy/scale_fact
            ax.scatter(cx, cy)


    def displayRtcMatrix(self):
        if not self.loaded:
            # print " widget not fully initialized"
            return

        data = None
        if(self.rtc):
            type_matrix = str(self.ui.wao_selectRtcMatrix.currentText())
            if(type_matrix == "imat"):
                data = self.rtc.get_imat(0)
            elif(type_matrix == "cmat"):
                data = self.rtc.get_cmat(0)
            elif(type_matrix == "Eigenvalues"):
                if(self.config.p_controllers[0].type_control == "ls" or self.config.p_controllers[0].type_control == "mv"):
                    data = self.rtc.getEigenvals(0)
            elif(type_matrix == "Cmm" and self.config.p_controllers[0].type_control == "mv"):
                tmp = self.rtc.get_cmm(0)
                ao.doTomoMatrices(0, self.rtc, self.config.p_wfss,
                                  self.dms, self.atm, self.wfs,
                                  self.config.p_rtc, self.config.p_geom,
                                  self.config.p_dms, self.config.p_tel, self.config.p_atmos)
                data = self.rtc.get_cmm(0)
                self.rtc.set_cmm(0, tmp)
            elif(type_matrix == "Cmm inverse" and self.config.p_controllers[0].type_control == "mv"):
                data = self.rtc.get_cmm(0)
            elif(type_matrix == "Cmm eigen" and self.config.p_controllers[0].type_control == "mv"):
                data = self.rtc.getCmmEigenvals(0)
            elif(type_matrix == "Cphim" and self.config.p_controllers[0].type_control == "mv"):
                data = self.rtc.get_cphim(0)

            if(data is not None):
                self.ui.wao_rtcWindow.canvas.axes.clear()
                ax = self.ui.wao_rtcWindow.canvas.axes
                if(len(data.shape) == 2):
                    self.ui.wao_rtcWindow.canvas.axes.matshow(
                        data, aspect="auto", origin="lower", cmap="gray")
                elif(len(data.shape) == 1):
                    # TODO : plot it properly, interactivity ?
                    self.ui.wao_rtcWindow.canvas.axes.plot(
                        range(len(data)), data, color="black")
                    ax.set_yscale('log')
                    if(type_matrix == "Eigenvalues"):
                        #    major_ticks = np.arange(0, 101, 20)
                        #    minor_ticks = np.arange(0, 101, 5)

                        #    self.ui.wao_rtcWindow.canvas.axes.set_xticks(major_ticks)
                        #    self.ui.wao_rtcWindow.canvas.axes.set_xticks(minor_ticks, minor=True)
                        #    self.ui.wao_rtcWindow.canvas.axes.set_yticks(major_ticks)
                        #    self.ui.wao_rtcWindow.canvas.axes.set_yticks(minor_ticks, minor=True)

                        # and a corresponding grid

                        self.ui.wao_rtcWindow.canvas.axes.grid(which='both')

                        # or if you want differnet settings for the grids:
                        self.ui.wao_rtcWindow.canvas.axes.grid(
                            which='minor', alpha=0.2)
                        self.ui.wao_rtcWindow.canvas.axes.grid(
                            which='major', alpha=0.5)
                        nfilt = self.rtc.get_nfiltered(0, self.config.p_rtc)
                        self.ui.wao_rtcWindow.canvas.axes.plot(
                            [nfilt - 0.5, nfilt - 0.5], [data.min(), data.max()], color="red", linestyle="dashed")
                        if(nfilt > 0):
                            self.ui.wao_rtcWindow.canvas.axes.scatter(np.arange(0, nfilt, 1), data[
                                                                      0:nfilt], color="red")  # TODO : plot it properly, interactivity ?
                            self.ui.wao_rtcWindow.canvas.axes.scatter(np.arange(nfilt, len(data), 1), data[
                                                                      nfilt:], color="blue")  # TODO : plot it properly, interactivity ?
                            tt = "%d modes Filtered" % nfilt
                            # ax.text(nfilt + 2, data[nfilt-1], tt)
                            ax.text(0.5, 0.2, tt, horizontalalignment='center',
                                    verticalalignment='center', transform=ax.transAxes)

                self.ui.wao_rtcWindow.canvas.draw()

    def circleCoords(self, ampli, npts, datashape0, datashape1):
        #ampli = self.config.p_geom.pupdiam/2
        #npts = 100
        cx = ampli * np.sin((np.arange(npts) + 1) * 2. * np.pi / npts) + datashape0/2
        cy = ampli * np.cos((np.arange(npts) + 1) * 2. * np.pi / npts) + datashape1/2
        return cx ,cy

    def computeDMrange(self, numdm, numwfs, push4imat=None, unitpervolt=None):
        i = numdm
        if(push4imat is None or push4imat==0):
            push4imat = self.config.p_dms[i].push4imat
        if(unitpervolt is None or unitpervolt==0):
            unitpervolt = self.config.p_dms[i].unitpervolt

        actuPushInMicrons = push4imat*unitpervolt
        coupling = self.config.p_dms[i].coupling
        a = coupling*actuPushInMicrons
        b = 0
        c = actuPushInMicrons
        d = coupling*actuPushInMicrons
        if(self.config.p_dms[i].type_dm is not "tt"):
            dist = self.config.p_tel.diam
        else:
            dist = self.config.p_tel.diam/self.config.p_wfss[numwfs].nxsub
        Delta = (1e-6*( ((c+d)/2) - ((a+b)/2)))
        actuPushInArcsecs = 206265.*Delta/dist
        return actuPushInArcsecs

    def updateDisplay(self):
        if not self.loaded:
            # print " widget not fully initialized"
            return

        data = None
        if not self.displayLock.acquire(False):
            # print " Display locked"
            return
        else:
            try:
                if(self.ui.wao_Display.isChecked()):
                    self.ui.wao_rtcWindowMPL.hide()
                    self.ui.wao_pgwindow.show()
                    if(self.SRCrossX and (self.imgType in ["Phase - Target", "Phase - DM", "Phase - Atmos", "Phase - WFS", "Spots - WFS", "Centroids - WFS", "Slopes - WFS"])):
                        self.SRCrossX.hide()
                        self.SRCrossY.hide()

                    #if(self.SRcircle and (self.imgType in ["Spots - WFS", "Centroids - WFS", "Slopes - WFS","PSF SE","PSF LE"])):
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
                            data = self.atm.get_screen(
                                self.config.p_atmos.alt[self.numberSelected])
                            if(self.imgType != self.currentViewSelected):
                                self.p1.setRange(
                                    xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType
                            self.SRcircleAtmos[self.numberSelected].show()

                    if(self.wfs):
                        if(self.imgType == "Phase - WFS"):
                            data = self.wfs.get_phase(self.numberSelected)
                            if(self.imgType != self.currentViewSelected):
                                self.p1.setRange(
                                    xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType
                            self.SRcircleWFS[self.numberSelected].show()

                        if(self.imgType == "Spots - WFS"):
                            if(self.config.p_wfss[self.numberSelected].type_wfs == "sh"):
                                data = self.wfs.get_binimg(self.numberSelected)
                            elif(self.config.p_wfss[self.numberSelected].type_wfs == "pyr"):
                                data = self.wfs.get_pyrimg(self.numberSelected)
                            if(self.imgType != self.currentViewSelected):
                                self.p1.setRange(
                                    xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType

                        if(self.imgType == "Pyrimg - LR"):
                            if(self.config.p_wfss[self.numberSelected].type_wfs == "pyrhr"):
                                data = self.wfs.get_pyrimg(self.numberSelected)
                            if(self.imgType != self.currentViewSelected):
                                self.p1.setRange(
                                    xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType

                        if(self.imgType == "Pyrimg - HR"):
                            if(self.config.p_wfss[self.numberSelected].type_wfs == "pyrhr"):
                                data = self.wfs.get_pyrimghr(self.numberSelected)
                            if(self.imgType != self.currentViewSelected):
                                self.p1.setRange(
                                    xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType

                        if(self.imgType == "Centroids - WFS"):
                            if(not self.ui.wao_rtcWindowMPL.isVisible()):
                                self.ui.wao_rtcWindowMPL.show()
                                self.ui.wao_pgwindow.hide()
                            self.ui.wao_rtcWindowMPL.canvas.axes.clear()
                            # retrieving centroids
                            centroids = self.rtc.getCentroids(0)
                            nvalid = [
                                2 * o._nvalid for o in self.config.p_wfss]
                            ind = np.sum(nvalid[:self.numberSelected])
                            x, y, vx, vy = tools.plsh(centroids[ind:ind + nvalid[self.numberSelected]], self.config.p_wfss[
                                                      self.numberSelected].nxsub, self.config.p_tel.cobs, returnquiver=True)  # Preparing mesh and vector for display
                            self.ui.wao_rtcWindowMPL.canvas.axes.quiver(
                                x, y, vx, vy, pivot='mid')
                            self.ui.wao_rtcWindowMPL.canvas.draw()
                            self.currentViewSelected = self.imgType

                            return
                        if(self.imgType == "Slopes - WFS"):
                            if(not self.ui.wao_rtcWindowMPL.isVisible()):
                                self.ui.wao_rtcWindowMPL.show()
                                self.ui.wao_pgwindow.hide()
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
                            data = self.tar.get_phase(self.numberSelected)
                            if(self.imgType != self.currentViewSelected):
                                self.p1.setRange(
                                    xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType
                            self.SRcircleTarget[self.numberSelected].show()

                        if(self.imgType == "PSF SE"):
                            data = self.tar.get_image(
                                self.numberSelected, "se")
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
                                self.SRCrossX.show()
                                self.SRCrossY.show()
                                self.p1.setRange(xRange=(data.shape[0] / 2 + 0.5 - zoom,
                                                         data.shape[0] / 2 + 0.5 + zoom),
                                                 yRange=(data.shape[1] / 2 + 0.5 - zoom,
                                                         data.shape[1] / 2 + 0.5 + zoom),)
                            self.currentViewSelected = self.imgType

                        if(self.imgType == "PSF LE"):
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
                self.displayLock.release()

    def mainLoop(self):
        if not self.loopLock.acquire(False):
            # print " Display locked"
            return
        else:
            try:
                start = time.time()
                self.atm.move_atmos()
                if(self.config.p_controllers[0].type_control == "geo"):
                    for t in range(self.config.p_target.ntargets):
                        self.tar.atmos_trace(t, self.atm, self.tel)
                        self.rtc.docontrol_geo(0, self.dms, self.tar, 0)
                        self.rtc.applycontrol(0, self.dms)
                        self.tar.dmtrace(0, self.dms)
                else:
                    for t in range(self.config.p_target.ntargets):
                        self.tar.atmos_trace(t, self.atm, self.tel)
                        self.tar.dmtrace(t, self.dms)
                    for w in range(len(self.config.p_wfss)):
                        self.wfs.sensors_trace(w, "all", self.tel, self.atm, self.dms)
                        self.wfs.sensors_compimg(w)

                    self.rtc.docentroids(0)
                    self.rtc.docontrol(0)
                    self.rtc.applycontrol(0, self.dms)

                signal_le = ""
                signal_se = ""
                for t in range(self.config.p_target.ntargets):
                    SR = self.tar.get_strehl(t)
                    signal_se += "%1.2f   " % SR[0]
                    signal_le += "%1.2f   " % SR[1]

                loopTime = time.time() - start
                if(self.RTDisplay):
                    t = 1 / float(self.RTDFreq) - loopTime  # Limit loop frequency
                    if t > 0:
                        time.sleep(t)  # Limit loop frequency
                    self.updateDisplay()  # Update GUI plots
                else:
                    freqLimit = 1 / 250.
                    if loopTime < freqLimit:
                        time.sleep(freqLimit - loopTime)  # Limit loop frequency
                currentFreq = 1 / (time.time() - start)

                if(self.RTDisplay):
                    self.ui.wao_strehlSE.setText(signal_se)
                    self.ui.wao_strehlLE.setText(signal_le)
                    self.ui.wao_currentFreq.setValue(currentFreq)
                if(time.time() - self.startTime > 1):
                    self.printInPlace("iter #%d SR: (L.E, S.E.)= %s, %srunning at %4.1fHz (real %4.1fHz)" % (
                        self.iter, signal_le, signal_se, currentFreq, 1 / loopTime))
                    self.startTime = time.time()
                    # This seems to trigger the GUI and keep it responsive
                    time.sleep(.01)
                self.iter += 1
            finally:
                self.loopLock.release()

    def printInPlace(self, text):
        # This seems to trigger the GUI and keep it responsive
        print "\r" + text,
        sys.stdout.flush()
        # sys.stdout.write(text)

    def run(self):
        # print "Loop started"
        self.c.set_activeDeviceForce(0)
        self.stop = False
        self.startTime = time.time()
        i = 0
        print "LOOP STARTED FOR %d iterations" % self.nbiters
        while i<=self.nbiters:
            self.mainLoop()
            i+=1
            if(self.stop):
                break
        self.ui.wao_run.setChecked(False)
#         print "Loop stopped"

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    wao = widgetAOWindow()
    wao.show()
    app.setStyle('cleanlooks')
    # ao.imat_init(0, wao.rtc, wao.config.p_rtc, wao.dms, wao.wfs, wao.config.p_wfss, wao.config.p_tel)
    # app.exec_()
