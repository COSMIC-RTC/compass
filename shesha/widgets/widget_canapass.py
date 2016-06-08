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
import glob
import tools
import hdf5_utils as h5u
import threading
from PyQt4.uic import loadUiType
from PyQt4 import QtGui
from functools import partial
import subprocess

sys.path.insert(0, os.environ["SHESHA_ROOT"] + "/data/par/")
WindowTemplate, TemplateBaseClass = loadUiType(
    os.environ["SHESHA_ROOT"] + "/widgets/widget_canapass.ui")


"""
low levels debugs:
gdb --args python -i widget_canapass.py

"""


class widgetAOWindow(TemplateBaseClass):

    def __init__(self):
        TemplateBaseClass.__init__(self)

        self.ui = WindowTemplate()
        self.ui.setupUi(self)

        ##############################################################
        #######       ATTRIBUTES   #######################
        #############################################################
        self.config = None
        self.c = None
        self.tel = None
        self.wfs = None
        self.rtc = None
        self.atm = None
        self.tar = None
        self.dms = None

        self.brama_flag = 1
        self.see_atmos = 0
        self.displayLock = threading.Lock()
        self.iter = 0
        self.loaded = False
        self.stop = False
        self.startTime = 0
        self.loop = None
        self.assistant = None

        ##############################################################
        #######       PYQTGRAPH WINDOW INIT   #######################
        #############################################################

        self.img = pg.ImageItem(border='w')  # create image area
        # self.p1 = self.ui.wao_pgwindow.addPlot() # create pyqtgraph plot area
        self.p1 = self.ui.wao_pgwindow.addViewBox()
        self.p1.setAspectLocked(True)

        self.p1.addItem(self.img)  # Put image in plot area
        self.hist = pg.HistogramLUTItem()  # Create an histogram
        self.hist.setImageItem(self.img)  # Compute histogram from img
        self.ui.wao_pgwindow.addItem(self.hist)
        self.hist.autoHistogramRange()  # init levels
        self.hist.setMaximumWidth(100)

        ##############################################################
        #######       CONNECTED BUTTONS  #######################
        #############################################################
        # Default path for config files
        self.defaultParPath = os.environ[
            "SHESHA_ROOT"] + "/data/par/par4bench/"
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
        self.ui.wao_Display.clicked.connect(self.updateFrameRate)
        self.ui.wao_frameRate.valueChanged.connect(self.updateFrameRate)
        self.ui.wao_rtcWindowMPL.hide()
        self.RTDisplay = self.ui.wao_Display.isChecked()
        self.RTDFreq = self.ui.wao_frameRate.value()
        self.ui.wao_PSFlogscale.clicked.connect(self.updateDisplay)
        self.ui.wao_atmosphere.setCheckable(True)
        self.ui.wao_atmosphere.clicked[bool].connect(self.set_atmos)

        self.SRcircleAtmos = {}
        self.SRcircleWFS = {}
        self.SRcircleDM = {}
        self.SRcircleTarget = {}

        self.addConfigFromFile(
            filepath=os.environ["SHESHA_ROOT"] + "/data/par/canapass.py")
        self.InitConfig()

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

        ##############################################################
        ##################       METHODS      #######################
        #############################################################

    def updateFrameRate(self):
        if(self.ui.wao_Display.isChecked()):
            self.RTDisplay = True
        else:
            self.RTDisplay = False
        self.RTDFreq = self.ui.wao_frameRate.value()

    def addConfigFromFile(self, filepath=None):
        if filepath is None:
            filepath = str(QtGui.QFileDialog(directory=self.defaultParPath).getOpenFileName(
                self, "Select parameter file", "", "parameters file (*.py);;hdf5 file (*.h5);;all files (*)"))
        self.configpath = filepath
        filename = filepath.split('/')[-1]
        if(filepath.split('.')[-1] == "py"):
            pathfile = filepath.split(filename)[-1]
            # if (pathfile not in sys.path):
            sys.path.insert(0, pathfile)
            exec("import %s as config" % filename.split(".py")[0])
            sys.path.remove(pathfile)
        elif(filepath.split('.')[-1] == "h5"):
            sys.path.insert(0, self.defaultParPath)
            import scao_16x16_8pix as config
            sys.path.remove(self.defaultParPath)
            h5u.configFromH5(filepath, config)
        else:
            raise ValueError("Parameter file extension must be .py or .h5")
        self.config = config
        self.updateNumberSelector(textType=self.imgType)

    def aoLoopClicked(self, pressed):
        if(pressed):
            self.loop = threading.Thread(target=self.run)
            self.loop.start()
        else:
            self.stop = True
            self.loop.join()

    def aoLoopFinished(self):
        if not self.ui.wao_run.isChecked:
            self.ui.wao_run.click()

    def setNumberSelection(self):
        if(self.ui.wao_selectNumber.currentIndex() > -1):
            self.numberSelected = self.ui.wao_selectNumber.currentIndex()
        else:
            self.numberSelected = 0
        self.updateDisplay()

    def updateNumberSelector(self, textType=None):
        if(textType == None):
            textType = str(self.ui.wao_selectScreen.currentText())
        self.imgType = textType
        self.ui.wao_selectNumber.clear()
        if(textType == "Phase - Atmos"):
            n = self.config.p_atmos.nscreens
        if(textType == "Phase - WFS" or textType == "Spots - WFS" or textType == "Centroids - WFS" or textType == "Slopes - WFS"):
            n = len(self.config.p_wfss)
        if(textType == "Phase - Target" or textType == "PSF LE" or textType == "PSF SE"):
            n = self.config.p_target.ntargets
        if(textType == "Phase - DM"):
            n = len(self.config.p_dms)
        self.ui.wao_selectNumber.addItems([str(i) for i in range(n)])
        self.updateDisplay()

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
        if(simul_name == ""):
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
            self.c = ch.naga_context(gpudevice)

        self.wfs, self.tel = ao.wfs_init(self.config.p_wfss, self.config.p_atmos, self.config.p_tel,
                                         self.config.p_geom, self.config.p_target, self.config.p_loop,
                                         self.config.p_dms)

        self.atm = ao.atmos_init(self.c, self.config.p_atmos, self.config.p_tel,
                                 self.config.p_geom, self.config.p_loop,
                                 self.config.p_wfss, self.config.p_target, clean=clean, load=matricesToLoad)

        self.dms = ao.dm_init(
            self.config.p_dms, self.config.p_wfss, self.wfs, self.config.p_geom, self.config.p_tel)

        self.tar = ao.target_init(self.c, self.tel, self.config.p_target, self.config.p_atmos,
                                  self.config.p_geom, self.config.p_tel, self.config.p_wfss,
                                  self.wfs, self.config.p_dms, brama=self.brama_flag)

        self.rtc = ao.rtc_init(self.tel, self.wfs, self.config.p_wfss, self.dms, self.config.p_dms,
                               self.config.p_geom, self.config.p_rtc, self.config.p_atmos,
                               self.atm, self.config.p_tel, self.config.p_loop, clean=clean, simul_name=simul_name, 
                               load=matricesToLoad, brama=self.brama_flag, brama_tar=self.tar)
        self.rtc.set_openloop(0, 1)

        if(simul_name is not ""):
            h5u.validDataBase(
                os.environ["SHESHA_ROOT"] + "/data/", matricesToLoad)

        for i in range(len(self.config.p_atmos.alt)):
            data = self.atm.get_screen(self.config.p_atmos.alt[i])
            cx, cy = self.circleCoords(self.config.p_geom.pupdiam / 2, 100, data.shape[0], data.shape[1])
            self.SRcircleAtmos[i] = pg.ScatterPlotItem(cx, cy, pen='r')
            self.p1.addItem(self.SRcircleAtmos[i])
            self.SRcircleAtmos[i].setPoints(cx, cy)
            self.SRcircleAtmos[i].hide()

        for i in range(len(self.config.p_wfss)):
            data = self.wfs.get_phase(i)
            cx, cy = self.circleCoords(self.config.p_geom.pupdiam / 2, 100, data.shape[0], data.shape[1])
            self.SRcircleWFS[i] = pg.ScatterPlotItem(cx, cy, pen='r')
            self.p1.addItem(self.SRcircleWFS[i])
            self.SRcircleWFS[i].setPoints(cx, cy)
            self.SRcircleWFS[i].hide()

        for i in range(len(self.config.p_dms)):
            dm_type = self.config.p_dms[i].type_dm
            alt = self.config.p_dms[i].alt
            data = self.dms.get_dm(dm_type, alt)
            cx, cy = self.circleCoords(self.config.p_geom.pupdiam / 2, 100, data.shape[0], data.shape[1])
            self.SRcircleDM[i] = pg.ScatterPlotItem(cx, cy, pen='r')
            self.p1.addItem(self.SRcircleDM[i])
            self.SRcircleDM[i].setPoints(cx, cy)
            self.SRcircleDM[i].hide()
        
        for i in range(self.config.p_target.ntargets):
            data = self.tar.get_phase(i)
            cx, cy = self.circleCoords(self.config.p_geom.pupdiam / 2, 100, data.shape[0], data.shape[1])
            self.SRcircleTarget[i] = pg.ScatterPlotItem(cx, cy, pen='r')
            self.p1.addItem(self.SRcircleTarget[i])
            self.SRcircleTarget[i].setPoints(cx, cy)
            self.SRcircleTarget[i].show()

        self.loaded = True
        self.updateDisplay()
        self.p1.autoRange()

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

        self.ui.wao_PSFlogscale.clicked.connect(self.updateDisplay)
        self.ui.wao_run.setDisabled(False)
        self.ui.wao_next.setDisabled(False)
        self.ui.wao_unzoom.setDisabled(False)
        self.ui.wao_resetSR.setDisabled(False)

    def circleCoords(self, ampli, npts, datashape0, datashape1):
        # ampli = self.config.p_geom.pupdiam/2
        # npts = 100
        cx = ampli * np.sin((np.arange(npts) + 1) * 2. * np.pi / npts) + datashape0 / 2
        cy = ampli * np.cos((np.arange(npts) + 1) * 2. * np.pi / npts) + datashape1 / 2
        return cx , cy

    def resetDM(self):
        if(self.dms):
            ndm = self.ui.wao_selectDM.currentIndex()
            if(ndm > -1):
                self.dms.resetdm(
                    str(self.ui.wao_dmTypeSelector.currentText()), self.ui.wao_dmAlt.value())
                self.updateDisplay()
            else:
                print "Invalid DM : please select a DM to reset"
        else:
            print "There is not any dm to reset"

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
        
                    # if(self.SRcircle and (self.imgType in ["Spots - WFS", "Centroids - WFS", "Slopes - WFS","PSF SE","PSF LE"])):
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
                        if(self.imgType == "Centroids - WFS"):
                            if(not self.ui.wao_rtcWindowMPL.isVisible()):
                                self.ui.wao_rtcWindowMPL.show()
                                self.ui.wao_pgwindow.hide()
                            self.ui.wao_rtcWindowMPL.canvas.axes.clear()
                            # retrieving centroids
                            centroids = self.rtc.getCentroids(0)
                            nvalid = [2 * o._nvalid for o in self.config.p_wfss]
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
                            wao.ui.wao_rtcWindowMPL.canvas.axes.quiver(
                                x, y, vx, vy, pivot='mid')
                            wao.ui.wao_rtcWindowMPL.canvas.draw()
                            self.currentViewSelected = self.imgType
                            return
        
                    if(self.dms):
                        if(self.imgType == "Phase - DM"):
                            dm_type = self.config.p_dms[self.numberSelected].type_dm
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
                            data = self.tar.get_image(self.numberSelected, "se")
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
                            data = self.tar.get_image(self.numberSelected, "le")
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
                    if (data is not None):
                        autoscale = self.ui.wao_autoscale.isChecked()
                        if(autoscale):
                            self.hist.setLevels(data.min(), data.max())  # inits levels
                        self.img.setImage(data, autoLevels=autoscale)
#                         self.p1.autoRange()
            finally:
                self.displayLock.release()

    def mainLoop(self):
        start = time.time()
        self.atm.move_atmos()
        for t in range(self.config.p_target.ntargets):
            if wao.see_atmos:
                self.tar.atmos_trace(t,self.atm,self.tel)
            else:
                self.tar.reset_phase(t)
            self.tar.dmtrace(t, self.dms)
        for w in range(len(self.config.p_wfss)):
            if wao.see_atmos:
                self.wfs.sensors_trace(w,"atmos",self.tel,self.atm,self.dms)
            else:
                self.wfs.reset_phase(w)
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

        if(wao.brama_flag):
            self.rtc.publish() #rtc_publish, g_rtc;
            self.tar.publish()

        if(self.RTDisplay):
            self.ui.wao_strehlSE.setText(signal_se)
            self.ui.wao_strehlLE.setText(signal_le)
            self.ui.wao_currentFreq.setValue(currentFreq)
        else:
            if(time.time() - self.startTime > 1):
                self.printInPlace("iter #%d SR: (L.E, S.E.)= %s, %srunning at %4.1fHz (real %4.1fHz)" % (
                    self.iter, signal_le, signal_se, currentFreq, 1 / loopTime))
                self.ui.wao_strehlSE.setText(signal_se)
                self.ui.wao_strehlLE.setText(signal_le)
                self.ui.wao_currentFreq.setValue(currentFreq)
                self.startTime = time.time()
                time.sleep(.01)
        self.iter += 1

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
        while True:
            self.mainLoop()
            if(self.stop):
                break
#         print "Loop stopped"

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    wao = widgetAOWindow()
    wao.show()
    # app.connect(wao.ui._quit,QtCore.SIGNAL("clicked()"),app,QtCore.SLOT("quit()"))
    app.setStyle('cleanlooks')
