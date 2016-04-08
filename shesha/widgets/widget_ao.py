
import cProfile
import pstats as ps

import sys, os
import numpy as np
import naga as ch
import shesha as ao
import time
import matplotlib.pyplot as pl
import pyqtgraph as pg
import glob
import tools
import hdf5_utils as h5u
sys.path.insert(0, os.environ["SHESHA_ROOT"]+"/data/par/")

from PyQt4.uic import loadUiType
from PyQt4 import QtCore, QtGui

from functools import partial
import time
WindowTemplate,TemplateBaseClass=loadUiType(os.environ["SHESHA_ROOT"]+"/widgets/widget_ao.ui")

import threading

"""
low levels debugs:
gdb --args python -i widget_ao.py

"""
class widgetAOWindow(TemplateBaseClass):
    def __init__(self, device=None):
        TemplateBaseClass.__init__(self)

        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        self.displayLock = threading.Lock()

        ##############################################################
        #######       ATTRIBUTES   #######################
        #############################################################
        self.configpath = None
        self.config = None
        self.c = None
        self.tel = None
        self.wfs = None
        self.rtc = None
        self.atm = None
        self.tar = None
        self.dms = None
        self.aoLoopThread = None
        self.mainLoop = [self.tel,self.atm,self.wfs,self.rtc,self.tar,self.dms]

        ##############################################################
        #######       PYQTGRAPH WINDOW INIT   #######################
        #############################################################

        self.img = pg.ImageItem(border='w') # create image area
        #self.p1 = self.ui.wao_pgwindow.addPlot() # create pyqtgraph plot area
        self.p1 = self.ui.wao_pgwindow.addViewBox()
        self.p1.setAspectLocked(True)

        self.p1.addItem(self.img) # Put image in plot area


        self.hist = pg.HistogramLUTItem() #Create an histogram
        self.hist.setImageItem(self.img) # Compute histogram from img
        self.ui.wao_pgwindow.addItem(self.hist)
        self.hist.autoHistogramRange() # init levels
        self.hist.setMaximumWidth(100)

        ##############################################################
        #######       CONNECTED BUTTONS  #######################
        #############################################################
        self.defaultParPath = os.environ["SHESHA_ROOT"]+"/data/par/par4bench/" # Default path for config files
        self.ui.wao_loadConfig.clicked.connect(self.loadConfig)
        self.loadDefaultConfig()
        self.ui.wao_init.clicked.connect(self.InitConfig)
        self.ui.wao_run.setCheckable(True)
        self.ui.wao_run.clicked[bool].connect(self.aoLoopClicked)
        self.ui.wao_openLoop.setCheckable(True)
        self.ui.wao_openLoop.clicked[bool].connect(self.aoLoopOpen)
        self.ui.wao_next.clicked.connect(self.nextClicked)
        self.imgType = str(self.ui.wao_selectScreen.currentText())
        self.ui.wao_configFromFile.clicked.connect(self.addConfigFromFile)
        self.ui.wao_unzoom.clicked.connect(self.p1.autoRange)
        self.ui.wao_selectScreen.currentIndexChanged.connect(partial(self.updateNumberSelector,textType=None))
        self.ui.wao_selectNumber.currentIndexChanged.connect(self.setNumberSelection)
        self.ui.wao_selectAtmosLayer.currentIndexChanged.connect(self.setLayerSelection)
        self.ui.wao_selectWfs.currentIndexChanged.connect(self.setWfsSelection)
        self.ui.wao_selectDM.currentIndexChanged.connect(self.setDmSelection)
        self.ui.wao_selectCentro.currentIndexChanged.connect(self.setCentroSelection)
        self.ui.wao_selectTarget.currentIndexChanged.connect(self.setTargetSelection)
        self.ui.wao_setAtmos.clicked.connect(self.setAtmosParams)
        self.ui.wao_setWfs.clicked.connect(self.setWfsParams)
        self.ui.wao_setDM.clicked.connect(self.setDmParams)
        self.ui.wao_setControl.clicked.connect(self.setRtcParams)
        self.ui.wao_setCentro.clicked.connect(self.setRtcParams)
        self.ui.wao_setTelescope.clicked.connect(self.setTelescopeParams)
        self.ui.wao_resetDM.clicked.connect(self.resetDM)
        self.ui.wao_selectRtcMatrix.currentIndexChanged.connect(self.displayRtcMatrix)
        self.ui.wao_rtcWindowMPL.hide()
        self.ui.wao_Display.clicked.connect(self.updateFrameRate)
        self.ui.wao_frameRate.valueChanged.connect(self.updateFrameRate)
        self.RTDisplay= self.ui.wao_Display.isChecked()
        self.RTDFreq = self.ui.wao_frameRate.value()
        self.ui.wao_PSFlogscale.clicked.connect(self.updateDisplay)
        self.ui.wao_resetSR.clicked.connect(self.resetSR)

        # Create Loop thread
        self.aoLoopThread = aoLoopThread(self.mainLoop, self.config, self.img, self.ui.wao_strehlSE, self.ui.wao_strehlLE, 1, self.hist, self.RTDisplay, self.RTDFreq, self.updateDisplay)
        self.connect(self.aoLoopThread, QtCore.SIGNAL("currentLoopFrequency(float)"), self.updateCurrentLoopFrequency)
        self.connect(self.aoLoopThread, QtCore.SIGNAL("currentSRSE(QString)"), self.updateSRSE)
        self.connect(self.aoLoopThread, QtCore.SIGNAL("currentSRLE(QString)"), self.updateSRLE)

        self.connect(self.aoLoopThread,QtCore.SIGNAL("finished()"),self.aoLoopFinished)

    def resetSR(self):
        tarnum = self.ui.wao_resetSR_tarNum.value()
        print "reset SR on target %d" % tarnum
        self.tar.reset_strehl(tarnum)

    def setDevice(self):
        self.c.set_activeDevice(self.ui.wao_deviceNumber.value())

    def manually_destroy(self):
        if(hasattr(self,"mainLoop")):
            del self.mainLoop
        if(hasattr(self,"atm")):
            del self.atm
            del self.aoLoopThread.atm
        if(hasattr(self,"tel")):
            del self.tel
            del self.aoLoopThread.tel
        if(hasattr(self,"wfs")):
            del self.wfs
            del self.aoLoopThread.wfs
        if(hasattr(self,"rtc")):
            del self.rtc
            del self.aoLoopThread.rtc
        if(hasattr(self,"tar")):
            del self.tar
            del self.aoLoopThread.tar
        if(hasattr(self,"dms")):
            del self.dms
            del self.aoLoopThread.dms
        self.atm = None
        self.tel = None
        self.wfs = None
        self.rtc = None
        self.tar = None
        self.dms = None
        self.mainLoop = [self.tel,self.atm,self.wfs,self.rtc,self.tar,self.dms]

    def updateSRSE(self, SRSE):
        self.ui.wao_strehlSE.setText(SRSE)

    def updateSRLE(self, SRLE):
        self.ui.wao_strehlLE.setText(SRLE)

    def updateCurrentLoopFrequency(self, freq):
        self.ui.wao_currentFreq.setValue(freq)

    def closeEvent(self, event):

        reply = QtGui.QMessageBox.question(self, 'Message',
            "Are you sure to quit?", QtGui.QMessageBox.Yes |
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
            self.aoLoopThread.terminate()
            self.destroy() # TODO: segfault

            quit()
            #sys.exit()
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
        if(self.aoLoopThread):
            self.aoLoopThread.RTDisplay = self.RTDisplay
            self.aoLoopThread.RTDFreq = self.RTDFreq

    def updateTelescopePanel(self):
        self.ui.wao_zenithAngle.setValue(self.config.p_geom.zenithangle)
        self.ui.wao_diamTel.setValue(self.config.p_tel.diam)
        self.ui.wao_cobs.setValue(self.config.p_tel.cobs)

    def updateDmPanel(self):
        ndm = self.ui.wao_selectDM.currentIndex()
        if(ndm < 0):
            ndm = 0
        self.ui.wao_numberofDMs.setText(str(len(self.config.p_dms)))
        self.ui.wao_dmTypeSelector.setCurrentIndex(self.ui.wao_dmTypeSelector.findText(self.config.p_dms[ndm].type_dm))
        self.ui.wao_dmAlt.setValue(self.config.p_dms[ndm].alt)
        self.ui.wao_dmNactu.setValue(self.config.p_dms[ndm].nact)
        self.ui.wao_dmUnitPerVolt.setValue(self.config.p_dms[ndm].unitpervolt)
        self.ui.wao_dmCoupling.setValue(self.config.p_dms[ndm].coupling)
        self.ui.wao_dmThresh.setValue(self.config.p_dms[ndm].thresh)

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

        # LGS panel
        if(self.config.p_wfss[nwfs].gsalt > 0):
            self.ui.wao_wfsIsLGS.setChecked(True)
            self.ui.wao_wfsGsAlt.setValue(self.config.p_wfss[nwfs].gsalt)
            self.ui.wao_wfsLLTx.setValue(self.config.p_wfss[nwfs].lltx)
            self.ui.wao_wfsLLTy.setValue(self.config.p_wfss[nwfs].llty)
            self.ui.wao_wfsLGSpower.setValue(self.config.p_wfss[nwfs].laserpower)
            self.ui.wao_wfsReturnPerWatt.setValue(self.config.p_wfss[nwfs].lgsreturnperwatt)
            self.ui.wao_wfsBeamSize.setValue(self.config.p_wfss[nwfs].beamsize)
            self.ui.wao_selectLGSProfile.setCurrentIndex(self.ui.wao_selectLGSProfile.findText(self.config.p_wfss[nwfs].proftype))

        else:
            self.ui.wao_wfsIsLGS.setChecked(False)

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
        self.ui.wao_windDirection.setValue(self.config.p_atmos.winddir[nscreen])
        if(self.config.p_atmos.dim_screens is not None):
            self.ui.wao_atmosDimScreen.setText(str(self.config.p_atmos.dim_screens[nscreen]))
        self.ui.wao_atmosWindow.canvas.axes.cla()
        width = (self.config.p_atmos.alt.max()/20.+0.1)/1000.
        self.ui.wao_atmosWindow.canvas.axes.barh(self.config.p_atmos.alt/1000.-width/2.,self.config.p_atmos.frac,width,color="blue")
        self.ui.wao_atmosWindow.canvas.draw()

    def updateRtcPanel(self):
        # Centroider panel
        ncentro = self.ui.wao_selectCentro.currentIndex()
        if(ncentro < 0):
            ncentro = 0
        self.ui.wao_centroTypeSelector.setCurrentIndex(self.ui.wao_centroTypeSelector.findText(self.config.p_centroiders[ncentro].type_centro))
        self.ui.wao_centroThresh.setValue(self.config.p_centroiders[ncentro].thresh)
        self.ui.wao_centroNbrightest.setValue(self.config.p_centroiders[ncentro].nmax)
        self.ui.wao_centroThresh.setValue(self.config.p_centroiders[ncentro].thresh)
        if(self.config.p_centroiders[ncentro].type_fct):
            self.ui.wao_centroFunctionSelector.setCurrentIndex(self.ui.wao_centroFunctionSelector.findText(self.config.p_centroiders[ncentro].type_fct))
        self.ui.wao_centroWidth.setValue(self.config.p_centroiders[ncentro].width)

        #Controller panel
        type_contro = self.config.p_controllers[0].type_control
        if(type_contro == "ls" and self.config.p_controllers[0].modopti==0):
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
        self.ui.wao_controlTTcond.setValue(self.config.p_controllers[0].maxcond) # TODO : TTcond

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
            self.ui.wao_targetWindow.canvas.axes.plot([w.xpos for w in self.config.p_wfss],[w.ypos for w in self.config.p_wfss],'o',color="green")
            xmax = np.max([xmax,np.max(np.abs([w.xpos for w in self.config.p_wfss]))])
            ymax = np.max([ymax,np.max(np.abs([w.ypos for w in self.config.p_wfss]))])
        self.ui.wao_targetWindow.canvas.axes.plot(self.config.p_target.xpos,self.config.p_target.ypos,'*',color="red")
        self.ui.wao_targetWindow.canvas.axes.set_xlim(-xmax-10,xmax+10)
        self.ui.wao_targetWindow.canvas.axes.set_ylim(-ymax-10,ymax+10)
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
        self.config.p_atmos.alt[nscreen]=self.ui.wao_atmosAlt.value()
        self.config.p_atmos.frac[nscreen]=self.ui.wao_atmosFrac.value()
        self.config.p_atmos.L0[nscreen]=self.ui.wao_atmosL0.value()
        self.config.p_atmos.windspeed[nscreen]=self.ui.wao_windSpeed.value()
        self.config.p_atmos.winddir[nscreen]=self.ui.wao_windDirection.value()
        print "New atmos parameters set"

    def setRtcParams(self):
        # Centroider params
        ncentro = self.ui.wao_selectCentro.currentIndex()
        if(ncentro < 0):
            ncentro = 0
        self.config.p_centroiders[ncentro].set_type(str(self.ui.wao_centroTypeSelector.currentText()))
        self.config.p_centroiders[ncentro].set_thresh(self.ui.wao_centroThresh.value())
        self.config.p_centroiders[ncentro].set_nmax(self.ui.wao_centroNbrightest.value())
        self.config.p_centroiders[ncentro].set_thresh(self.ui.wao_centroThresh.value())
        self.config.p_centroiders[ncentro].set_type_fct(str(self.ui.wao_centroFunctionSelector.currentText()))
        self.config.p_centroiders[ncentro].set_width(self.ui.wao_centroWidth.value())

        #Controller panel
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

        self.config.p_controllers[0].set_maxcond(self.ui.wao_controlCond.value())
        self.config.p_controllers[0].set_delay(self.ui.wao_controlDelay.value())
        self.config.p_controllers[0].set_gain(self.ui.wao_controlGain.value())
        #self.config.p_controllers[0].set_TTcond(self.ui.wao_controlTTcond.value()) # TODO : TTcond
        print "New rtc parameters set"

    def setWfsParams(self):
        nwfs = self.ui.wao_selectWfs.currentIndex()
        if(nwfs < 0):
            nwfs = 0
        self.config.p_wfss[nwfs].set_nxsub(self.ui.wao_wfsNxsub.value())
        self.config.p_wfss[nwfs].set_npix( self.ui.wao_wfsNpix.value())
        self.config.p_wfss[nwfs].set_pixsize( self.ui.wao_wfsPixSize.value())
        self.config.p_wfss[nwfs].set_xpos( self.ui.wao_wfsXpos.value())
        self.config.p_wfss[nwfs].set_ypos( self.ui.wao_wfsYpos.value())
        self.config.p_wfss[nwfs].set_fracsub( self.ui.wao_wfsFracsub.value())
        self.config.p_wfss[nwfs].set_Lambda( self.ui.wao_wfsLambda.value())
        self.config.p_wfss[nwfs].set_gsmag( self.ui.wao_wfsMagnitude.value())
        #TODO: find a way to correctly set zerop (limited by the maximum value allowed by the double spin box)
        self.config.p_wfss[nwfs].set_zerop(10**(self.ui.wao_wfsZp.value()))
        self.config.p_wfss[nwfs].set_optthroughput( self.ui.wao_wfsThrough.value())
        self.config.p_wfss[nwfs].set_noise( self.ui.wao_wfsNoise.value())

        # LGS params
        if(self.ui.wao_wfsIsLGS.isChecked()):
            self.config.p_wfss[nwfs].set_gsalt(self.ui.wao_wfsGsAlt.value())
            self.config.p_wfss[nwfs].set_lltx(self.ui.wao_wfsLLTx.value())
            self.config.p_wfss[nwfs].set_llty(self.ui.wao_wfsLLTy.value())
            self.config.p_wfss[nwfs].set_laserpower(self.ui.wao_wfsLGSpower.value())
            self.config.p_wfss[nwfs].set_lgsreturnperwatt(self.ui.wao_wfsReturnPerWatt.value())
            self.config.p_wfss[nwfs].set_beamsize(self.ui.wao_wfsBeamSize.value())
            self.config.p_wfss[nwfs].set_proftype(str(self.ui.wao_selectLGSProfile.currentText()))
        print "New wfs parameters set"

    def setDmParams(self):
        ndm = self.ui.wao_selectDM.currentIndex()
        if(ndm < 0):
            ndm = 0
        self.config.p_dms[ndm].set_type(str(self.ui.wao_dmTypeSelector.currentText()))
        self.config.p_dms[ndm].set_alt(self.ui.wao_dmAlt.value())
        self.config.p_dms[ndm].set_nact(self.ui.wao_dmNactu.value())
        self.config.p_dms[ndm].set_unitpervolt(self.ui.wao_dmUnitPerVolt.value())
        self.config.p_dms[ndm].set_coupling(self.ui.wao_dmCoupling.value())
        self.config.p_dms[ndm].set_thresh(self.ui.wao_dmThresh.value())
        print "New DM parameters set"

    def updateLayerSelection(self):
        self.ui.wao_selectAtmosLayer.clear()
        self.ui.wao_selectAtmosLayer.addItems([str(i) for i in range(self.config.p_atmos.nscreens)])

    def updateTargetSelection(self):
        self.ui.wao_selectTarget.clear()
        self.ui.wao_selectTarget.addItems([str(i) for i in range(self.config.p_target.ntargets)])

    def updateWfsSelection(self):
        self.ui.wao_selectWfs.clear()
        self.ui.wao_selectWfs.addItems([str(i) for i in range(len(self.config.p_wfss))])

    def updateDmSelection(self):
        self.ui.wao_selectDM.clear()
        self.ui.wao_selectDM.addItems([str(i) for i in range(len(self.config.p_dms))])

    def updateCentroSelection(self):
        self.ui.wao_selectCentro.clear()
        self.ui.wao_selectCentro.addItems([str(i) for i in range(len(self.config.p_centroiders))])

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
        filepath = str(QtGui.QFileDialog(directory=self.defaultParPath).getOpenFileName(self,"Select parameter file","","parameters file (*.py);;hdf5 file (*.h5);;all files (*)"))
        self.configpath = filepath
        filename = filepath.split('/')[-1]
        if(filepath.split('.')[-1] == "py"):
            self.ui.wao_selectConfig.addItem(filename,0)
            pathfile = filepath.split(filename)[0]
            if (pathfile not in sys.path):
                sys.path.insert(0, pathfile)
            exec("import %s as config" % filename.split(".py")[0])
            sys.path.remove(pathfile)
        elif(filepath.split('.')[-1] == "h5"):
            sys.path.insert(0,self.defaultParPath)
            import scao_16x16_8pix as config
            sys.path.remove(self.defaultParPath)
            h5u.configFromH5(filepath,config)
        else:
            print "Parameter file extension must be .py or .h5"
            return
        self.config = config
        self.aoLoopThread.config = config
        self.ui.wao_selectConfig.clear()
        self.ui.wao_selectConfig.addItem(filename)
        self.updateNumberSelector(textType=self.imgType)
        self.updatePanels()

    def nextClicked(self):
#        self.aoLoopThread = aoLoopThread(self.mainLoop, self.config, self.img, self.ui.wao_strehlSE, self.ui.wao_strehlLE, 1, self.imgType, self.numberSelected)
#        self.connect(self.aoLoopThread,QtCore.SIGNAL("finished()"),self.aoLoopFinished)
        self.ui.wao_strehlSE.clear()
        self.ui.wao_strehlLE.clear()
        self.aoLoopThread.framebyframe=1
        self.aoLoopThread.start()

    def aoLoopClicked(self,pressed):
        if(pressed):
#            self.aoLoopThread = aoLoopThread(self.mainLoop, self.config, self.img, self.ui.wao_strehlSE, self.ui.wao_strehlLE, 0, self.imgType, self.numberSelected)
#            self.connect(self.aoLoopThread,QtCore.SIGNAL("finished()"),self.aoLoopFinished)
            self.aoLoopThread.go = True
            self.aoLoopThread.framebyframe=0
            self.ui.wao_strehlSE.clear()
            self.ui.wao_strehlLE.clear()
            self.aoLoopThread.start()
        else:
            self.aoLoopThread.go = False

    def aoLoopFinished(self):
        if not self.ui.wao_run.isChecked:
            self.ui.wao_run.click()

    def aoLoopOpen(self,pressed):
        if(pressed):
            self.rtc.set_openloop(0,1)
        else:
            self.rtc.set_openloop(0,0)

    def loadConfig(self):
        self.configpath = self.defaultParPath+str(self.ui.wao_selectConfig.currentText())
        sys.path.insert(0, self.defaultParPath)
        exec("import %s as config" % str(wao.ui.wao_selectConfig.currentText()).split(".py")[0])
        self.config = config
        self.aoLoopThread.config = config
        self.updateNumberSelector(textType=self.imgType)
        self.updatePanels()

    def setNumberSelection(self):
        if(self.ui.wao_selectNumber.currentIndex() > -1):
            self.numberSelected = self.ui.wao_selectNumber.currentIndex()
        else:
            self.numberSelected = 0
        self.aoLoopThread.numberSelected = self.numberSelected
        self.updateDisplay()

    def updateNumberSelector(self,textType=None):
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
        if(self.aoLoopThread):
            self.aoLoopThread.imgType = self.imgType

    def loadDefaultConfig(self):
        parlist = sorted(glob.glob(self.defaultParPath+"*.py"))
        self.ui.wao_selectConfig.clear()
        self.ui.wao_selectConfig.addItems([parlist[i].split('/')[-1] for i in range(len(parlist))])

    def updateDisplay(self):

        data = None
        if not self.displayLock.acquire(False):
            #print " Display locked"
            return
        else:
            try:
                if(self.ui.wao_Display.isChecked()):
                    self.ui.wao_rtcWindowMPL.hide()
                    self.ui.wao_pgwindow.show()




                    if(self.atm):
                        if(self.SRCrossX and (self.imgType in ["Phase - Target","Phase - DM","Phase - Atmos","Phase - WFS","Spots - WFS","Centroids - WFS","Slopes - WFS"])):
                            self.SRCrossX.hide()
                            self.SRCrossY.hide()

                        if(self.imgType == "Phase - Atmos"):
                            data = self.atm.get_screen(self.config.p_atmos.alt[self.numberSelected])
                            if(self.imgType!=self.currentViewSelected):
                                self.p1.setRange(xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType

                    if(self.wfs):
                        if(self.imgType == "Phase - WFS"):
                            data =self.wfs.get_phase(self.numberSelected)
                            if(self.imgType!=self.currentViewSelected):
                                self.p1.setRange(xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType

                        if(self.imgType == "Spots - WFS"):
                            if(self.config.p_wfss[self.numberSelected].type_wfs == "sh"):
                                data =self.wfs.get_binimg(self.numberSelected)
                            elif(self.config.p_wfss[self.numberSelected].type_wfs == "pyr"):
                                data =self.wfs.get_pyrimg(self.numberSelected)
                            if(self.imgType!=self.currentViewSelected):
                                self.p1.setRange(xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType

                        if(self.imgType == "Centroids - WFS"):
                            if(not self.ui.wao_rtcWindowMPL.isVisible()):
                                self.ui.wao_rtcWindowMPL.show()
                                self.ui.wao_pgwindow.hide()
                            self.ui.wao_rtcWindowMPL.canvas.axes.clear()
                            centroids = self.rtc.getCentroids(0) # retrieving centroids
                            nvalid = [2*o._nvalid for o in self.config.p_wfss]
                            ind = np.sum(nvalid[:self.numberSelected])
                            x, y, vx, vy = tools.plsh(centroids[ind:ind+nvalid[self.numberSelected]], self.config.p_wfss[self.numberSelected].nxsub,self.config.p_tel.cobs , returnquiver=True) # Preparing mesh and vector for display
                            self.ui.wao_rtcWindowMPL.canvas.axes.quiver(x, y, vx, vy, pivot='mid')
                            self.ui.wao_rtcWindowMPL.canvas.draw()
                            self.currentViewSelected = self.imgType

                            return
                        if(self.imgType == "Slopes - WFS"):
                            if(not self.ui.wao_rtcWindowMPL.isVisible()):
                                self.ui.wao_rtcWindowMPL.show()
                                self.ui.wao_pgwindow.hide()
                            self.ui.wao_rtcWindowMPL.canvas.axes.clear()
                            self.wfs.slopes_geom(self.numberSelected,0)
                            slopes = self.wfs.get_slopes(self.numberSelected)
                            x, y, vx, vy = tools.plsh(slopes, self.config.p_wfss[self.numberSelected].nxsub,self.config.p_tel.cobs , returnquiver=True) # Preparing mesh and vector for display
                            wao.ui.wao_rtcWindowMPL.canvas.axes.quiver(x, y, vx, vy, pivot='mid')
                            wao.ui.wao_rtcWindowMPL.canvas.draw()
                            self.currentViewSelected = self.imgType

                            return

                    if(self.dms):
                        if(self.imgType == "Phase - DM"):
                            dm_type = self.config.p_dms[self.numberSelected].type_dm
                            alt = self.config.p_dms[self.numberSelected].alt
                            data =self.dms.get_dm(dm_type,alt)
                            if(self.imgType!=self.currentViewSelected):
                                self.p1.setRange(xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType
                    if(self.tar):
                        if(self.imgType == "Phase - Target"):
                            data =self.tar.get_phase(self.numberSelected)
                            if(self.imgType!=self.currentViewSelected):
                                self.p1.setRange(xRange=(0, data.shape[0]), yRange=(0, data.shape[1]))
                            self.currentViewSelected = self.imgType

                        if(self.imgType == "PSF SE"):
                            data =self.tar.get_image(self.numberSelected,"se")
                            if(self.ui.wao_PSFlogscale.isChecked()):
                                data = np.log10(data)

                            if (not self.SRCrossX):
                                Delta = 5
                                self.SRCrossX = pg.PlotCurveItem(np.array([data.shape[0]/2+0.5-Delta, data.shape[0]/2+0.5+Delta]), np.array([data.shape[1]/2+0.5, data.shape[1]/2+0.5]), pen='r')
                                self.SRCrossY = pg.PlotCurveItem(np.array([data.shape[0]/2+0.5, data.shape[0]/2+0.5]), np.array([data.shape[1]/2+0.5-Delta, data.shape[1]/2+0.5+Delta]), pen='r')
                                self.p1.addItem(self.SRCrossX) # Put image in plot area
                                self.p1.addItem(self.SRCrossY) # Put image in plot area

                            if(self.imgType!=self.currentViewSelected):
                                zoom = 50
                                self.SRCrossX.show()
                                self.SRCrossY.show()
                                self.p1.setRange(xRange=(data.shape[0]/2+0.5-zoom, data.shape[0]/2+0.5+zoom), yRange=(data.shape[1]/2+0.5-zoom, data.shape[1]/2+0.5+zoom),)
                            self.currentViewSelected = self.imgType




                        if(self.imgType == "PSF LE"):
                            data =self.tar.get_image(self.numberSelected,"le")
                            if(self.ui.wao_PSFlogscale.isChecked()):
                                data = np.log10(data)
                            if (not self.SRCrossX):
                                Delta = 5
                                self.SRCrossX = pg.PlotCurveItem(np.array([data.shape[0]/2+0.5-Delta, data.shape[0]/2+0.5+Delta]), np.array([data.shape[1]/2+0.5, data.shape[1]/2+0.5]), pen='r')
                                self.SRCrossY = pg.PlotCurveItem(np.array([data.shape[0]/2+0.5, data.shape[0]/2+0.5]), np.array([data.shape[1]/2+0.5-Delta, data.shape[1]/2+0.5+Delta]), pen='r')

                                self.p1.addItem(self.SRCrossX) # Put image in plot area
                                self.p1.addItem(self.SRCrossY) # Put image in plot area
                            if(self.imgType!=self.currentViewSelected):
                                zoom = 50
                                self.p1.setRange(xRange=(data.shape[0]/2+0.5-zoom, data.shape[0]/2+0.5+zoom), yRange=(data.shape[1]/2+0.5-zoom, data.shape[1]/2+0.5+zoom))
                                self.SRCrossX.show()
                                self.SRCrossY.show()

                            self.currentViewSelected = self.imgType



                    if (data is not None):
                        autoscale = self.ui.wao_autoscale.isChecked()
                        if(autoscale):
                            self.hist.setLevels(data.min(), data.max()) # inits levels
                        self.img.setImage(data, autoLevels=autoscale)
                        #self.p1.autoRange()



            finally:
                self.displayLock.release()

    def InitConfig(self):
        self.currentViewSelected = None
        self.SRCrossX = None
        self.SRCrossY = None

        self.manually_destroy()
        #set simulation name
        if(hasattr(self.config,"simul_name")):
            if(self.config.simul_name is None):
                simul_name=""
            else:
                simul_name=self.config.simul_name
        else:
            simul_name=""
        matricesToLoad={}
        if(simul_name=="" or not self.ui.wao_useDatabase.isChecked()):
            clean=1
        else:
            clean=0
            param_dict = h5u.params_dictionary(self.config)
            matricesToLoad = h5u.checkMatricesDataBase(os.environ["SHESHA_ROOT"]+"/data/",self.config,param_dict)

        self.c = ch.naga_context(self.ui.wao_deviceNumber.value())
        self.wfs,self.tel=ao.wfs_init(self.config.p_wfss,self.config.p_atmos,self.config.p_tel,
                             self.config.p_geom,self.config.p_target,self.config.p_loop,
                             1,0,self.config.p_dms)

        self.atm=ao.atmos_init(self.c, self.config.p_atmos, self.config.p_tel,
                                                self.config.p_geom,self.config.p_loop,
                                                self.config.p_wfss,self.config.p_target,clean=clean,load=matricesToLoad)
        self.ui.wao_atmosDimScreen.setText(str(self.config.p_atmos.dim_screens[0]))

        self.dms=ao.dm_init(self.config.p_dms,self.config.p_wfss,self.config.p_geom,self.config.p_tel)

        self.tar=ao.target_init(self.c, self.tel, self.config.p_target,self.config.p_atmos,
                                                  self.config.p_geom,self.config.p_tel,self.config.p_wfss,
                                                  self.wfs,self.config.p_dms)

        self.rtc=ao.rtc_init(self.tel,self.wfs,self.config.p_wfss,self.dms,self.config.p_dms,
                             self.config.p_geom,self.config.p_rtc,self.config.p_atmos,
                             self.atm,self.config.p_tel,self.config.p_loop,self.tar,
                             self.config.p_target,clean=clean,simul_name=simul_name, load=matricesToLoad)

        if(not clean):
            h5u.validDataBase(os.environ["SHESHA_ROOT"]+"/data/",matricesToLoad)
        self.mainLoop = [self.tel,self.atm,self.wfs,self.rtc,self.tar,self.dms]
        self.aoLoopThread.wfs = self.wfs
        self.aoLoopThread.atm = self.atm
        self.aoLoopThread.tar = self.tar
        self.aoLoopThread.rtc = self.rtc
        self.aoLoopThread.dms = self.dms
        self.aoLoopThread.tel = self.tel
        self.updateDisplay()
        self.displayRtcMatrix()
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

    def resetDM(self):
        if(self.dms):
            ndm = self.ui.wao_selectDM.currentIndex()
            if(ndm > -1):
                self.dms.resetdm(str(self.ui.wao_dmTypeSelector.currentText()),self.ui.wao_dmAlt.value())
                self.updateDisplay()
            else:
                print "Invalid DM : please select a DM to reset"
        else:
            print "There is not any dm to reset"

    def displayRtcMatrix(self):
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
                self.rtc.set_cmm(0,tmp)
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
                    self.ui.wao_rtcWindow.canvas.axes.matshow(data, aspect="auto", origin="lower", cmap="gray")
                elif(len(data.shape) == 1):
                    self.ui.wao_rtcWindow.canvas.axes.plot(range(len(data)),data, color="black") # TODO : plot it properly, interactivity ?
                    ax.set_yscale('log')
                    if(type_matrix == "Eigenvalues"):
#                        major_ticks = np.arange(0, 101, 20)
#                        minor_ticks = np.arange(0, 101, 5)
#
#                        self.ui.wao_rtcWindow.canvas.axes.set_xticks(major_ticks)
#                        self.ui.wao_rtcWindow.canvas.axes.set_xticks(minor_ticks, minor=True)
#                        self.ui.wao_rtcWindow.canvas.axes.set_yticks(major_ticks)
#                        self.ui.wao_rtcWindow.canvas.axes.set_yticks(minor_ticks, minor=True)
#
#                        # and a corresponding grid

                        self.ui.wao_rtcWindow.canvas.axes.grid(which='both')

                        # or if you want differnet settings for the grids:
                        self.ui.wao_rtcWindow.canvas.axes.grid(which='minor', alpha=0.2)
                        self.ui.wao_rtcWindow.canvas.axes.grid(which='major', alpha=0.5)
                        nfilt = self.rtc.get_nfiltered(0, self.config.p_rtc)
                        self.ui.wao_rtcWindow.canvas.axes.plot([nfilt-0.5,nfilt-0.5],[data.min(), data.max()], color="red", linestyle="dashed")
                        if(nfilt>0):
                            self.ui.wao_rtcWindow.canvas.axes.scatter(np.arange(0,nfilt,1),data[0:nfilt], color="red") # TODO : plot it properly, interactivity ?
                            self.ui.wao_rtcWindow.canvas.axes.scatter(np.arange(nfilt,len(data),1),data[nfilt:], color="blue") # TODO : plot it properly, interactivity ?
                            tt = "%d modes Filtered" % nfilt
                            #ax.text(nfilt + 2, data[nfilt-1], tt)
                            ax.text(0.5, 0.2, tt, horizontalalignment='center',
     verticalalignment='center', transform = ax.transAxes)

                self.ui.wao_rtcWindow.canvas.draw()




class aoLoopThread(QtCore.QThread):
    def __init__(self,LoopParams,config, img, strehlSE, strehlLE, framebyframe, histo, RTDisplay, RTDFreq, updateDisplay, imgType=None, numberSelected=None):
        QtCore.QThread.__init__(self)

        self.tel = LoopParams[0]
        self.wfs = LoopParams[2]
        self.atm = LoopParams[1]
        self.tar = LoopParams[4]
        self.rtc = LoopParams[3]
        self.dms = LoopParams[5]
        self.config = config
        self.go = True
        self.img = img
        self.imgType = imgType
        self.numberSelected = numberSelected
        self.framebyframe = framebyframe
        self.strehlSE = strehlSE
        self.strehlLE = strehlLE
        self.histo = histo
        self.RTDFreq = RTDFreq
        self.RTDisplay = RTDisplay
        self.updateDisplay = updateDisplay

    def __del__(self):
        self.wait()

    """
    def updateDisplay(self):
        data = None
        if((not wao.ui.wao_pgwindow.isVisible()) & (self.imgType != "Centroids - WFS") & (self.imgType != "Slopes - WFS")):
            wao.ui.wao_pgwindow.show()
            wao.ui.wao_rtcWindowMPL.hide()

        if(self.atm):
            if(self.imgType == "Phase - Atmos"):
                data = self.atm.get_screen(self.config.p_atmos.alt[self.numberSelected])
        if(self.wfs):
            if(self.imgType == "Phase - WFS"):
                data =self.wfs.get_phase(self.numberSelected)
            if(self.imgType == "Spots - WFS"):
                if(self.config.p_wfss[self.numberSelected].type_wfs == "sh"):
                    data =self.wfs.get_binimg(self.numberSelected)
                elif(self.config.p_wfss[self.numberSelected].type_wfs == "pyr"):
                    data =self.wfs.get_pyrimg(self.numberSelected)
            if(self.imgType == "Centroids - WFS"):
                if(not wao.ui.wao_rtcWindowMPL.isVisible()):
                    wao.ui.wao_rtcWindowMPL.show()
                    wao.ui.wao_pgwindow.hide()
                wao.ui.wao_rtcWindowMPL.canvas.axes.clear()
                centroids = self.rtc.getCentroids(0) # retrieving centroids
                nvalid = [2*o._nvalid for o in self.config.p_wfss]
                ind = np.sum(nvalid[:self.numberSelected])
                x, y, vx, vy = tools.plsh(centroids[ind:ind+nvalid[self.numberSelected]], self.config.p_wfss[self.numberSelected].nxsub,self.config.p_tel.cobs , returnquiver=True) # Preparing mesh and vector for display
                wao.ui.wao_rtcWindowMPL.canvas.axes.quiver(x, y, vx, vy, pivot='mid')
                wao.ui.wao_rtcWindowMPL.canvas.draw()
            if(self.imgType == "Slopes - WFS"):
                if(not wao.ui.wao_rtcWindowMPL.isVisible()):
                    wao.ui.wao_rtcWindowMPL.show()
                    wao.ui.wao_pgwindow.hide()
                wao.ui.wao_rtcWindowMPL.canvas.axes.clear()
                self.wfs.slopes_geom(self.numberSelected,0)
                slopes = self.wfs.get_slopes(self.numberSelected)
                x, y, vx, vy = tools.plsh(slopes, self.config.p_wfss[self.numberSelected].nxsub,self.config.p_tel.cobs , returnquiver=True) # Preparing mesh and vector for display
                wao.ui.wao_rtcWindowMPL.canvas.axes.quiver(x, y, vx, vy, pivot='mid')
                wao.ui.wao_rtcWindowMPL.canvas.draw()

        if(self.dms):
            if(self.imgType == "Phase - DM"):
                dm_type = self.config.p_dms[self.numberSelected].type_dm
                alt = self.config.p_dms[self.numberSelected].alt
                data =self.dms.get_dm(dm_type,alt)
        if(self.tar):
            if(self.imgType == "Phase - Target"):
                data =self.tar.get_phase(self.numberSelected)
            if(self.imgType == "PSF SE"):
                data =self.tar.get_image(self.numberSelected,"se")
                if(wao.ui.wao_PSFlogscale.isChecked()):
                    data = np.log10(data)
            if(self.imgType == "PSF LE"):
                data =self.tar.get_image(self.numberSelected,"le")
                if(wao.ui.wao_PSFlogscale.isChecked()):
                    data = np.log10(data)
        if (data is not None):
            autoscale = wao.ui.wao_autoscale.isChecked()
            if(autoscale):
                wao.hist.setLevels(data.min(), data.max()) # inits levels
            wao.img.setImage(data, autoLevels=autoscale)
            # wao.p1.autoRange()
    """
    def run(self):
        i=0
        if(self.framebyframe):
            self.mainLoop(1)
        else:
            print "Starting loop"
            while(self.go):
                i +=1
                #print i
                self.mainLoop(i)

            print "Loop stopped"

    def mainLoop(self, itnum):
        start = time.time()
        if(self.atm):
                #print itnum
                self.printInPlace("")

                self.atm.move_atmos()

                self.printInPlace("")

        if(self.config.p_controllers[0].type_control == "geo"):
            if(self.tar):
                for t in range(self.config.p_target.ntargets):
                    self.printInPlace("")
                    self.tar.atmos_trace(t,self.atm,self.tel)
                    self.rtc.docontrol_geo(0, self.dms, self.tar, 0)
                    self.rtc.applycontrol(0,self.dms)
                    self.printInPlace("")
                    self.tar.dmtrace(0,self.dms)
        else:
            if(self.tar):
                for t in range(self.config.p_target.ntargets):
                    self.printInPlace("")
                    self.tar.atmos_trace(t,self.atm,self.tel)
                    self.tar.dmtrace(t,self.dms)
            if(self.wfs):
                for w in range(len(self.config.p_wfss)):
                    self.printInPlace("")
                    self.wfs.sensors_trace(w,"all",self.tel,self.atm,self.dms)
                    self.printInPlace("")
                    self.wfs.sensors_compimg(w)
            if(self.rtc):
                self.printInPlace("")
                self.rtc.docentroids(0)
                self.rtc.docontrol(0)
                self.rtc.applycontrol(0,self.dms)
                self.printInPlace("")
        if(self.tar):
            signal_le = ""
            signal_se= ""
            for t in range(self.config.p_target.ntargets):
                SR = self.tar.get_strehl(t)
                signal_se += "%1.2f   "%SR[0]
                signal_le += "%1.2f   "%SR[1]

        if(self.RTDisplay):
            self.updateDisplay()# Update GUI plots
            t = 1/float(self.RTDFreq) # Limit loop frequency
            time.sleep(t)# Limit loop frequency
        #self.loopFreq.setValue(CurrentFreq) #
        CurrentFreq = 1/(time.time() - start)
        if(self.RTDisplay):
            self.emit(QtCore.SIGNAL('currentSRSE(QString)'), signal_se)#"%1.2f"%SR[0])#str(SR[0]))
            self.emit(QtCore.SIGNAL('currentSRLE(QString)'), signal_le)#"%1.2f"%SR[1])#str(SR[1]))
            self.emit(QtCore.SIGNAL('currentLoopFrequency(float)'), CurrentFreq)

        #print CurrentFreq
        self.printInPlace("iter #%d SR: (L.E, S.E.)= %s, %srunning at %4.1fHz" % (itnum, signal_le, signal_se, CurrentFreq))
        #sys.stdout.flush()

    def printInPlace(self, text):
        print "\r" + text ,;# This seems to trigger the GUI and keep it responsive
        #sys.stdout.flush()
        #sys.stdout.write(text)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    wao = widgetAOWindow()
    wao.show()
    #app.connect(wao.ui._quit,QtCore.SIGNAL("clicked()"),app,QtCore.SLOT("quit()"))
    app.setStyle('cleanlooks')
