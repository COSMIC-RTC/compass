
import cProfile
import pstats as ps

import sys
sys.path.insert(0, "./lib") # for SESAME lib

import numpy as np
import chakra as ch
import chakra_ao as ao
import time
import matplotlib.pyplot as pl
import os
import pyqtgraph as pg
import glob
sys.path.insert(0, os.environ["CHAKRA_AO"]+"/data/par/") # for SESAME lib

from PyQt4.uic import loadUiType
from PyQt4 import QtCore, QtGui

from functools import partial

WindowTemplate,TemplateBaseClass=loadUiType("widget_ao.ui")



class widgetAOWindow(TemplateBaseClass):  
    def __init__(self):
        TemplateBaseClass.__init__(self)
        
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
        ##############################################################
        #######       ATTRIBUTES   #######################
        #############################################################
        self.configpath = None
        self.config = None
        self.c = ch.chakra_context()
        self.wfs = None
        self.rtc = None
        self.atm = None
        self.tar = None
        self.dms = None
        self.aoLoopThread = None
        self.mainLoop = [self.atm,self.wfs,self.rtc,self.tar,self.dms]
        
        ##############################################################
        #######       PYQTGRAPH WINDOW INIT   #######################
        #############################################################
        
        self.img = pg.ImageItem() # create image area
        self.p1 = self.ui.wao_pgwindow.addPlot() # create pyqtgraph plot area
        self.p1.addItem(self.img) # Put image in plot area
        self.hist = pg.HistogramLUTItem() #Create an histogram
        self.hist.setImageItem(self.img) # Compute histogram from img
        self.ui.wao_pgwindow.addItem(self.hist)
        self.hist.autoHistogramRange() # init levels
        self.hist.setMaximumWidth(100)      
        
        ##############################################################
        #######       CONNECTED BUTTONS  #######################
        #############################################################
        self.defaultParPath = os.environ["CHAKRA_AO"]+"/data/par/par4bench/" # Default path for config files
        self.ui.wao_loadConfig.clicked.connect(self.loadConfig)
        self.loadDefaultConfig()
        self.ui.wao_init.clicked.connect(self.InitConfig)
        self.ui.wao_run.setCheckable(True)
        self.ui.wao_run.clicked[bool].connect(self.aoLoopClicked)
        self.ui.wao_next.clicked.connect(self.nextClicked)
        self.imgType = str(self.ui.wao_selectScreen.currentText())
        self.ui.wao_configFromFile.clicked.connect(self.addConfigFromFile)
        self.ui.wao_selectScreen.currentIndexChanged.connect(partial(self.updateNumberSelector,textType=None))
        self.ui.wao_selectNumber.currentIndexChanged.connect(self.setNumberSelection)
        self.ui.wao_selectAtmosLayer.currentIndexChanged.connect(self.setLayerSelection)
        self.ui.wao_setAtmos.clicked.connect(self.setAtmosParams)
 
        self.aoLoopThread = aoLoopThread(self.mainLoop, self.config, self.img, self.ui.wao_strehlSE, self.ui.wao_strehlLE, 1)
        self.connect(self.aoLoopThread,QtCore.SIGNAL("finished()"),self.aoLoopFinished)

        ##############################################################
        ##################       METHODS      #######################
        #############################################################
        
    def updateTelescopePanel(self):
        self.ui.wao_zenithAngle.setValue(self.config.p_geom.zenithangle)
        self.ui.wao_diamTel.setValue(self.config.p_tel.diam)
        self.ui.wao_cobs.setValue(self.config.p_tel.cobs)
    
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
        if(self.config.p_atmos.dim_screens):
            self.ui.wao_atmosDimScreen.setText(str(self.config.p_atmos.dim_screens[nscreen]))
        
    def updatePanels(self):
        self.updateTelescopePanel()
        self.updateAtmosPanel()
        self.updateLayerSelection()
    
    def setAtmosParams(self):
        nscreen = self.ui.wao_selectAtmosLayer.currentIndex()
        if(nscreen < 0):
            nscreen = 0
        self.config.p_atmos.alt[nscreen]=self.ui.wao_atmosAlt.value()
        self.config.p_atmos.frac[nscreen]=self.ui.wao_atmosFrac.value()
        self.config.p_atmos.L0[nscreen]=self.ui.wao_atmosL0.value()
        self.config.p_atmos.windspeed[nscreen]=self.ui.wao_windSpeed.value()
        self.config.p_atmos.winddir[nscreen]=self.ui.wao_windDirection.value()
        
    def updateLayerSelection(self):
        self.ui.wao_selectAtmosLayer.addItems([str(i) for i in range(self.config.p_atmos.nscreens)])
        
    def setLayerSelection(self):
        self.updateAtmosPanel()
        
    def addConfigFromFile(self):
        filepath = str(QtGui.QFileDialog.getOpenFileName(self,"Select parameter file","","parameters file (*.py);;all files (*)"))
        self.configpath = filepath
        filename = filepath.split('/')[-1]
        self.ui.wao_selectConfig.addItem(filename,0)
        pathfile = filepath.split(filename)[-1]    
        #if (pathfile not in sys.path):
        sys.path.insert(0, pathfile)      
        exec("import %s as config" % filename.split(".py")[0])
        sys.path.remove(pathfile)
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
        if(textType == "Phase - WFS" or textType == "Spots - WFS" or textType == "Slopes - WFS"):
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
        if(self.atm):
            if(self.imgType == "Phase - Atmos"):
                self.img.setImage(self.atm.get_screen(self.config.p_atmos.alt[self.numberSelected]))
        if(self.wfs):
            if(self.imgType == "Phase - WFS"):
                self.img.setImage(self.wfs.get_phase(self.numberSelected))
            if(self.imgType == "Spots - WFS"):
                if(self.config.p_wfss[self.numberSelected].type_wfs == "sh"):
                    self.img.setImage(self.wfs.get_binimg(self.numberSelected))
                elif(self.config.p_wfss[self.numberSelected].type_wfs == "pyr"):
                    self.img.setImage(self.wfs.get_pyrimg(self.numberSelected))
            if(self.imgType == "Slopes - WFS"):
                pass
        if(self.dms):
            if(self.imgType == "Phase - DM"):
                dm_type = self.config.p_dms[self.numberSelected].type_dm
                alt = self.config.p_dms[self.numberSelected].alt
                self.img.setImage(self.dms.get_dm(dm_type,alt))
        if(self.tar):
            if(self.imgType == "Phase - Target"):
                self.img.setImage(self.tar.get_phase(self.numberSelected))
            if(self.imgType == "PSF SE"):
                self.img.setImage(self.tar.get_image(self.numberSelected,"se"))
            if(self.imgType == "PSF LE"):
                self.img.setImage(self.tar.get_image(self.numberSelected,"le"))
                
            
    def InitConfig(self):
        self.wfs=ao.wfs_init(self.config.p_wfss,self.config.p_atmos,self.config.p_tel,
                             self.config.p_geom,self.config.p_target,self.config.p_loop,
                             1,0,self.config.p_dms)

        self.atm=self.config.p_atmos.atmos_init(self.c,self.config.p_tel,
                                                self.config.p_geom,self.config.p_loop)
        self.ui.wao_atmosDimScreen.setText(str(self.config.p_atmos.dim_screens[0]))

        self.dms=ao.dm_init(self.config.p_dms,self.config.p_wfs0,self.config.p_geom,self.config.p_tel)

        self.tar=self.config.p_target.target_init(self.c,self.config.p_atmos,
                                                  self.config.p_geom,self.config.p_tel,self.config.p_wfss,
                                                  self.wfs,self.config.p_dms)

        self.rtc=ao.rtc_init(self.wfs,self.config.p_wfss,self.dms,self.config.p_dms,
                             self.config.p_geom,self.config.p_rtc,self.config.p_atmos,
                             self.atm,self.config.p_tel,self.config.p_loop,self.tar,
                             self.config.p_target)
        self.mainLoop = [self.atm,self.wfs,self.rtc,self.tar,self.dms]
        self.aoLoopThread.wfs = self.wfs
        self.aoLoopThread.atm = self.atm
        self.aoLoopThread.tar = self.tar
        self.aoLoopThread.rtc = self.rtc
        self.aoLoopThread.dms = self.dms
        print "Inits done"


class aoLoopThread(QtCore.QThread):
    def __init__(self,LoopParams,config, img, strehlSE, strehlLE, framebyframe, imgType=None, numberSelected=None):
        QtCore.QThread.__init__(self)
        
        self.wfs = LoopParams[1]
        self.atm = LoopParams[0]
        self.tar = LoopParams[3]
        self.rtc = LoopParams[2]
        self.dms = LoopParams[4]
        self.config = config
        self.go = True
        self.img = img
        self.imgType = imgType
        self.numberSelected = numberSelected
        self.framebyframe = framebyframe
        self.strehlSE = strehlSE
        self.strehlLE = strehlLE
        
    def __del__(self):
        self.wait()
    
    def updateDisplay(self):
        if(self.atm):
            if(self.imgType == "Phase - Atmos"):
                self.img.setImage(self.atm.get_screen(self.config.p_atmos.alt[self.numberSelected]))
        if(self.wfs):
            if(self.imgType == "Phase - WFS"):
                self.img.setImage(self.wfs.get_phase(self.numberSelected))
            if(self.imgType == "Spots - WFS"):
                if(self.config.p_wfss[self.numberSelected].type_wfs == "sh"):
                    self.img.setImage(self.wfs.get_binimg(self.numberSelected))
                elif(self.config.p_wfss[self.numberSelected].type_wfs == "pyr"):
                    self.img.setImage(self.wfs.get_pyrimg(self.numberSelected))
            if(self.imgType == "Slopes - WFS"):
                pass
        if(self.dms):
            if(self.imgType == "Phase - DM"):
                dm_type = self.config.p_dms[self.numberSelected].type_dm
                alt = self.config.p_dms[self.numberSelected].alt
                self.img.setImage(self.dms.get_dm(dm_type,alt))
        if(self.tar):
            if(self.imgType == "Phase - Target"):
                self.img.setImage(self.tar.get_phase(self.numberSelected))
            if(self.imgType == "PSF SE"):
                self.img.setImage(self.tar.get_image(self.numberSelected,"se"))
            if(self.imgType == "PSF LE"):
                self.img.setImage(self.tar.get_image(self.numberSelected,"le"))

    def run(self):
        i=0
        if(self.framebyframe):
            self.mainLoop()
        else:
            print "Starting loop"
            while(self.go):
                i +=1
                self.mainLoop()
#            if(self.imgType == "Phase Atmos"):
#                self.img.setImage(self.atm.get_screen())
#            self.img.setImage(self.tar.get_image(0,"se"))
            
            #print "SR : ",self.tar.get_strehl(0)[0]
            print "Loop stopped"
    
    def mainLoop(self):
        if(self.atm):
                self.atm.move_atmos()
        if(self.config.p_controllers[0].type_control == "geo"):
            if(self.tar):
                for t in range(self.config.p_target.ntargets):
                    self.tar.atmos_trace(t,self.atm)
                    self.rtc.docontrol_geo(0, self.dms, self.tar, 0)
                    self.rtc.applycontrol(0,self.dms)
                    self.tar.dmtrace(0,self.dms)
        else:
            if(self.tar):
                for t in range(self.config.p_target.ntargets):
                    self.tar.atmos_trace(t,self.atm)
                    self.tar.dmtrace(t,self.dms)
            if(self.wfs):
                for w in range(len(self.config.p_wfss)):
                    self.wfs.sensors_trace(w,"all",self.atm,self.dms)
                    self.wfs.sensors_compimg(w)
            if(self.rtc):
                self.rtc.docentroids(0)
                self.rtc.docontrol(0)
                self.rtc.applycontrol(0,self.dms)
        
        if(self.tar):
            SR = self.tar.get_strehl(0)
            self.strehlSE.setText("%.2f"%(SR[0]))
            self.strehlLE.setText("%.2f"%(SR[1]))
        
        self.updateDisplay()
        time.sleep(0.01)
        
        a1 = pg.ArrowItem(angle=-160, tipAngle=60, headLen=40, tailLen=40, tailWidth=20, pen={'color': 'w', 'width': 3})

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    wao = widgetAOWindow()
    wao.show()
