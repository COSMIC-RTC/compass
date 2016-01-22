
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
WindowTemplate,TemplateBaseClass=loadUiType(os.environ["SHESHA_ROOT"]+"/widgets/widget_canapass.ui")


"""
low levels debugs:
gdb --args python -i widget_ao.py

"""
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
        self.c = ch.naga_context()
        self.tel = None
        self.wfs = None
        self.rtc = None
        self.atm = None
        self.tar = None
        self.dms = None
        self.aoLoopThread = None
        self.mainLoop = [self.tel,self.atm,self.wfs,self.rtc,self.tar,self.dms]
        self.brama_flag=1;

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
        self.ui.wao_loadConfig.clicked.connect(self.addConfigFromFile)
        self.ui.wao_init.clicked.connect(self.InitConfig)
        self.ui.wao_run.setCheckable(True)
        self.ui.wao_run.clicked[bool].connect(self.aoLoopClicked)
        self.ui.wao_next.clicked.connect(self.nextClicked)
        self.imgType = str(self.ui.wao_selectScreen.currentText())
        self.ui.wao_unzoom.clicked.connect(self.p1.autoRange)
        self.ui.wao_selectScreen.currentIndexChanged.connect(partial(self.updateNumberSelector,textType=None))
        self.ui.wao_selectNumber.currentIndexChanged.connect(self.setNumberSelection)
        self.ui.wao_Display.clicked.connect(self.updateFrameRate)
        self.ui.wao_frameRate.valueChanged.connect(self.updateFrameRate)
        self.ui.wao_rtcWindowMPL.hide()
        self.RTDisplay= self.ui.wao_Display.isChecked()
        self.RTDFreq = self.ui.wao_frameRate.value()
        self.ui.wao_PSFlogscale.clicked.connect(self.updateDisplay)
        self.ui.wao_setActiveDevice.clicked.connect(self.setDevice)
        
        
        # Create Loop thread
        self.aoLoopThread = aoLoopThread(self.mainLoop, self.config, self.img, self.ui.wao_strehlSE, self.ui.wao_strehlLE, 1, self.hist, self.RTDisplay, self.RTDFreq)
        self.connect(self.aoLoopThread, QtCore.SIGNAL("currentLoopFrequency(float)"), self.updateCurrentLoopFrequency)
        self.connect(self.aoLoopThread, QtCore.SIGNAL("currentSRSE(QString)"), self.updateSRSE)
        self.connect(self.aoLoopThread, QtCore.SIGNAL("currentSRLE(QString)"), self.updateSRLE)
        
        self.connect(self.aoLoopThread,QtCore.SIGNAL("finished()"),self.aoLoopFinished)
        
    def setDevice(self):
        self.c.set_activeDevice(self.ui.wao_deviceNumber.value())
        
    def manually_destroy(self):
        if(self.atm):
            del self.atm
            del self.aoLoopThread.atm           
        if(self.tel):
            del self.tel
            del self.aoLoopThread.tel
        if(self.wfs):
            del self.wfs
            del self.aoLoopThread.wfs
        if(self.rtc):
            del self.rtc
            del self.aoLoopThread.rtc
        if(self.tar):
            del self.tar
            del self.aoLoopThread.tar
        if(self.dms):
            del self.dms
            del self.aoLoopThread.dms
            
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
            self.destroy()
            
            #quit()
            #sys.exit(app.exec_())
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
            
              
    def addConfigFromFile(self):
        filepath = str(QtGui.QFileDialog(directory=self.defaultParPath).getOpenFileName(self,"Select parameter file","","parameters file (*.py);;hdf5 file (*.h5);;all files (*)"))
        self.configpath = filepath
        filename = filepath.split('/')[-1]
        if(filepath.split('.')[-1] == "py"):
            pathfile = filepath.split(filename)[-1]    
            #if (pathfile not in sys.path):
            sys.path.insert(0, pathfile)      
            exec("import %s as config" % filename.split(".py")[0])
            sys.path.remove(pathfile)
        elif(filepath.split('.')[-1] == "h5"):
            sys.path.insert(0,self.defaultParPath)
            import scao_16x16_8pix as config
            sys.path.remove(self.defaultParPath)
            h5u.configFromH5(filepath,config)
        else:
            raise ValueError("Parameter file extension must be .py or .h5")
        self.config = config
        self.aoLoopThread.config = config
        self.updateNumberSelector(textType=self.imgType)
        
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
        
    def updateDisplay(self):
        data = None
        if(wao.ui.wao_Display.isChecked()):
            self.ui.wao_rtcWindowMPL.hide()
            self.ui.wao_pgwindow.show()
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
                    return
                    
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
                    if(self.ui.wao_PSFlogscale.isChecked()):
                        data = np.log10(data)
                if(self.imgType == "PSF LE"):
                    data =self.tar.get_image(self.numberSelected,"le")
                    if(self.ui.wao_PSFlogscale.isChecked()):
                        data = np.log10(data)
            if (data is not None):
                self.img.setImage(data)
                self.hist.setLevels(np.min(data),np.max(data))
                self.p1.autoRange()
                
            
    def InitConfig(self):
        
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
        if(simul_name==""):
            clean=1
        else:
            clean=0
            param_dict = h5u.params_dictionary(self.config)
            matricesToLoad = h5u.checkMatricesDataBase(os.environ["SHESHA_ROOT"]+"/data/",self.config,param_dict)
        self.wfs,self.tel=ao.wfs_init(self.config.p_wfss,self.config.p_atmos,self.config.p_tel,
                             self.config.p_geom,self.config.p_target,self.config.p_loop,
                             1,0,self.config.p_dms)

        self.atm=ao.atmos_init(self.c, self.config.p_atmos, self.config.p_tel,
                                                self.config.p_geom,self.config.p_loop,
                                                self.config.p_wfss,self.config.p_target,clean=clean,load=matricesToLoad)

        self.dms=ao.dm_init(self.config.p_dms,self.config.p_wfss,self.config.p_geom,self.config.p_tel)

        self.tar=ao.target_init(self.c, self.tel, self.config.p_target,self.config.p_atmos,
                                                  self.config.p_geom,self.config.p_tel,self.config.p_wfss,
                                                  self.wfs,self.config.p_dms, brama=self.brama_flag)

        self.rtc=ao.rtc_init(self.tel,self.wfs,self.config.p_wfss,self.dms,self.config.p_dms,
                             self.config.p_geom,self.config.p_rtc,self.config.p_atmos,
                             self.atm,self.config.p_tel,self.config.p_loop,self.tar,
                             self.config.p_target,clean=clean,simul_name=simul_name, load=matricesToLoad, brama=self.brama_flag)
                             
        if(simul_name is not ""):
            h5u.validDataBase(os.environ["SHESHA_ROOT"]+"/data/",matricesToLoad)
        self.mainLoop = [self.tel,self.atm,self.wfs,self.rtc,self.tar,self.dms]
        self.aoLoopThread.wfs = self.wfs
        self.aoLoopThread.atm = self.atm
        self.aoLoopThread.tar = self.tar
        self.aoLoopThread.rtc = self.rtc
        self.aoLoopThread.dms = self.dms
        self.aoLoopThread.tel = self.tel
        self.updateDisplay()
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
            
class aoLoopThread(QtCore.QThread):
    def __init__(self,LoopParams,config, img, strehlSE, strehlLE, framebyframe, histo, RTDisplay, RTDFreq, imgType=None, numberSelected=None):
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

    def __del__(self):
        self.wait()
    
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
            self.img.setImage(data, autoLevels=False)
            #self.histo.setLevels(np.min(data),np.max(data))

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
                for t in range(self.config.p_target.ntargets):
                    self.printInPlace("") 
                    self.tar.atmos_trace(t,self.atm,self.tel)
                    self.tar.dmtrace(t,self.dms)
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

        if(wao.brama_flag):
            self.rtc.publish() #rtc_publish, g_rtc;
            self.tar.publish()
            

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
