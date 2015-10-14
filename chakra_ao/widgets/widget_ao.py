
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
        self.c = ch.chakra_context()
        self.wfs = None
        self.rtc = None
        self.atm = None
        self.tar = None
        self.dms = None
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
        
        
        
        ##############################################################
        #######       CONNECTED BUTTONS  #######################
        #############################################################
        self.defaultParPath = os.environ["CHAKRA_AO"]+"/data/par/par4bench/" # Default path for config files
        self.ui.wao_loadConfig.clicked.connect(self.loadConfig)
        self.loadDefaultConfig()
        self.ui.wao_init.clicked.connect(self.InitConfig)
        self.ui.wao_run.setCheckable(True)
        self.ui.wao_run.clicked[bool].connect(self.aoLoopClicked)
        self.imgType = str(self.ui.wao_selectScreen.currentText())
        self.ui.wao_configFromFile.clicked.connect(self.addConfigFromFile)
        self.ui.wao_selectScreen.currentIndexChanged.connect(partial(self.updateNumberSelector,textType=None))
        
        ##############################################################
        ##################       METHODS      #######################
        #############################################################
        
    def addConfigFromFile(self):
        filepath = str(QtGui.QFileDialog.getOpenFileName(self,"Select parameter file","","parameters file (*.py);;all files (*)"))
        self.configpath = filepath
        filename = filepath.split('/')[-1]
        self.ui.wao_selectConfig.addItem(filename,0)
        pathfile = filepath.split(filename)[-1]    
        if (pathfile not in sys.path):
            sys.path.insert(0, pathfile)      
        exec("import %s as config" % filename.split(".py")[0])
        sys.path.remove(pathfile)
        self.config = config
        self.updateNumberSelector(textType=self.imgType)
        
    def aoLoopClicked(self,pressed):
        if(pressed):
            self.aoLoopThread = aoLoopThread(self.mainLoop, self.config, self.img, self.imgType)
            self.connect(self.aoLoopThread,QtCore.SIGNAL("finished()"),self.aoLoopFinished)
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
        self.updateNumberSelector(textType=self.imgType)
        
    def updateNumberSelector(self,textType=None):
        if(textType == None):
            print "coucou"
            textType = str(self.ui.wao_selectScreen.currentText())
        self.imgType = textType
        self.ui.wao_selectNumber.clear()
        print textType
        if(textType == "Phase - Atmos"):
            n = self.config.p_atmos.nscreens
        if(textType == "Phase - WFS"):
            n = len(self.config.p_wfss)
        self.ui.wao_selectNumber.addItems([str(i) for i in range(n)])
        
    def loadDefaultConfig(self):
        parlist = sorted(glob.glob(self.defaultParPath+"*.py"))
        self.ui.wao_selectConfig.clear()
        self.ui.wao_selectConfig.addItems([parlist[i].split('/')[-1] for i in range(len(parlist))])
    
    def InitConfig(self):
        self.wfs=ao.wfs_init(self.config.p_wfss,self.config.p_atmos,self.config.p_tel,
                             self.config.p_geom,self.config.p_target,self.config.p_loop,
                             1,0,self.config.p_dms)

        self.atm=self.config.p_atmos.atmos_init(self.c,self.config.p_tel,
                                                self.config.p_geom,self.config.p_loop)

        self.dms=ao.dm_init(self.config.p_dms,self.config.p_wfs0,self.config.p_geom,self.config.p_tel)

        self.tar=self.config.p_target.target_init(self.c,self.config.p_atmos,
                                                  self.config.p_geom,self.config.p_tel,self.config.p_wfss,
                                                  self.wfs,self.config.p_dms)

        self.rtc=ao.rtc_init(self.wfs,self.config.p_wfss,self.dms,self.config.p_dms,
                             self.config.p_geom,self.config.p_rtc,self.config.p_atmos,
                             self.atm,self.config.p_tel,self.config.p_loop,self.tar,
                             self.config.p_target)
        self.mainLoop = [self.atm,self.wfs,self.rtc,self.tar,self.dms]
        print "Inits done"


class aoLoopThread(QtCore.QThread):
    def __init__(self,LoopParams,config, img, imgType):
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
        print "thread created"
        
    def __del__(self):
        self.wait()
        
    def run(self):
        i=0
        print "Starting loop"
        while(self.go):
            i +=1
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
            
            if(self.imgType == "Phase Atmos"):
                self.img.setImage(self.atm.get_screen())
            self.img.setImage(self.tar.get_image(0,"se"))
            
            #print "SR : ",self.tar.get_strehl(0)[0]
        print "Loop stopped"

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    wao = widgetAOWindow()
    wao.show()
#
#if(len(sys.argv)!=2):
#    error= 'command line should be:"python -i test.py parameters_filename"\n with "parameters_filename" the path to the parameters file'
#    raise StandardError(error)
#
##get parameters from file
#param_file=sys.argv[1]
#execfile(param_file)
#start=param_file.rindex("/")
#end=param_file.rindex(".")
##simul_name=param_file[start+1:end]
#print "param_file is",param_file
##print "simul name is",simul_name
#
#
#
##initialisation:
##   context
#c=ch.chakra_context()
#c.set_activeDevice(0)
#
##    wfs
#print "->wfs"
#wfs=ao.wfs_init(p_wfss,p_atmos,p_tel,p_geom,p_target,p_loop, 1,0,p_dms)
#
##   atmos
#print "->atmos"
#atm=p_atmos.atmos_init(c,p_tel,p_geom,p_loop)
#
##   dm 
#print "->dm"
#dms=ao.dm_init(p_dms,p_wfs0,p_geom,p_tel)
#
##   target
#print "->target"
#tar=p_target.target_init(c,p_atmos,p_geom,p_tel,p_wfss,wfs,p_dms)
#
#print "->rtc"
##   rtc
#rtc=ao.rtc_init(wfs,p_wfss,dms,p_dms,p_geom,p_rtc,p_atmos,atm,p_tel,p_loop,tar,p_target)#,simul_name=simul_name)
#
#print "===================="
#print "init done"
#print "===================="
#print "objects initialzed on GPU:"
#print "--------------------------------------------------------"
#print atm
#print wfs
#print dms
#print tar
#print rtc
#
#mimg = 0.# initializing average image
#print "----------------------------------------------------";
#print "iter# | S.E. SR | L.E. SR | Est. Rem. | framerate";
#print "----------------------------------------------------";
#
#def loop( n):
#    #fig,((turbu,image),(shak,defMir))=pl.subplots(2,2, figsize=(15,15))
#    #pl.ion()
#    #pl.show()
#    for i in range(n):
#        atm.move_atmos()
#        #print "atmos"
#        if(p_controller0.type_control == "geo"):
#            for t in range(p_target.ntargets):
#                tar.atmos_trace(t,atm)
#                rtc.docontrol_geo(0, dms, tar, 0)
#                rtc.applycontrol(0,dms)
#                tar.dmtrace(0,dms)
#        else:
#            for t in range(p_target.ntargets):
#                tar.atmos_trace(t,atm)
#                #print "atmos_trace : ",t
#                tar.dmtrace(t,dms)
#                #print "dm_trace : ",t
#            for w in range(len(p_wfss)):
#                wfs.sensors_trace(w,"all",atm,dms)
#                #print "sensors_trace : ",w
#                wfs.sensors_compimg(w)
#                #print "compimg : ",w
#
#            rtc.docentroids(0)
#            #print "do_centroids"
#            rtc.docontrol(0)
#            #print "do_control"
#            rtc.applycontrol(0,dms)
#            #print "apply_control"
#        
#
#        if((i+1)%100==0):
#            """
#            turbu.clear()
#            image.clear()
#            shak.clear()
#            defMir.clear()
#
#            screen=atm.get_screen(p_atmos.alt)
#            f1=turbu.matshow(screen,cmap='Blues_r')
#
#            im=tar.get_image(0,"se")
#            im=np.roll(im,im.shape[0]/2,axis=0)
#            im=np.roll(im,im.shape[1]/2,axis=1)
#            f2=image.matshow(im,cmap='Blues_r')
#
#            sh=wfs.get_binimg(0)
#            f3=shak.matshow(sh,cmap='Blues_r')
#
#            dm=dms.get_dm("pzt",0.)
#            f4=defMir.matshow(dm)
#
#            pl.draw()
#            """
#            strehltmp = tar.get_strehl(0)
#            print i+1,"\t",strehltmp[0],"\t",strehltmp[1]
#
#
#
##loop(p_loop.niter)
