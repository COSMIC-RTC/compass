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
from tools import plsh, plpyr
from tqdm import trange
import astropy.io.fits as pfits
from PyQt5 import QtWidgets
from supervisor.canapassSupervisor import CanapassSupervisor

from typing import Any, Dict, Tuple, Callable, List
from docopt import docopt

from widget_ao import widgetAOWindow

class widgetCanapassWindowPyro(widgetAOWindow):

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

        self.uiAO.actionShow_Pyramid_Tools.toggled.connect(self.showPyrTools)
        self.wpyrNbBuffer = 1
        #############################################################
        #                       METHODS                             #
        #############################################################

    def initConfig(self) -> None:
        self.supervisor.clearInitSim()
        super().initConfig()

    def loadConfig(self) -> None:
        '''
            Callback when 'LOAD' button is hit
        '''
        super().loadConfig(ISupervisor=CanapassSupervisor)

    def loopOnce(self) -> None:
        super().loopOnce()
        if (self.uiAO.actionShow_Pyramid_Tools.isChecked()):  # PYR only
            self.wpyr.Fe = 1 / self.config.p_loop.ittime  # needs Fe for PSD...
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
        for i in trange(nbIters):
            self.supervisor.singleNext(showAtmos=see_atmos)

    def initPyrTools(self):
        ADOPTPATH = os.getenv("ADOPTPATH")
        sys.path.append(ADOPTPATH + "/widgets")
        from pyrStats import widget_pyrStats
        print("OK Pyramid Tools Widget initialized")
        self.wpyr = widget_pyrStats()
        self.wpyrNbBuffer = self.wpyr.CBNumber
        self.wpyr.show()

    def setPyrToolsParams(self, ai):
        self.wpyr.pup = self.supervisor.getSpupil()
        self.wpyr.phase = self.supervisor.getTargetPhase(0)
        self.wpyr.updateResiduals(ai)
        if (self.supervisor.ph2modes is None):
            print('computing phase 2 Modes basis')
            self.supervisor.computePh2Modes()
        self.wpyr.ph2modes = self.supervisor.ph2modes

    def showPyrTools(self):
        if (self.wpyr is None):
            try:
                print("Lauching pyramid widget...")
                self.initPyrTools()
                print("Done")
            except:
                raise ValueError("ERROR: ADOPT  not found. Cannot launch Pyramid tools")
        else:
            if (self.uiAO.actionShow_Pyramid_Tools.isChecked()):
                self.wpyr.show()
            else:
                self.wpyr.hide()

    def getAi(self):
        return self.wpyr.ai

    def updateSRSE(self, SRSE):
        self.uiAO.wao_strehlSE.setText(SRSE)

    def updateSRLE(self, SRLE):
        self.uiAO.wao_strehlLE.setText(SRLE)

    def updateCurrentLoopFrequency(self, freq):
        self.uiAO.wao_currentFreq.setValue(freq)

if __name__ == '__main__':
    arguments = docopt(__doc__)
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('cleanlooks')
    wao = widgetCanapassWindowPyro(arguments["<parameters_filename>"], BRAMA=True,
                             expert=arguments["--expert"])
    wao.show()
