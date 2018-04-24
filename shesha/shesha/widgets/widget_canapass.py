"""Widget to simulate a closed loop

Usage:
  widget_canapass.py [<parameters_filename>]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -i, --interactive  keep the script interactive
"""

import os, sys
import numpy as np
import time

import pyqtgraph as pg
from shesha.util.tools import plsh, plpyr
from tqdm import trange
import astropy.io.fits as pfits
from PyQt5 import QtWidgets
from shesha.supervisor.canapassSupervisor import CanapassSupervisor

from typing import Any, Dict, Tuple, Callable, List
from docopt import docopt

from .widget_ao import widgetAOWindow

global server
server = None


class widgetCanapassWindowPyro(widgetAOWindow):

    def __init__(self, configFile: Any=None, BRAMA: bool=False,
                 expert: bool=False) -> None:
        widgetAOWindow.__init__(self, configFile, BRAMA, hideHistograms=True)
        #Pyro.core.ObjBase.__init__(self)

        self.CB = {}
        self.wpyr = None
        self.currentBuffer = 1
        #############################################################
        #                 CONNECTED BUTTONS                         #
        #############################################################
        # Default path for config files
        #self.uiAO.wao_openLoop.setChecked(False)
        #self.uiAO.wao_openLoop.setText("Close Loop")
        self.uiAO.actionShow_Pyramid_Tools.toggled.connect(self.showPyrTools)
        self.wpyrNbBuffer = 1
        #############################################################
        #                       METHODS                             #
        #############################################################

    def initConfig(self) -> None:
        self.supervisor.clearInitSim()
        WidgetBase.initConfig(self)
        global server
        server = self.startPyroServer()

    def loadConfig(self) -> None:
        '''
            Callback when 'LOAD' button is hit
        '''
        WidgetBase.loadConfig(self, ISupervisor=CanapassSupervisor)

    def loopOnce(self) -> None:
        WidgetBase.loopOnce(self)
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

    def startPyroServer(self):
        try:
            from subprocess import Popen, PIPE
            from hraa.server.pyroServer import PyroServer

            # Init looper
            wao_loop = loopHandler(self)

            # Find username
            p = Popen("whoami", shell=True, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            if (err != b''):
                print(err)
                raise ValueError("ERROR CANNOT RECOGNIZE USER")
            else:
                user = out.split(b"\n")[0].decode("utf-8")
                print("User is " + user)

            server = PyroServer()
            server.add_device(self.supervisor, "waoconfig_" + user)
            server.add_device(wao_loop, "waoloop_" + user)
            server.start()
        except:
            raise EnvironmentError("Missing dependencies (code HRAA, Pyro4)")

        return server


class loopHandler:

    def __init__(self, wao):
        self.wao = wao

    def start(self):
        self.wao.aoLoopClicked(True)
        self.wao.uiAO.wao_run.setChecked(True)

    def stop(self):
        self.wao.aoLoopClicked(False)
        self.wao.uiAO.wao_run.setChecked(False)

    def alive(self):
        return "alive"


if __name__ == '__main__':
    arguments = docopt(__doc__)
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('cleanlooks')
    wao = widgetCanapassWindowPyro(arguments["<parameters_filename>"], BRAMA=True)
    wao.show()
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())
