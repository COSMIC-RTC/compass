## @package   shesha.widgets.widget_canapass
## @brief     Widget to simulate a closed loop using CANAPASS
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   4.3.0
## @date      2011/01/28
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
#  All rights reserved.
#  Distributed under GNU - LGPL
#
#  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
#  General Public License as published by the Free Software Foundation, either version 3 of the License,
#  or any later version.
#
#  COMPASS: End-to-end AO simulation tool using GPU acceleration
#  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
#
#  The final product includes a software package for simulating all the critical subcomponents of AO,
#  particularly in the context of the ELT and a real-time core based on several control approaches,
#  with performances consistent with its integration into an instrument. Taking advantage of the specific
#  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
#  conduct large simulation campaigns called to the ELT.
#
#  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
#  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
#  various systems configurations such as multi-conjugate AO.
#
#  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
#  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
#  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
"""
Widget to simulate a closed loop using CANAPASS

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

from shesha.widgets.widget_ao import widgetAOWindow, widgetAOWindow

global server
server = None


class widgetCanapassWindowPyro(widgetAOWindow):

    def __init__(self, configFile: Any = None, cacao: bool = False,
                 expert: bool = False) -> None:
        widgetAOWindow.__init__(self, configFile, cacao, hideHistograms=True)
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
        widgetAOWindow.initConfig(self)
        global server
        server = self.startPyroServer()

    def loadConfig(self, *args, configFile=None, supervisor=None) -> None:
        '''
            Callback when 'LOAD' button is hit
        '''
        if supervisor is None:
            supervisor = CanapassSupervisor(configFile)

        widgetAOWindow.loadConfig(self, args, configFile=configFile,
                                  supervisor=supervisor)

    def loopOnce(self) -> None:
        widgetAOWindow.loopOnce(self)
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
    wao = widgetCanapassWindowPyro(arguments["<parameters_filename>"], cacao=True)
    wao.show()
    # if arguments["--interactive"]:
    #     from shesha.util.ipython_embed import embed
    #     from os.path import basename
    #     embed(basename(__file__), locals())
