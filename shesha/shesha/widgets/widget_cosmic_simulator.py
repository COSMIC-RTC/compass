## @package   shesha.widgets.widget_cosmic
## @brief     Widget to simulate a closed loop using cosmic
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @version   4.3.0
## @date      2022/01/24
## @copyright 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>
#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>

# COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation, either version 3 of the 
# License, or any later version.

# COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
# If not, see <https://www.gnu.org/licenses/>

# Copyright (C) 2011-2024 COSMIC Team <https//://github.com/COSMIC-RTC/compass>

"""
Widget to simulate a closed loop using cosmic

Usage:
  widget_cosmic_simulator.py [<parameters_filename>]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -i, --interactive  keep the script interactive
"""

import sys


try:
    from PyQt5 import QtWidgets
except ModuleNotFoundError:
    try:    
        from PySide2 import QtWidgets
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError("No module named 'PyQt5' or PySide2', please install one of them\nException raised: "+e.msg)

from shesha.supervisor.cosmicSupervisor import CosmicSupervisor
from typing import Any
from docopt import docopt

from shesha.widgets.widget_base import WidgetBase
from shesha.widgets.widget_ao import widgetAOWindow

global server
server = None


class widgetCosmic(widgetAOWindow):

    def __init__(self, config_file: Any = None, expert: bool = False) -> None:
        widgetAOWindow.__init__(self, config_file, hide_histograms=True)

        #############################################################
        #                       METHODS                             #
        #############################################################
    def init_config(self) -> None:
        self.supervisor = CosmicSupervisor(self.config)
        WidgetBase.init_config(self)


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
    wao = widgetCosmic(arguments["<parameters_filename>"])
    wao.show()
    # if arguments["--interactive"]:
    #     from shesha.util.ipython_embed import embed
    #     from os.path import basename
    #     embed(basename(__file__), locals())
