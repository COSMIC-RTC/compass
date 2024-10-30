#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
#
# COMPASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COMPASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COMPASS. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2024 COSMIC Team
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import (
    NavigationToolbar2QT as Navigationtoolbar,
)
from matplotlib.figure import Figure
from matplotlib import gridspec

try:
    from PyQt6 import QtWidgets
except ModuleNotFoundError:
    try:
        from PySide2 import QtWidgets
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "No module named 'PyQt6' or PySide2', please install one of them\nException raised: "
            + e.msg
        )

# import matplotlib
# matplotlib.use('Qt5Agg')

# matplotlib.rcParams['backend.qt4']='PySide'
# matplotlib.style.use('ggplot')
# matplotlib.style.use('seaborn-muted')

"""
Created on Tue Jun 24 00:27:01 2014

@author: fvidal
"""
"""
from PyQt4 import QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

from matplotlib.figure import Figure

class MplCanvas(FigureCanvas):

    def __init__(self):
        self.fig = Figure(facecolor='white')
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Policy.Expanding,QtGui.QSizePolicy.Policy.Expanding)
        FigureCanvas.updateGeometry(self)


class MatplotlibWidget(QtGui.QWidget):

    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)
"""


# Embeddable matplotlib figure/canvas
class MplCanvas(FigureCanvas):
    def __init__(self):
        self.fig = Figure(frameon=True)
        self.gs1 = gridspec.GridSpec(1, 1)
        self.axes = self.fig.add_subplot(self.gs1[0], aspect="auto")

        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(
            self,
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding,
        )
        FigureCanvas.updateGeometry(self)


# creates embeddable matplotlib figure/canvas with toolbar
class MatplotlibWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.create_framentoolbar()

    def create_framentoolbar(self):
        self.frame = QtWidgets.QWidget()
        self.canvas = MplCanvas()
        self.canvas.setParent(self.frame)
        self.mpltoolbar = Navigationtoolbar(self.canvas, self.frame)
        self.vbl = QtWidgets.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.vbl.addWidget(self.mpltoolbar)
        self.setLayout(self.vbl)
