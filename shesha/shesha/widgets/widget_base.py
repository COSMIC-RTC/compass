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


import os
import threading
import warnings
from typing import Any, Callable

import numpy as np
import pyqtgraph as pg

try:
    from PyQt6 import QtWidgets, QtCore
    from PyQt6.QtCore import QThread, QTimer
    from PyQt6.uic import loadUiType
except ModuleNotFoundError:
    try:
        from PySide2 import QtWidgets
        from PySide2.QtCore import QThread, QTimer
        from PySide2.QtUiTools import loadUiType
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "No module named 'PyQt6' or PySide2', please install one of them\nException raised: "
            + e.msg
        )

from pyqtgraph.dockarea import Dock, DockArea

from shesha.util.matplotlibwidget import MatplotlibWidget


def uiLoader(moduleName):
    return loadUiType(os.path.dirname(__file__) + f"/{moduleName}.ui")  # type: type, type


BaseWidgetTemplate, BaseClassTemplate = uiLoader("widget_base")


class PupilBoxes(QtWidgets.QGraphicsPathItem):
    def __init__(self, x, y):
        """x and y are 2D arrays of shape (Nplots, Nsamples)"""
        connect = np.ones(x.shape, dtype=bool)
        connect[:, -1] = 0  # don't draw the segment between each trace
        self.path = pg.arrayToQPath(x.flatten(), y.flatten(), connect.flatten())
        pg.QtGui.QGraphicsPathItem.__init__(self, self.path)
        self.setPen(pg.mkPen("r"))

    def shape(
        self,
    ):  # override because QGraphicsPathItem.shape is too expensive.
        return pg.QtGui.QGraphicsItem.shape(self)

    def boundingRect(self):
        return self.path.boundingRect()


class WidgetBase(BaseClassTemplate):
    def __init__(self, parent=None, hide_histograms=False) -> None:
        self.app = QtCore.QCoreApplication.instance()
        BaseClassTemplate.__init__(self, parent=parent)
        self.uiBase = BaseWidgetTemplate()
        print("UI SETUP")
        self.uiBase.setupUi(self)
        #                   ATTRIBUTES                              #
        print("QTimer")

        self.gui_timer = QTimer()  # type: QTimer
        self.gui_timer.timeout.connect(self.updateDisplay)
        if self.uiBase.wao_Display.isChecked():
            self.gui_timer.start(1000 // int(self.uiBase.wao_frameRate.value()))
        self.loopLock = threading.Lock()  # type: Threading.Lock # Asynchronous loop / display safe-threading
        self.hide_histograms = hide_histograms
        #               PYQTGRAPH DockArea INIT                     #
        print("DockArea")

        self.area = DockArea()
        self.uiBase.wao_DisplayDock.setWidget(self.area)
        self.gridSH = []

        #                 CONNECTED BUTTONS                         #
        # Default path for config files
        print("Buttons")
        self.defaultParPath = "."
        self.defaultAreaPath = "."
        self.uiBase.wao_load_config.clicked.connect(self.load_config)
        self.uiBase.wao_loadArea.clicked.connect(self.loadArea)
        self.uiBase.wao_saveArea.clicked.connect(self.saveArea)
        self.uiBase.wao_init.clicked.connect(self.init_config)
        self.uiBase.wao_configFromFile.clicked.connect(self.addConfigFromFile)

        self.uiBase.wao_Display.stateChanged.connect(self.gui_timer_config)
        self.uiBase.wao_frameRate.setValue(2)

        self.uiBase.wao_load_config.setDisabled(False)
        self.uiBase.wao_init.setDisabled(True)

        self.disp_checkboxes = []
        self.docks = {}  # type: Dict[str, pg.dockarea.Dock]
        self.viewboxes = {}  # type: Dict[str, pg.ViewBox]
        self.imgs = {}  # type: Dict[str, pg.ImageItem]
        self.hists = {}  # type: Dict[str, pg.HistogramLUTItem]

        self.PupilLines = None
        self.adjustSize()

    def gui_timer_config(self, state) -> None:
        self.uiBase.wao_frameRate.setDisabled(state)
        if state:
            self.gui_timer.start(1000 // int(self.uiBase.wao_frameRate.value()))
        else:
            self.gui_timer.stop()

    def closeEvent(self, event: Any) -> None:
        self.quitGUI(event)

    def quitGUI(self, event: Any = None) -> None:
        reply = QtWidgets.QMessageBox.question(
            self,
            "Message",
            "Are you sure to quit?",
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No,
        )

        if reply == QtWidgets.QMessageBox.StandardButton.Yes:
            if event:
                event.accept()
            quit()
        else:
            if event:
                event.ignore()

        #                       METHODS                             #

    def saveArea(self, widget, filename=None):
        """
        Callback when a area layout file is double clicked in the file browser
        Place the selected file name in the browsing drop-down menu,
        the call the self.load_config callback of the load button.
        """
        if filename is None:
            filepath = QtWidgets.QFileDialog(directory=self.defaultAreaPath).getSaveFileName(
                self,
                "Select area layout file",
                "",
                "area layout file (*.area);;all files (*)",
            )
            filename = filepath[0]

        try:
            with open(filename, "w+") as f:
                st = self.area.saveState()
                f.write(str(st))
        except FileNotFoundError as err:
            warnings.warn(filename + " not loaded: " + err)

    def showDock(self, name):
        for disp_checkbox in self.disp_checkboxes:
            if disp_checkbox.text() == name:
                disp_checkbox.setChecked(True)
                break
        if name in self.docks.keys():
            self.area.addDock(self.docks[name])

    def restoreMyState(self, state):
        typ, contents, _ = state

        if typ == "dock":
            self.showDock(contents)
        else:
            for o in contents:
                self.restoreMyState(o)

    def loadArea(self, widget=None, filename=None):
        # close all docks
        for disp_checkbox in self.disp_checkboxes:
            disp_checkbox.setChecked(False)
        for dock in self.docks.values():
            if dock.isVisible():
                dock.close()

        if filename is None:
            filepath = QtWidgets.QFileDialog(directory=self.defaultAreaPath).getOpenFileName(
                self,
                "Select area layout file",
                "",
                "area layout file (*.area);;all files (*)",
            )
            filename = filepath[0]

        try:
            with open(filename, "r") as f:
                st = eval(f.readline())

                # restore docks from main area
                if st["main"] is not None:
                    self.restoreMyState(st["main"])

                # restore docks from floating area
                for win in st["float"]:
                    self.restoreMyState(win[0]["main"])

                # rearange dock s as in stored state
                self.area.restoreState(st)
        except FileNotFoundError as err:
            warnings.warn(filename + "not loaded: " + err)

    def addConfigFromFile(self) -> None:
        """
        Callback when a config file is double clicked in the file browser
        Place the selected file name in the browsing drop-down menu,
        the call the self.load_config callback of the load button.
        """
        filepath = QtWidgets.QFileDialog(directory=self.defaultParPath).getOpenFileName(
            self,
            "Select parameter file",
            "",
            "parameters file (*.py);;hdf5 file (*.h5);;all files (*)",
        )

        self.uiBase.wao_selectConfig.clear()
        self.uiBase.wao_selectConfig.addItem(str(filepath[0]))

        self.load_config(config_file=self.uiBase.wao_selectConfig.currentText())

    def update_displayDock(self):
        guilty_guy = self.sender().text()
        state = self.sender().isChecked()
        if state:
            self.area.addDock(self.docks[guilty_guy])
        elif self.docks[guilty_guy].isVisible():
            self.docks[guilty_guy].close()

    def add_dispDock(self, name: str, parent, type: str = "pg_image") -> Dock:
        checkBox = QtWidgets.QCheckBox(name, parent)
        checkBox.clicked.connect(self.update_displayDock)
        checkableAction = QtWidgets.QWidgetAction(parent)
        checkableAction.setDefaultWidget(checkBox)
        parent.addAction(checkableAction)
        self.disp_checkboxes.append(checkBox)

        d = Dock(name)  # , closable=True)
        self.docks[name] = d
        if type == "pg_image":
            img = pg.ImageItem(border="w", image=np.zeros((2, 2)))
            self.imgs[name] = img

            viewbox = pg.ViewBox()

            viewbox.setAspectLocked(True)
            viewbox.addItem(img)  # Put image in plot area
            self.viewboxes[name] = viewbox
            viewbox.invertY(False)

            iv = pg.ImageView(view=viewbox, imageItem=img)
            try:
                cmap = pg.colormap.get("viridis")  # prepare a viridis color map
                iv.setColorMap(cmap)
            except Exception:
                pass
            if self.hide_histograms:
                iv.ui.histogram.hide()
            iv.ui.histogram.autoHistogramRange()  # init levels
            iv.ui.histogram.setMaximumWidth(100)
            iv.ui.menuBtn.hide()
            iv.ui.roiBtn.hide()
            d.addWidget(iv)
        elif type == "pg_plot":
            img = pg.PlotItem(border="w")
            self.imgs[name] = img
            d.addWidget(img)
        elif type == "MPL":
            img = MatplotlibWidget()
            self.imgs[name] = img
            d.addWidget(img)
        # elif type == "SR":
        #     d.addWidget(self.uiBase.wao_Strehl)
        return d

    def load_config(self, *args, **kwargs) -> None:
        """
        Callback when 'LOAD' button is hit
        """
        for groupbox in [
            self.uiBase.wao_phasesgroup_tb,
            self.uiBase.wao_imagesgroup_tb,
            self.uiBase.wao_graphgroup_tb,
        ]:
            layout = groupbox.menu()
            while layout and not layout.isEmpty():
                w = layout.children()[0]
                layout.removeAction(w)
                w.setParent(None)
        self.disp_checkboxes.clear()

        # TODO: remove self.imgs, self.viewboxes and self.docks children
        for _, dock in self.docks.items():
            if dock.isVisible():
                dock.close()

        self.docks.clear()
        self.imgs.clear()
        self.viewboxes.clear()

        self.wao_phasesgroup_cb = QtWidgets.QMenu(self)
        self.uiBase.wao_phasesgroup_tb.setMenu(self.wao_phasesgroup_cb)
        self.uiBase.wao_phasesgroup_tb.setText("Select")
        self.uiBase.wao_phasesgroup_tb.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)

        self.wao_graphgroup_cb = QtWidgets.QMenu(self)
        self.uiBase.wao_graphgroup_tb.setMenu(self.wao_graphgroup_cb)
        self.uiBase.wao_graphgroup_tb.setText("Select")
        self.uiBase.wao_graphgroup_tb.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)

        self.uiBase.wao_imagesgroup_tb.setText("Select")
        self.wao_imagesgroup_cb = QtWidgets.QMenu(self)
        self.uiBase.wao_imagesgroup_tb.setMenu(self.wao_imagesgroup_cb)
        self.uiBase.wao_imagesgroup_tb.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)

        # self.uiBase.wao_init.setDisabled(False)
        #
        # if (hasattr(self.sim.config, "layout")):
        #     area_filename = self.defaultAreaPath + "/" + self.sim.config.layout + ".area"
        #     self.loadArea(filename=area_filename)
        #
        # self.adjustSize()

    def loadDefaultConfig(self) -> None:
        import glob

        parlist = sorted(glob.glob(self.defaultParPath + "/*.py"))
        self.uiBase.wao_selectConfig.clear()
        self.uiBase.wao_selectConfig.addItems(
            [parlist[i].split("/")[-1] for i in range(len(parlist))]
        )

    def init_config(self) -> None:
        self.loopLock.acquire(True)
        self.uiBase.wao_load_config.setDisabled(True)
        self.uiBase.wao_init.setDisabled(True)
        self.thread = WorkerThread(self.init_configThread)
        self.thread.finished.connect(self.init_configFinished)
        self.thread.start()

    def init_configThread(self) -> None:
        pass

    def init_configFinished(self) -> None:
        self.uiBase.wao_load_config.setDisabled(False)
        self.uiBase.wao_init.setDisabled(False)
        self.loopLock.release()

    def updateDisplay(self) -> None:
        if not self.loopLock.acquire(False):
            return
        else:
            try:
                pass
            finally:
                self.loopLock.release()

    def addSHGrid(self, pg_image, valid_sub, sspsize, pitch):
        # First remove the old grid, if any
        if self.PupilLines is not None:
            pg_image.removeItem(self.PupilLines)

        nssp_tot = valid_sub[0].size
        connect = np.ones((nssp_tot, 5), dtype=bool)
        connect[:, -1] = 0  # don't draw the segment between each trace
        roi_x = np.ones((nssp_tot, 5), dtype=int)
        roi_y = np.ones((nssp_tot, 5), dtype=int)
        for idx_ssp in range(nssp_tot):
            (x, y) = (valid_sub[0][idx_ssp], valid_sub[1][idx_ssp])
            roi_x[idx_ssp, :] = [x, x, x + sspsize, x + sspsize, x]
            roi_y[idx_ssp, :] = [y, y + sspsize, y + sspsize, y, y]
        self.PupilLines = PupilBoxes(roi_x, roi_y)
        pg_image.addItem(self.PupilLines)

    def printInPlace(self, text: str) -> None:
        # This seems to trigger the GUI and keep it responsive
        print(text, end="\r", flush=True)

    def run(self):
        self.loop_once()
        if not self.stop:
            QTimer.singleShot(0, self.run)  # Update loop


class WorkerThread(QThread):
    def __init__(self, loopFunc: Callable) -> None:
        QThread.__init__(self)
        self.loopFunc = loopFunc

    def run(self) -> None:
        self.loopFunc()

    def stop(self) -> None:
        pass

    def cleanUp(self) -> None:
        pass
