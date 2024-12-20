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

"""
Widget to use on a bench

Usage:
  widget_bench.py [<parameters_filename>] [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -d, --devices devices      Specify the devices
  -i, --interactive  keep the script interactive
"""

import os
import sys
import time

try:
    from PyQt6 import QtWidgets
    from PyQt6.QtCore import Qt
except ModuleNotFoundError:
    try:
        from PySide2 import QtWidgets
        from PySide2.QtCore import Qt
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "No module named 'PyQt6' or PySide2', please install one of them\nException raised: "
            + e.msg
        )

from typing import Any

from docopt import docopt

from shesha.widgets.widget_base import WidgetBase, uiLoader

from shesha.supervisor.benchSupervisor import (
    WFSType,
    BenchSupervisor as Supervisor,
)

BenchWindowTemplate, BenchClassTemplate = uiLoader("widget_bench")

# For debug
# from IPython.core.debugger import Pdb
# then add this line to create a breakpoint
# Pdb().set_trace()


class widgetBenchWindow(BenchClassTemplate, WidgetBase):
    def __init__(self, config_file: Any = None, devices: str = None) -> None:
        WidgetBase.__init__(self)
        BenchClassTemplate.__init__(self)

        self.devices = devices

        self.uiBench = BenchWindowTemplate()
        self.uiBench.setupUi(self)

        #                   ATTRIBUTES                              #

        self.supervisor = None
        self.config = None
        self.stop = False  # type: bool  # Request quit

        self.uiBench.wao_nbiters.setValue(1000)  # Default GUI nIter box value
        self.nbiter = self.uiBench.wao_nbiters.value()
        self.refreshTime = 0  # type: float  # System time at last display refresh
        self.loopThread = None  # type: QThread
        self.assistant = None  # type: Any

        #                 CONNECTED BUTTONS                         #
        # Default path for config files
        self.defaultParPath = os.environ["SHESHA_ROOT"] + "/data/par/bench-config"
        self.defaultAreaPath = os.environ["SHESHA_ROOT"] + "/data/layouts"
        self.loadDefaultConfig()

        self.uiBench.wao_run.setCheckable(True)
        self.uiBench.wao_run.clicked[bool].connect(self.aoLoopClicked)
        self.uiBench.wao_open_loop.setCheckable(True)
        self.uiBench.wao_open_loop.clicked[bool].connect(self.aoLoopOpen)
        self.uiBench.wao_next.clicked.connect(self.loop_once)

        self.uiBench.wao_forever.stateChanged.connect(self.updateForever)

        self.dispStatsInTerminal = False

        self.uiBench.wao_run.setDisabled(True)
        self.uiBench.wao_next.setDisabled(True)
        self.uiBench.wao_unzoom.setDisabled(True)

        self.addDockWidget(Qt.DockWidgetArea(1), self.uiBase.wao_ConfigDock)
        self.addDockWidget(Qt.DockWidgetArea(1), self.uiBase.wao_DisplayDock)
        self.uiBase.wao_ConfigDock.setFloating(False)
        self.uiBase.wao_DisplayDock.setFloating(False)

        self.adjustSize()

        if config_file is not None:
            self.uiBase.wao_selectConfig.clear()
            self.uiBase.wao_selectConfig.addItem(config_file)
            self.load_config()
            self.init_config()

        #                       METHODS                             #

    # def updateStatsInTerminal(self, state):
    #     self.dispStatsInTerminal = state

    def updateForever(self, state):
        self.uiBench.wao_nbiters.setDisabled(state)

    def add_dispDock(self, name: str, parent, type: str = "pg_image") -> None:
        d = WidgetBase.add_dispDock(self, name, parent, type)
        if type == "SR":
            d.addWidget(self.uiBench.wao_Strehl)

    def load_config(self) -> None:
        """
        Callback when 'LOAD' button is hit
        """
        WidgetBase.load_config(self)
        config_file = str(self.uiBase.wao_selectConfig.currentText())
        sys.path.insert(0, self.defaultParPath)

        self.supervisor = Supervisor(config_file)

        self.config = self.supervisor.get_config()

        # if self.devices:
        #     self.config.p_loop.set_devices([
        #             int(device) for device in self.devices.split(",")
        #     ])

        try:
            sys.path.remove(self.defaultParPath)
        except Exception:
            pass

        try:
            self.nwfs = len(self.config.p_wfss)
        except Exception:
            self.nwfs = 1  # Default. config may very well not define p_wfss
        for wfs in range(self.nwfs):
            name = "slpComp_%d" % wfs
            self.add_dispDock(name, self.wao_graphgroup_cb, "MPL")
            name = "wfs_%d" % wfs
            self.add_dispDock(name, self.wao_imagesgroup_cb)

        self.uiBench.wao_run.setDisabled(True)
        self.uiBench.wao_next.setDisabled(True)
        self.uiBench.wao_unzoom.setDisabled(True)

        self.uiBase.wao_init.setDisabled(False)

        if hasattr(self.config, "layout"):
            area_filename = self.defaultAreaPath + "/" + self.config.layout + ".area"
            self.loadArea(filename=area_filename)

        self.adjustSize()

    def aoLoopClicked(self, pressed: bool) -> None:
        if pressed:
            self.stop = False
            self.refreshTime = time.time()
            self.nbiter = self.uiBench.wao_nbiters.value()
            if self.dispStatsInTerminal:
                if self.uiBench.wao_forever.isChecked():
                    print("LOOP STARTED")
                else:
                    print("LOOP STARTED FOR %d iterations" % self.nbiter)
            self.run()
        else:
            self.stop = True

    def aoLoopOpen(self, pressed: bool) -> None:
        if pressed:
            self.supervisor.open_loop()
            self.uiAO.wao_open_loop.setText("Close Loop")
        else:
            self.supervisor.close_loop()
            self.uiAO.wao_open_loop.setText("Open Loop")

    def init_config(self) -> None:
        WidgetBase.init_config(self)

    def init_configThread(self) -> None:
        # self.uiBench.wao_deviceNumber.setDisabled(True)
        # self.config.p_loop.devices = self.uiBench.wao_deviceNumber.value()  # using GUI value
        # gpudevice = "ALL"  # using all GPU avalaible
        # gpudevice = np.array([2, 3], dtype=np.int32)
        # gpudevice = np.arange(4, dtype=np.int32) # using 4 GPUs: 0-3
        # gpudevice = 0  # using 1 GPU : 0
        self.supervisor.init_config()

    def init_configFinished(self) -> None:
        # Thread carma context reload:
        try:
            self.supervisor.force_context()
        except Exception:
            print("Warning: could not call supervisor.force_context().")

        for i in range(self.nwfs):
            if self.config.p_wfss[i].type == WFSType.SH:
                key = "wfs_%d" % i
                self.addSHGrid(
                    self.docks[key].widgets[0],
                    self.config.p_wfss[i].get_validsub(),
                    self.config.p_wfss[i].npix,
                    self.config.p_wfss[i].npix,
                )

        self.updateDisplay()

        self.uiBench.wao_run.setDisabled(False)
        self.uiBench.wao_next.setDisabled(False)
        self.uiBench.wao_open_loop.setDisabled(False)
        self.uiBench.wao_unzoom.setDisabled(False)

        WidgetBase.init_configFinished(self)

    def updateDisplay(self) -> None:
        if (
            (self.supervisor is None)
            or (not self.supervisor.is_init())
            or (not self.uiBase.wao_Display.isChecked())
        ):
            # print("Widget not fully initialized")
            return
        if not self.loopLock.acquire(False):
            return
        else:
            try:
                for key, dock in self.docks.items():
                    if key == "Strehl":
                        continue
                    elif dock.isVisible():
                        index = int(key.split("_")[-1])
                        data = None
                        if "wfs" in key:
                            data = self.supervisor.wfs.get_wfs_image(index)
                            if data is not None:
                                autoscale = True  # self.uiBench.actionAuto_Scale.isChecked()
                                # if (autoscale):
                                #     # inits levels
                                #     self.hist.setLevels(data.min(), data.max())
                                self.imgs[key].setImage(data, autoLevels=autoscale)
                                # self.p1.autoRange()

                        elif "slp" in key:  # Slope display
                            if (
                                self.config.p_wfss[index].type == WFSType.PYRHR
                                or self.config.p_wfss[index].type == WFSType.PYRLR
                            ):
                                raise RuntimeError("PYRHR not usable")
                            self.imgs[key].canvas.axes.clear()
                            x, y = self.supervisor.config.p_wfss[index].get_validsub()

                            nssp = x.size
                            centroids = self.supervisor.rtc.get_slopes()
                            vx = centroids[:nssp]
                            vy = centroids[nssp:]

                            offset = (self.supervisor.config.p_wfss[index].npix - 1) / 2
                            self.imgs[key].canvas.axes.quiver(
                                x + offset,
                                y + offset,
                                vy,
                                vx,
                                angles="xy",
                                scale_units="xy",
                                scale=1,
                            )  # wao.supervisor.config.p_wfss[0].pixsize)
                            self.imgs[key].canvas.draw()

            finally:
                self.loopLock.release()

    def loop_once(self) -> None:
        if not self.loopLock.acquire(False):
            print("Display locked")
            return
        else:
            try:
                start = time.time()
                self.supervisor.single_next()
                loopTime = time.time() - start
                refreshDisplayTime = 1.0 / self.uiBase.wao_frameRate.value()

                if time.time() - self.refreshTime > refreshDisplayTime:
                    currentFreq = 1 / loopTime
                    self.uiBench.wao_currentFreq.setValue(currentFreq)
                    self.refreshTime = start
            except Exception:
                pass
            finally:
                self.loopLock.release()

    def run(self):
        WidgetBase.run(self)
        if not self.uiBench.wao_forever.isChecked():
            self.nbiter -= 1

        if self.nbiter <= 0:
            self.stop = True
            self.uiBench.wao_run.setChecked(False)


if __name__ == "__main__":
    arguments = docopt(__doc__)
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("cleanlooks")
    wao = widgetBenchWindow(arguments["<parameters_filename>"], devices=arguments["--devices"])
    wao.show()
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename

        embed(basename(__file__), locals())
