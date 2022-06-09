## @package   naga.context
## @brief     Documentation for naga
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   5.3.0
## @date      2022/01/24
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

from carmaWrap import context
import numpy as np


class Context:
    """
    Python class for wrapping a CarmaContext
    """

    def __init__(self, devices=None):
        """ Initializes a Context object

    Args:
        devices: (list): (optionnal) list of devices to use. Default is 0
        """
        if devices is None:
            self.context = context.get_instance_1gpu(0)
            self.devices = 0
        else:
            if isinstance(devices, list):
                devices = np.array(devices, dtype=np.int32)
            else:
                raise TypeError("Devices must be a list of integer")
            self.devices = devices
            self.context = context.get_instance_ngpu(devices)

    def getActiveDevice(self):
        """ Return the index of the current active device """
        return self.context.active_device

    def setActiveDevice(self, index: int):
        """ Set the device index as the active device """
        if index in self.devices:
            self.context.set_active_device(index)
        else:
            raise ValueError("Index given is not valid")

    def enableTensorCores(self):
        """ Enable the tensor cores math mode """
        self.context.activate_tensor_cores(True)

    def disableTensorCores(self):
        """ Disable the tensor cores math mode """
        self.context.activate_tensor_cores(True)
