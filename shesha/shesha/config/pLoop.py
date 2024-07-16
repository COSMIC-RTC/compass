## @package   shesha.config.pLoop
## @brief     Class that contains the configuration for the loop module of the COMPASS simulator.
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
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


import numpy as np
import shesha.config.config_setter_utils as csu


class ParamLoop:
    """Class that contains the configuration for the loop module of the COMPASS simulator.

    This class contains the configuration for the loop module of the COMPASS simulator.

    Attributes:
        niter (long): Number of iteration.
        ittime (float): Iteration time.
        devices (np.ndarray[ndim=1, dtype=np.int32_t]): List of GPU devices.
    """

    def __init__(self):
        self.__niter = 0
        self.__ittime = 0.0
        self.__devices = np.array([0], dtype=np.int32)

    def get_devices(self):
        """Get the list of GPU devices used

        :return: (np.ndarray[ndim=1, dtype=np.int32_t]) : list of GPU devices
        """
        return self.__devices

    def set_devices(self, devices):
        """Set the list of GPU devices used

        Args:
            devices: (np.ndarray[ndim=1, dtype=np.int32_t]) : list of GPU devices
        """
        self.__devices = csu.enforce_array(
            devices, len(devices), dtype=np.int32, scalar_expand=False
        )

    devices = property(get_devices, set_devices)

    def get_niter(self):
        """Get the number of iteration

        :return: (long) : number of iteration
        """
        return self.__niter

    def set_niter(self, n):
        """Set the number of iteration

        Args:
            n: (long) : number of iteration
        """
        self.__niter = csu.enforce_int(n)

    niter = property(get_niter, set_niter)

    def get_ittime(self):
        """Get iteration time

        :return: (float) :iteration time
        """
        return self.__ittime

    def set_ittime(self, t):
        """Set iteration time

        Args:
            t: (float) :iteration time
        """
        self.__ittime = csu.enforce_float(t)

    ittime = property(get_ittime, set_ittime)
