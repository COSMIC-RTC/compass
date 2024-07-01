## @package   shesha.config.pHrtc
## @brief     Class that contains the configuration for the HRTC module of the COMPASS simulator.
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   5.4.4
## @date      2022/01/24
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
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

import numpy as np
from . import config_setter_utils as csu
import shesha.constants as scons


class ParamHrtc:
    """Class that contains the configuration for the HRTC module of the COMPASS simulator.

    This class contains the configuration for the HRTC module of the COMPASS simulator.

    Attributes:
        hrtc_host (str): HRTC host in format 'IP:port'.
        local_host (str): Local host in format 'IP:port'.
        wfs_payload_size (int): MUDPI payload size for WFS frame.
        com_payload_size (int): MUDPI payload size for DM command.
        framesize (int): H-RTC WFS frame size.
        frame_shm (str): H-RTC WFS frame SHM name for localhost link.
        com_shm (str): H-RTC DM commands SHM name for localhost link.
    """

    def __init__(self):
        self.__wfs_payload_size = 1000  # MUDPI payload size for WFS frame
        self.__com_payload_size = 1000  # MUDPI payload size for DM command
        self.__hrtc_host = None  # H-RTC host in format 'IP:port'
        self.__local_host = None  # local host in format 'IP:port'
        self.__framesize = None  # H-RTC WFS frame size
        self.__frame_shm = None  # H-RTC WFS frame SHM name for localhost link
        self.__com_shm = None  # H-RTC DM commands SHM name for localhost link

    @property
    def hrtc_host(self):
        """Get or set the HRTC host.

        Returns:
            str: HRTC host.
        """
        return self.__hrtc_host

    @hrtc_host.setter
    def hrtc_host(self, host):
        self.__hrtc_host = host

    @property
    def local_host(self):
        """Get or set the local host.

        Returns:
            str: Local host.
        """
        return self.__local_host

    @local_host.setter
    def local_host(self, host):
        self.__local_host = host

    @property
    def wfs_payload_size(self):
        """Get or set the payload size for WFS frame.

        Returns:
            int: WFS payload size.
        """
        return self.__wfs_payload_size

    @wfs_payload_size.setter
    def wfs_payload_size(self, ps):
        self.__wfs_payload_size = csu.enforce_int(ps)

    @property
    def com_payload_size(self):
        """Get or set the payload size for DM command.

        Returns:
            int: DM payload size.
        """
        return self.__com_payload_size

    @com_payload_size.setter
    def com_payload_size(self, ps):
        self.__com_payload_size = csu.enforce_int(ps)

    @property
    def framesize(self):
        """Get or set the HRTC WFS frame size.

        Returns:
            int: HRTC WFS frame size.
        """
        return self.__framesize

    @framesize.setter
    def framesize(self, fs):
        self.__framesize = csu.enforce_int(fs)

    @property
    def frame_shm(self):
        """Get or set the HRTC WFS Frame SHM name.

        Returns:
            str: HRTC WFS Frame SHM name.
        """
        return self.__frame_shm

    @frame_shm.setter
    def frame_shm(self, name):
        self.__frame_shm = name

    @property
    def com_shm(self):
        """Get or set the HRTC DM commands SHM name.

        Returns:
            str: HRTC DM commands SHM name.
        """
        return self.__com_shm

    @com_shm.setter
    def com_shm(self, name):
        self.__com_shm = name
