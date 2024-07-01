## @package   shesha.config.pHrtc
## @brief     ParamHrtc class definition
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


#################################################
# P-Class (parametres) ParamHrtc
#################################################
class ParamHrtc:

    def __init__(self):
        self.__wfs_payload_size = 1000 
        """ MUDPI payload size for WFS frame"""
        self.__com_payload_size = 1000 
        """ MUDPI payload size for DM command"""
        self.__hrtc_host = None 
        """ H-RTC host in format 'IP:port'"""
        self.__local_host = None 
        """ local host in format 'IP:port'"""
        self.__framesize = None 
        """ H-RTC WFS frame size"""
        self.__frame_shm = None 
        """ H-RTC WFS frame SHM name for localhost link"""
        self.__com_shm = None 
        """ H-RTC DM commands SHM name for localhost link"""

    def get_hrtc_host(self):
        """ Get the hrtc host

        :return: (str) : hrtc host
        """
        return self.__hrtc_host

    def set_hrtc_host(self, host):
        """ Set the hrtc host

        Args:
            t: (str) : hrtc host
        """
        self.__hrtc_host = host

    hrtc_host = property(get_hrtc_host, set_hrtc_host)

    def get_local_host(self):
        """ Get the local host

        :return: (str) : local host
        """
        return self.__local_host

    def set_local_host(self, host):
        """ Set the local host

        Args:
            t: (str) : local host
        """
        self.__local_host = host

    local_host = property(get_local_host, set_local_host)

    def get_wfs_payload_size(self):
        """ Get the payload size for WFS frame

        :return: (int) : WFS payload size
        """
        return self.__wfs_payload_size

    def set_wfs_payload_size(self, ps):
        """ set the payload size for WFS frame

        :return: (int) : WFS payload size
        """
        self.__wfs_payload_size = csu.enforce_int(ps)

    wfs_payload_size = property(get_wfs_payload_size, set_wfs_payload_size)

    def get_com_payload_size(self):
        """ Get the payload size for DM com

        :return: (int) : com payload size
        """
        return self.__com_payload_size

    def set_com_payload_size(self, ps):
        """ set the payload size for DM com

        :return: (int) : com payload size
        """
        self.__com_payload_size = csu.enforce_int(ps)

    com_payload_size = property(get_com_payload_size, set_com_payload_size)    
    
    def get_framesize(self):
        """ Get the hrtc wfs frame size

        :return: (int) : com payload size
        """
        return self.__framesize

    def set_framesize(self, fs):
        """ set the hrtc wfs frame size

        :return: (int) : com payload size
        """
        self.__framesize = csu.enforce_int(fs)

    framesize = property(get_framesize, set_framesize)
    
    def get_frame_shm(self):
        """ Get the hrtc WFS Frame SHM name

        :return: (str) : hrtc WFS Frame SHM name
        """
        return self.__frame_shm

    def set_frame_shm(self, name):
        """ Set the hrtc WFS Frame SHM name

        Args:
            t: (str) : hrtc WFS Frame SHM name
        """
        self.__frame_shm = name

    frame_shm = property(get_frame_shm, set_frame_shm)

    def get_com_shm(self):
        """ Get the hrtc DM commands SHM name

        :return: (str) : hrtc DM commands SHM name
        """
        return self.__com_shm

    def set_com_shm(self, name):
        """ Set the hrtc DM commands SHM name

        Args:
            t: (str) : hrtc DM commands SHM name
        """
        self.__com_shm = name

    com_shm = property(get_com_shm, set_com_shm)
