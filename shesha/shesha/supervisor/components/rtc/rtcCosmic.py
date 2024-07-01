## @package   shesha.supervisor
## @brief     User layer for initialization and execution of a COMPASS simulation
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

from shesha.supervisor.components.sourceCompass import SourceCompass
import shesha.constants as scons
import numpy as np
from typing import Union

import tides.streamers as ts

class RtcCosmic():
    """ RTC handler for compass simulation
    """

    def __init__(self, config, wfs, dms):
        """ Initialize a RtcCompass component for rtc related supervision

        Args:

            config : (config module) : Parameters configuration structure module

            wfs: (Sensors) : Sensors object

            dms: (Dms) : Dms object

        """
        self.framesize = config.p_hrtc.framesize
        self._wfs = wfs
        self._dms = dms
        self.config = config

        com_shape = np.sum([p_dm._ntotact for p_dm in config.p_dms])
        self.hrtc_host = config.p_hrtc.hrtc_host
        self.local_host = config.p_hrtc.local_host

        self.frame = ts.MudpiFrame(np.zeros((self.framesize, self.framesize), dtype=np.uint16), config.p_hrtc.wfs_payload_size)
        self.com = ts.MudpiFrame(np.zeros(com_shape, dtype=np.float32), config.p_hrtc.com_payload_size)

        self.publisher = ts.Streamer("rtms")
        self.publisher.configure(self.hrtc_host, self.frame)

        self.subscriber = ts.Streamer("rtms")
        self.subscriber.configure(self.local_host, self.com)

        self.framecounter = 2
        img_shape = wfs.get_wfs_image(0).shape
        self.crop = (img_shape[0] - self.framesize) // 2

    def do_control(self):
        """ Send WFS frame to H-RTC and receipt DM command
        """
        wfs_cropped = self._wfs.get_wfs_image(0)[self.crop:self.crop + self.framesize, self.crop:self.crop + self.framesize]
        self.frame = ts.MudpiFrame(wfs_cropped.astype(np.uint16), self.config.p_hrtc.wfs_payload_size)
        self.frame.framecounter = self.framecounter
        self.publisher.publish(self.hrtc_host, self.frame)
        self.com = self.subscriber.getFrame({self.local_host:1})
        self.framecounter += 1
    
    def apply_control(self):
        """ Apply DM command to DM
        """
        self._dms.set_command(np.array(self.com[self.local_host][0]))