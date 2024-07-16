## @package   shesha.supervisor
## @brief     User layer for initialization and execution of a COMPASS simulation
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @version   5.4.4
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

import CacaoInterfaceWrap as ciw
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
        self.publisher = None
        self.subscriber = None
        self.frame_shm = None
        self.com_shm = None

        if self.config.p_hrtc.frame_shm is None:
            com_shape = np.sum([p_dm._ntotact for p_dm in config.p_dms])
            self.hrtc_host = config.p_hrtc.hrtc_host
            self.local_host = config.p_hrtc.local_host

            self.frame = ts.MudpiFrame(np.zeros((self.framesize, self.framesize), dtype=np.uint16), config.p_hrtc.wfs_payload_size)
            self.com = ts.MudpiFrame(np.zeros(com_shape, dtype=np.float32), config.p_hrtc.com_payload_size)

            self.publisher = ts.Streamer("rtms")
            self.publisher.configure(self.hrtc_host, self.frame)

            self.subscriber = ts.Streamer("rtms")
            self.subscriber.configure(self.local_host, self.com)
        else:
            self.frame_shm = ciw.CacaoInterfaceWrap(self.config.p_hrtc.frame_shm)
            self.com_shm = ciw.CacaoInterfaceWrap(self.config.p_hrtc.com_shm)
            
        img_shape = wfs.get_wfs_image(0).shape
        self.crop = (img_shape[0] - self.framesize) // 2
        self.framecounter = 2
        
    def do_control(self):
        """ Send WFS frame to H-RTC and receipt DM command
        """
        wfs_cropped = self._wfs.get_wfs_image(0)[self.crop:self.crop + self.framesize, self.crop:self.crop + self.framesize]
        if self.publisher is not None:
            self.frame = ts.MudpiFrame(wfs_cropped.astype(np.uint16), self.config.p_hrtc.wfs_payload_size)
            self.frame.framecounter = self.framecounter
            self.publisher.publish(self.hrtc_host, self.frame)
            tmp = self.subscriber.getFrame({self.local_host:1})
            self.com = np.array(tmp[self.local_host][0])
        else:
            self.frame_shm.send(wfs_cropped)
            self.com = self.com_shm.recv()
        self.framecounter += 1
    
    def apply_control(self):
        """ Apply DM command to DM
        """
        self._dms.set_command(self.com)