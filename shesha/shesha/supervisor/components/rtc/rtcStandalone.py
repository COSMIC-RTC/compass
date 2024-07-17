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

from shesha.init.rtc_init import rtc_standalone

from shesha.supervisor.components.rtc.rtcAbstract import (
    RtcAbstract,
    carma_context,
)


class RtcStandalone(RtcAbstract):
    """RTC handler for compass standalone"""

    def __init__(
        self,
        context: carma_context,
        config,
        nwfs: int,
        nvalid: list,
        nactu: int,
        centroider_type: list,
        delay: list,
        offset: list,
        scale: list,
    ):
        """Initialize a RtcStandalone component for rtc related supervision

        Args:
            context : (carmaContext) : CarmaContext instance

            config : (config module) : Parameters configuration structure module

            nwfs: (int): number of wavefront sensors

            nvalid: (int): number of valid measures as input

            nactu: (int): number of actuators as output

            centroider_type: (list): type of centroiders

            delay: (list): delay of each controller

            offset: (list): offset added in the cog computation of each WFS

            scale: (list): scale factor used in the cog computation of each WFS
        """
        RtcAbstract.__init__(self, context, config)

        self.rtc_init(nwfs, nvalid, nactu, centroider_type, delay, offset, scale)

    def rtc_init(
        self,
        nwfs: int,
        nvalid: list,
        nactu: int,
        centroider_type: list,
        delay: list,
        offset: list,
        scale: list,
    ):
        """Initialize a RtcStandalone component for rtc related supervision

        Args:
            nwfs: (int): number of wavefront sensors

            nvalid: (int): number of valid measures as input

            nactu: (int): number of actuators as output

            centroider_type: (list): type of centroiders

            delay: (list): delay of each controller

            offset: (list): offset added in the cog computation of each WFS

            scale: (list): scale factor used in the cog computation of each WFS
        """
        self._rtc = rtc_standalone(
            self._context,
            nwfs,
            nvalid,
            nactu,
            centroider_type,
            delay,
            offset,
            scale,
        )
