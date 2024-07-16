## @package   shesha.supervisor
## @brief     User layer for initialization and execution of a COMPASS simulation
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @version   5.5.0
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


from shesha.init.rtc_init import rtc_init

from shesha.supervisor.components.rtc.rtcAbstract import RtcAbstract, carmaWrap_context


class RtcCompass(RtcAbstract):
    """ RTC handler for compass simulation
    """

    def __init__(self, context: carmaWrap_context, config, tel, wfs, dms, atm):
        """ Initialize a RtcCompass component for rtc related supervision

        Args:
            context : (carmaContext) : CarmaContext instance

            config : (config module) : Parameters configuration structure module

            tel: (Telescope) : Telescope object

            wfs: (Sensors) : Sensors object

            dms: (Dms) : Dms object

            atm: (Atmos) : Atmos object
        """
        RtcAbstract.__init__(self, context, config)
        self.rtc_init(tel, wfs, dms, atm)

    def rtc_init(self, tel, wfs, dms, atm):
        """ Initialize a RtcCompass component for rtc related supervision

        Args:
            tel: (Telescope) : Telescope object

            wfs: (Sensors) : Sensors object

            dms: (Dms) : Dms object

            atm: (Atmos) : Atmos object
        """
        self._rtc = rtc_init(self._context, tel._tel, wfs._wfs, dms._dms, atm._atmos,
                             self._config.p_wfss, self._config.p_tel,
                             self._config.p_geom, self._config.p_atmos,
                             self._config.p_loop.ittime, self._config.p_centroiders,
                             self._config.p_controllers, self._config.p_dms)
