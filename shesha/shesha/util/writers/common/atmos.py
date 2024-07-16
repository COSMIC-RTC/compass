## @package   shesha.util
## @brief     Shesha utilities
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
def atmos_to_json(atmos, name=""):
    """return a json description of a the atmos

    Args:
        atmos : (ParamAtmos) : compass atmospheric parameters
    """
    json_atm = {
        "nLayer" : atmos.get_nscreens(),
        "r0" : atmos.get_r0(),
        "h" : atmos.get_alt().tolist(),
        "fracCn2" : atmos.get_frac().tolist(),
        "L0" : atmos.get_L0().tolist(),
        "windDir" : atmos.get_winddir().tolist(),
        "windSpeed" : atmos.get_windspeed().tolist()
    }
    if(name != ""):
        json_atm["name"] = name
    return json_atm


def atmos_json_notice():
    notice = {
        "name": "          : profile name",
        "nLayer": "          : number of layers in the turbulent profile",
        "r0": "r0 at 500 nm (fried parameter)",
        "h": " meter    : (list) altitude of each layer",
        "fracCn2": " percent  : (list) cn2 fraction of each layer",
        "L0":  "meter    : (list) outer scale of each layer",
        "windDir": " degree   : (list) wind sirection of each layer",
        "windSpeed": " meter/s  : (list) wind speed of each layer"
    }
    return notice