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
import numpy as np

def write_atm(file_name, atm, screen_file, zenithangle):
    """Write (append) atmospheric parameters to file for YAO use

    Args:
        file_name : (str) : name of the file to append the parameter to

        atm : (ParamAtmos) : compass atmospheric parameters. Note that
            atm.winddir is transformed

        screen_file : (str) : path to the yao turbulent screen files. Note
            that the string is passed through raw (without quotes around it)
            in order to use yorick variables in the path name (e.g., Y_USER).
    """
    f = open(file_name,"a+")
    f.write("\n\n//------------------------------")
    f.write("\n//ATM  parameters")
    f.write("\n//------------------------------")

    f.write("\nr0              =" + str(atm.r0) + "; //qt 500 nm")
    f.write("\natm.dr0at05mic  = tel.diam/r0;")

    indexList = '"1"'
    for i in range(2, atm.nscreens + 1):
        indexList += ',"' + str(i) + '"'
    f.write("\natm.screen = &(" + screen_file + "+["+indexList + \
            "]+\".fits\")")
    f.write("\natm.layerspeed  = &(" + np.array2string(atm.windspeed / np.cos(np.pi*zenithangle/180), \
            separator=',', max_line_width=300) + ");")
    f.write("\natm.layeralt    = &(" + np.array2string(atm.alt * np.cos(np.pi*zenithangle/180), \
            separator=',', max_line_width=300) + ");")
    f.write("\natm.layerfrac   = &(" + np.array2string(atm.frac, \
            separator=',', max_line_width=300) + ");")
    f.write("\natm.winddir     = &(" + np.array2string(-(atm.winddir+90)%360, \
            separator=',', max_line_width=300) + ");")
    f.close()
