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
def write_targets(file_name, tars, *, sub_system=1):
    """Write (append) target parameter to file for YAO use for a single dm

    Args:
        file_name : (str) : name of the file to append the parameter to

        tars : (list[ParamTarget]) : compass target parameters list

    Kwargs:
        sub_system : (int) : (optional), default 1 yao sub system index
    """
    f=open(file_name,"a+")
    f.write("\n\n//------------------------------")
    f.write("\n//TAR  parameters")
    f.write("\n//------------------------------")

    f.write("\ntarget.lambda       = &(" + np.array2string(np.array( \
            [t.Lambda for t in tars]), separator=',', max_line_width=300) + \
            ");") #&([0.55]);
    f.write("\ntarget.xposition    = &(" + np.array2string(np.array(\
            [t.xpos for t in tars]), separator=',', max_line_width=300) + \
            ");") # &mmsepos_asec1;
    f.write("\ntarget.yposition    = &(" + np.array2string(np.array( \
            [t.ypos for t in tars]), separator=',', max_line_width=300) + \
            ");") # &mmsepos_asec2;
    dispzoom = np.ones((len(tars)))
    f.write("\ntarget.dispzoom     = &(" + np.array2string(dispzoom, \
            separator=',',max_line_width=300) + ") ; // not set by compass")
            #+ np.array2string(np.array([t.mag for t in tars]),separator=',',max_line_width=300)+";)") #  &array(5.0,numberof(mmsepos_asec1));
