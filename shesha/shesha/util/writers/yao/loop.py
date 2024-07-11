## @package   shesha.util
## @brief     Shesha utilities
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
def write_loop(file_name, loop, controller):
    """Write (append) AO loop parameters to file for YAO

    Args:
        file_name (str) : yao parameter file name

        loop : (ParamLoop) : compass loop parameters

        controller : (ParamController) : compass controller parameters
    """
    f=open(file_name,"a+")
    f.write("\n\n//------------------------------")
    f.write("\n//LOOP  parameters")
    f.write("\n//------------------------------")
    f.write("\nloop.method     = " + "\"none\"" + ";")
    f.write("\nloop.leak       = " + str(0.001) + ";")
    f.write("\nloop.gain       = " + str(controller.gain) + ";")
    f.write("\nloop.framedelay = " + str(controller.delay+1) + ";") # delay_yao = delay_compass + 1
    f.write("\nloop.niter      = " + str(loop.niter) + ";")
    f.write("\nloop.ittime     = " + str(loop.ittime) + ";")
    f.write("\nloop.skipevery  = " + str(100000) + ";")
    f.write("\nloop.startskip  = " + str(30) + ";")
    f.write("\nloop.skipby     = " + str(5000) + ";")

    f.close()
