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
def write_gs(file_name, zero_point, lgs_return_per_watt, zenith_angle):
    """Write (append) guide stars parameters to file for YAO

    Args:
        file_name : (str) : name of the file to append the parameter to

        zero_point : (float) : flux for magnitude 0 (ph/m²/s)

        lgs_return_per_watt : (float) : return per watt factor (ph/cm²/s/W)

        zenith_angle : (float) : zenithal angle (degree)
    """
    f=open(file_name,"a+")
    f.write("\n\n//------------------------------")
    f.write("\n//GS parameters")
    f.write("\n//------------------------------")

    # N.B. YAO zeropoint is in flux/telescope aperture, so this value is scaled
    # in outer function
    f.write("\ngs.zeropoint         = " + str(zero_point)+";")
    f.write("\ngs.lgsreturnperwatt  = " + str(lgs_return_per_watt)+";")
    f.write("\ngs.zenithangle       = " + str(zenith_angle) + ";")
