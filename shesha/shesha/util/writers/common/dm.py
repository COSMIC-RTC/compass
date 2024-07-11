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
def get_actu_pos_pixel(dm):
    """return the coordinates in pixel of a given DM actuators

    Args:
        dm : (ParamDm) : Dm to get the actuators position from

    Returns:
        xpos : (np.ndarray[ndim=1, dtype=np.float32]) : actuators positions along axis x

        ypos : (np.ndarray[ndim=1, dtype=np.float32]) : actuators positions along axis y
    """

    return dm._xpos+1, dm._ypos+1


def get_actu_pos_meter(sup, dm_id):
    """return the coordinates in meters of a given DM actuators

    Args:
        sup : (compasSSupervisor) : supervisor

        dm_id : (int) : index of the DM

    Returns:
        xpos : (np.ndarray[ndim=1, dtype=np.float32]) : actuators positions along axis x

        ypos : (np.ndarray[ndim=1, dtype=np.float32]) : actuators positions along axis y
    """

    config = sup.config
    dm=config.p_dms[dm_id]
    geom = config.p_geom
    valid_X = ( dm._xpos - geom.get_cent() ) * geom.get_pixsize()
    valid_Y = ( dm._ypos - geom.get_cent() ) * geom.get_pixsize()
    return valid_X, valid_Y


def dm_to_json(dm, geom):
    """return a json description of a dm

    Args:
        dm : (ParamDm) : dm to represent as json
    """
    dm_json = {
        "n_actu" : dm.get_nact(),
        "h" : dm.get_alt(),
        "coupling" : dm.get_coupling(),
        "shift_x" : dm.get_dx() * geom.get_pixsize(),
        "shift_y" : dm.get_dy() * geom.get_pixsize(),
        "theta" : dm.get_theta()
    }
    return dm_json