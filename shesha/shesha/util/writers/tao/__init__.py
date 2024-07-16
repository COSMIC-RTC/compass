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
from shesha.util.writers.tao.sysParams import write_json_sys_param
from shesha.util.writers.tao.atmParams import write_json_atm_param
from shesha.util.writers import common

def write_parfiles(sup, *, file_name_sys="./sysParams.json",
    file_name_atm="./atmParams.json", file_name_data="sys-input.fits",
    wfss_indices=None, ts=False, dms_indices=None, imat_type="controller",
    controller_id=-1,influ_index=0):
    """write the parameter files for SHIPS

    Args:
        sup : (CompassSupervisor) : supervisor

    Kargs:
        file_name_sys : (str) : AO system parameter file name (default = sysParams.json)

        file_name_atm : (str) : atmospheric parameter file name (default = atmParams.json)

        file_name_data : (str) : data fits file name (default = sys-input.fits), contains sub-apertures and actuator position etc

        wfss_indices : (list(int)) : list of wfs to write to file

        ts : (bool) : write truth sensor to file

        dms_indices : (list(int)) : list of dm to write to file

        imat_type : (str) : (optional), default "controller" use of regular controller or split tomography (among "controller", "splitTomo")

        controller_id : (int) : index of te controller (default : all)

        influ_index : (int) : actuator index to get the influence function from
    """

    write_json_sys_param(sup, wfss_indices=wfss_indices, ts=ts,
        dms_indices=dms_indices,file_name=file_name_sys)

    write_json_atm_param(sup, file_name=file_name_atm)

    common.write_data(file_name_data, sup, wfss_indices=wfss_indices,
        dms_indices=dms_indices, controller_id=controller_id,
        influ=influ_index, compose_type="controller")