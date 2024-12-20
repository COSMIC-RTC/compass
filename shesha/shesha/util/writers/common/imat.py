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
import numpy as np


def compose_imat(sup, *, compose_type="controller", controller_id=0):
    """Return an interaction matrix

    return either the specified controller interaction matrix (if compose_type="controller")
    or an imat composed of all controller interaction matrices (if compose_type="splitTomo")

    Args:
        sup : (compasSSupervisor) : supervisor

    Kargs:
        compose_type : (str) : (optional), default "controller" possibility to specify split tomography case ("controller" or "splitTomo")

        controller_id : (int) : (optional), default 0 controller index

    Returns:
        imat : (np.ndarray[ndim=1, dtype=np.float32]) : interaction matrix
    """
    if compose_type == "controller":
        return sup.rtc.get_interaction_matrix(controller_id)
    elif compose_type == "splitTomo":
        n_actu = 0
        n_meas = 0
        for c in range(len(sup.config.p_controllers)):
            imShape = sup.rtc.get_interaction_matrix(c).shape
            n_meas += imShape[0]
            n_actu += imShape[1]
        imat = np.zeros((n_meas, n_actu))
        n_meas = 0
        n_actu = 0
        for c in range(len(sup.config.p_controllers)):
            im = sup.rtc.get_interaction_matrix(c)
            imat[n_meas : n_meas + im.shape[0], n_actu : n_actu + im.shape[1]] = np.copy(im)
            n_meas += im.shape[0]
            n_actu += im.shape[1]
        return imat

    else:
        print("Unknown composition type")
