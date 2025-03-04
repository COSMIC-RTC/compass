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

from shesha.init.dm_init import dm_init
import numpy as np
from typing import Tuple


class DmCompass(object):
    """DM handler for compass simulation

    Attributes:
        _dms : (sutra.Dms) : Sutra dms instance

        _context : (carmaContext) : CarmaContext instance

        _config : (config module) : Parameters configuration structure module
    """

    def __init__(self, context, config):
        """Initialize a DmCompass component for DM related supervision

        Args:
            context : (carmaContext) : CarmaContext instance

            config : (config module) : Parameters configuration structure module
        """
        self._context = context
        self._config = config  # Parameters configuration coming from supervisor init
        print("->dms init")
        self._dms = dm_init(
            self._context,
            self._config.p_dms,
            self._config.p_tel,
            self._config.p_geom,
            self._config.p_wfss,
        )

    def set_command(
        self,
        commands: np.ndarray,
        *,
        dm_index: int = None,
        shape_dm: bool = True,
    ) -> None:
        """Immediately sets provided command to DMs - does not affect integrator

        Args:
            commands : (np.ndarray) : commands vector to apply

        Kwargs:
            dm_index : (int) : Index of the DM to set. If None (default), set all the DMs.
                               In that case, provided commands vector must have a size equal
                               to the sum of all the DMs actuators.
                               If index_dm is set, size must be equal to the number of actuator
                               of the specified DM.

            shape_dm : (bool) : If True (default), immediately apply the given commands on the DMs
        """
        if dm_index is None:
            self._dms.set_full_com(commands, shape_dm)
        else:
            self._dms.d_dms[dm_index].set_com(commands, shape_dm)

    def set_one_actu(self, dm_index: int, nactu: int, *, ampli: float = 1) -> None:
        """Push the selected actuator

        Args:
            dm_index : (int) : DM index

            nactu : (int) : actuator index to push

        Kwargs:
            ampli : (float) : amplitude to apply. Default is 1 volt
        """
        self._dms.d_dms[dm_index].comp_oneactu(nactu, ampli)

    def get_influ_function(self, dm_index: int) -> np.ndarray:
        """Returns the influence function cube for the given dm

        Args:
            dm_index : (int) : index of the DM

        Returns:
            influ : (np.ndarray) : Influence functions of the DM dm_index
        """
        return self._config.p_dms[dm_index]._influ

    def get_influ_function_ipupil_coords(self, dm_index: int) -> Tuple[np.ndarray, np.ndarray]:
        """Returns the lower left coordinates of the influ function support in the ipupil coord system

        Args:
            dm_index : (int) : index of the DM

        Returns:
            coords : (tuple) : (i, j)
        """
        i1 = self._config.p_dms[dm_index]._i1  # i1 is in the dmshape support coords
        j1 = self._config.p_dms[dm_index]._j1  # j1 is in the dmshape support coords
        ii1 = i1 + self._config.p_dms[dm_index]._n1  # in  ipupil coords
        jj1 = j1 + self._config.p_dms[dm_index]._n1  # in  ipupil coords
        return ii1, jj1

    def reset_dm(self, dm_index: int = -1) -> None:
        """Reset the specified DM or all DMs if dm_index is -1

        Kwargs:
            dm_index : (int) : Index of the DM to reset
                                         Default is -1, i.e. all DMs are reset
        """
        if dm_index == -1:  # All Dms reset
            for dm in self._dms.d_dms:
                dm.reset_shape()
        else:
            self._dms.d_dms[dm_index].reset_shape()

    def get_dm_shape(self, indx: int) -> np.ndarray:
        """Return the current phase shape of the selected DM

        Args:
            indx : (int) : Index of the DM

        Returns:
            dm_shape : (np.ndarray) : DM phase screen

        """
        return np.array(self._dms.d_dms[indx].d_shape)

    def set_dm_registration(
        self,
        dm_index: int,
        *,
        dx: float = None,
        dy: float = None,
        theta: float = None,
        G: float = None,
    ) -> None:
        """Set the registration parameters for DM #dm_index

        Args:
            dm_index : (int) : DM index

        Kwargs:
            dx : (float) : X axis registration parameter [meters]. If None, re-use the last one

            dy : (float) : Y axis registration parameter [meters]. If None, re-use the last one

            theta : (float) : Rotation angle parameter [rad]. If None, re-use the last one

            G : (float) : Magnification factor. If None, re-use the last one
        """
        if dx is not None:
            self._config.p_dms[dm_index].set_dx(dx)
        if dy is not None:
            self._config.p_dms[dm_index].set_dy(dy)
        if theta is not None:
            self._config.p_dms[dm_index].set_theta(theta)
        if G is not None:
            self._config.p_dms[dm_index].set_G(G)

        self._dms.d_dms[dm_index].set_registration(
            self._config.p_dms[dm_index].dx / self._config.p_geom._pixsize,
            self._config.p_dms[dm_index].dy / self._config.p_geom._pixsize,
            self._config.p_dms[dm_index].theta,
            self._config.p_dms[dm_index].G,
        )
