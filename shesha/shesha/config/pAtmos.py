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
from typing import List
import shesha.config.config_setter_utils as csu


class ParamAtmos:
    """Class that contains the configuration for the atmosphere module of the COMPASS simulator.

    This class contains the configuration for the atmosphere module of the COMPASS simulator.

    Attributes:
        nscreens (int): Number of turbulent layers.
        r0 (float): Global r0.
        pupixsize (float): Pupil pixel size (in meters).
        L0 (List[float]): L0 per layers in meters.
        dim_screens (List[float]): Linear size of phase screens.
        alt (List[float]): Altitudes of each layer.
        winddir (List[float]): Wind directions of each layer.
        windspeed (List[float]): Wind speeds of each layer.
        frac (List[float]): Fraction of r0 for each layer.
        deltax (List[float]): x translation speed (in pix / iteration) for each layer.
        deltay (List[float]): y translation speed (in pix / iteration) for each layer.
        seeds (List[int]): RNG Seeds for each layer.
    """

    def __init__(self):
        """Constructor of the ParamAtmos class."""
        self.__nscreens = 0  # Number of turbulent layers
        self.__r0 = None  # Global r0
        self.__pupixsize = None  # Pupil pixel size (in meters)
        self.__L0 = None  # L0 per layers in meters
        self.__dim_screens = None  # Linear size of phase screens
        self.__alt = None  # Altitudes of each layer
        self.__winddir = None  # Wind directions of each layer
        self.__windspeed = None  # Wind speeds of each layer
        self.__frac = None  # Fraction of r0 for each layer
        self.__deltax = None  # x translation speed (in pix / iteration) for each layer
        self.__deltay = None  # y translation speed (in pix / iteration) for each layer
        self.__seeds = None  # RNG Seeds for each layer

    def get_nscreens(self) -> int:
        """Get the number of turbulent layers.

        Returns:
            int: Number of turbulent layers.
        """
        return self.__nscreens

    def set_nscreens(self, n: int) -> None:
        """Set the number of turbulent layers.

        Args:
            n (int): Number of turbulent layers.
        """
        self.__nscreens = csu.enforce_int(n)

    nscreens: int = property(get_nscreens, set_nscreens)

    def get_r0(self) -> float:
        """Get the global r0.

        Returns:
            float: Global r0.
        """
        return self.__r0

    def set_r0(self, r: float) -> None:
        """Set the global r0.

        Args:
            r (float): Global r0.
        """
        self.__r0 = csu.enforce_float(r)

    r0: float = property(get_r0, set_r0)

    def get_pupixsize(self) -> float:
        """Get the pupil pixel size.

        Returns:
            float: Pupil pixel size.
        """
        return self.__pupixsize

    def set_pupixsize(self, xsize: float) -> None:
        """Set the pupil pixel size.

        Args:
            xsize (float): Pupil pixel size.
        """
        self.__pupixsize = csu.enforce_float(xsize)

    pupixsize: float = property(get_pupixsize, set_pupixsize)

    def get_L0(self) -> List[float]:
        """Get the L0 per layers.

        Returns:
            List[float]: L0 for each layer.
        """
        return self.__L0

    def set_L0(self, l0_layers: List[float]) -> None:
        """Set the L0 per layers.

        Args:
            l0_layers (List[float]): L0 for each layer.
        """
        self.__L0 = csu.enforce_array(
            l0_layers, size=self.nscreens, dtype=np.float32, scalar_expand=True
        )

    L0: List[float] = property(get_L0, set_L0)

    def get_dim_screens(self) -> List[float]:
        """Get the size of the phase screens.

        Returns:
            List[float]: Phase screens sizes.
        """
        return self.__dim_screens

    def set_dim_screens(self, size_layers: List[float]) -> None:
        """Set the size of the phase screens.

        Args:
            size_layers (List[float]): Phase screens sizes.
        """
        self.__dim_screens = csu.enforce_array(
            size_layers, size=self.nscreens, dtype=np.int64, scalar_expand=False
        )

    dim_screens: List[float] = property(get_dim_screens, set_dim_screens)

    def get_alt(self) -> List[float]:
        """Get the altitudes of each layer.

        Returns:
            List[float]: Altitudes.
        """
        return self.__alt

    def set_alt(self, h: List[float]) -> None:
        """Set the altitudes of each layer.

        Args:
            h (List[float]): Altitudes.
        """
        self.__alt = csu.enforce_array(h, size=self.nscreens, dtype=np.float32, scalar_expand=False)

    alt: List[float] = property(get_alt, set_alt)

    def get_winddir(self) -> List[float]:
        """Get the wind direction for each layer.

        Returns:
            List[float]: Wind directions.
        """
        return self.__winddir

    def set_winddir(self, wind_layers: List[float]) -> None:
        """Set the wind direction for each layer.

        Args:
            wind_layers (List[float]): Wind directions.
        """
        self.__winddir = csu.enforce_array(
            wind_layers,
            size=self.nscreens,
            dtype=np.float32,
            scalar_expand=True,
        )

    winddir: List[float] = property(get_winddir, set_winddir)

    def get_windspeed(self) -> List[float]:
        """Get the wind speed for each layer.

        Returns:
            List[float]: Wind speeds.
        """
        return self.__windspeed

    def set_windspeed(self, windspeed_layers: List[float]) -> None:
        """Set the wind speed for each layer.

        Args:
            windspeed_layers (List[float]): Wind speeds.
        """
        self.__windspeed = csu.enforce_array(
            windspeed_layers,
            size=self.nscreens,
            dtype=np.float32,
            scalar_expand=True,
        )

    windspeed: List[float] = property(get_windspeed, set_windspeed)

    def get_frac(self) -> List[float]:
        """Get the fraction of r0 for each layer.

        Returns:
            List[float]: Fraction of r0.
        """
        return self.__frac

    def set_frac(self, r0_layers: List[float]) -> None:
        """Set the fraction of r0 for each layer.

        Args:
            r0_layers (List[float]): Fraction of r0.
        """
        self.__frac = csu.enforce_array(
            r0_layers, size=self.nscreens, dtype=np.float32, scalar_expand=True
        )

    frac: List[float] = property(get_frac, set_frac)

    def get_deltax(self) -> List[float]:
        """Get the translation speed on the x-axis for each layer.

        Returns:
            List[float]: The translation speed on the x-axis for each layer.
        """
        return self.__deltax

    def set_deltax(self, deltax_layers: List[float]) -> None:
        """Set the translation speed on the x-axis for each layer.

        Args:
            deltax_layers (List[float]): The translation speed on the x-axis for each layer.
        """
        self.__deltax = csu.enforce_array(
            deltax_layers,
            size=self.nscreens,
            dtype=np.float32,
            scalar_expand=True,
        )

    deltax: List[float] = property(get_deltax, set_deltax)

    def get_deltay(self) -> List[float]:
        """Get the translation speed on the y-axis for each layer.

        Returns:
            List[float]: The translation speed on the y-axis for each layer.
        """
        return self.__deltay

    def set_deltay(self, deltay_layers: List[float]) -> None:
        """Set the translation speed on the y-axis for each layer.

        Args:
            deltay_layers (List[float]): The translation speed on the y-axis for each layer.
        """
        self.__deltay = csu.enforce_array(
            deltay_layers,
            size=self.nscreens,
            dtype=np.float32,
            scalar_expand=True,
        )

    deltay: List[float] = property(get_deltay, set_deltay)

    def get_seeds(self) -> List[int]:
        """Get the seed for each layer.

        Returns:
            List[int]: The seed for each layer.
        """
        return self.__seeds

    def set_seeds(self, seed_layers: List[int]) -> None:
        """Set the seed for each layer.

        Args:
            seed_layers (List[int]): The seed for each layer.
        """
        self.__seeds = csu.enforce_array(
            seed_layers, size=self.nscreens, dtype=np.int64, scalar_expand=True
        )

    seeds: List[int] = property(get_seeds, set_seeds)
