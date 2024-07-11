## @package   shesha.config.pCentroider
## @brief     Class definition for the ParamCentroider class
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


import numpy as np
import shesha.config.config_setter_utils as csu
import shesha.constants as scons


class ParamCentroider:
    """Class that defines the centroider parameters

    Attributes:
        nwfs (int): Index of wfs in y_wfs structure on which we want to do centroiding
        type (str): Type of centroiding cog, tcog, bpcog, wcog, corr
        nslope (int): Number of output slopes
        type_fct (str): Type of ref function gauss, file, model
        weights (np.ndarray): Optional reference function(s) used for centroiding
        nmax (int): Number of brightest pixels
        thresh (float): Threshold
        width (float): Width of the Gaussian
        sizex (int): x-size for inter mat (corr)
        sizey (int): y-size for inter mat (corr)
        interpmat (np.ndarray): Optional reference function(s) used for corr centroiding
        method (int): Optional method used in the pyrhr centroider (
                    0: nosinus global
                    1: sinus global
                    2: nosinus local
                    3: sinus local)
        pyrscale (float): pyrscale = (p_wfs.Lambda * 1e-6 / sim.config.p_tel.diam) * p_wfs.pyr_ampl * CONST.RAD2ARCSEC
        filter_TT (bool): filter_TT = False
    """

    def __init__(self):
        self.__nwfs = None  # Index of wfs in y_wfs structure on which we want to do centroiding
        self.__type = None  # Type of centroiding cog, tcog, bpcog, wcog, corr
        self.__nslope = 0  # Number of output slopes
        self.__type_fct = scons.CentroiderFctType.GAUSS  # Type of ref function gauss, file, model
        self.__weights = None  # Optional reference function(s) used for centroiding
        self.__nmax = None  # Number of brightest pixels
        self.__thresh = 1.0e-4  # Threshold
        self.__width = None  # Width of the Gaussian
        self.__sizex = None  # x-size for inter mat (corr)
        self.__sizey = None  # y-size for inter mat (corr)
        self.__interpmat = None  # Optional reference function(s) used for corr centroiding
        self.__method = 1  # Optional method used in the pyrhr centroider
        self.__pyrscale = 0  # pyrscale = (p_wfs.Lambda * 1e-6 / sim.config.p_tel.diam) * p_wfs.pyr_ampl * CONST.RAD2ARCSEC
        self.__filter_TT = False  # Filter out TT from the slopes (LGS mode)

    def get_nwfs(self) -> int:
        """Get the index of the WFS handled by the centroider

        Returns:
            int: WFS index
        """
        return self.__nwfs

    def set_nwfs(self, n: int) -> None:
        """Set the index of the WFS handled by the centroider

        Args:
            n (int): WFS index
        """
        self.__nwfs = csu.enforce_int(n)

    nwfs: int = property(get_nwfs, set_nwfs)

    def get_nslope(self) -> int:
        """Get the number of slope

        Returns:
            int: Number of slope
        """
        return self.__nslope

    def set_nslope(self, n: int) -> None:
        """Set the number of slope

        Args:
            n (int): Number of slope
        """
        self.__nslope = csu.enforce_int(n)

    nslope: int = property(get_nslope, set_nslope)

    def get_type(self) -> str:
        """Get the centroider type

        Returns:
            str: Type
        """
        return self.__type

    def set_type(self, t: str) -> None:
        """Set the centroider type

        Args:
            t (str): Type
        """
        self.__type = scons.check_enum(scons.CentroiderType, t)

    type: str = property(get_type, set_type)

    def get_type_fct(self) -> str:
        """Get the type of reference function used by the centroider

        Returns:
            str: Type
        """
        return self.__type_fct

    def set_type_fct(self, t: str) -> None:
        """Set the type of reference function used by the centroider

        Args:
            t (str): Type
        """
        self.__type_fct = scons.check_enum(scons.CentroiderFctType, t)

    type_fct: str = property(get_type_fct, set_type_fct)

    def get_weights(self) -> np.ndarray:
        """Get the weights used by a wcog centroider

        Returns:
            np.ndarray: Weights
        """
        return self.__weights

    def set_weights(self, w: np.ndarray) -> None:
        """Set the weights used by a wcog centroider

        Args:
            w (np.ndarray): Weights
        """
        self.__weights = csu.enforce_arrayMultiDim(w, w.shape, dtype=np.float32)

    weights: np.ndarray = property(get_weights, set_weights)

    def get_nmax(self) -> int:
        """Get the nmax pixels used by a bpcog centroider

        Returns:
            int: nmax
        """
        return self.__nmax

    def set_nmax(self, n: int) -> None:
        """Set the nmax pixels used by a bpcog centroider

        Args:
            n (int): nmax
        """
        self.__nmax = csu.enforce_int(n)

    nmax: int = property(get_nmax, set_nmax)

    def get_thresh(self) -> float:
        """Get the threshold used by a tcog centroider

        Returns:
            float: Threshold
        """
        return self.__thresh

    def set_thresh(self, t: float) -> None:
        """Set the threshold used by a tcog centroider

        Args:
            t (float): Threshold
        """
        self.__thresh = csu.enforce_float(t)

    thresh: float = property(get_thresh, set_thresh)

    def get_width(self) -> float:
        """Get the width of the gaussian used by a corr centroider

        Returns:
            float: Width
        """
        return self.__width

    def set_width(self, t: float) -> None:
        """Set the width of the gaussian used by a corr centroider

        Args:
            t (float): Width
        """
        self.__width = csu.enforce_float(t)

    width: float = property(get_width, set_width)

    def get_sizex(self) -> int:
        """Get the x size of the interpolation matrix for corr centroider.

        Returns:
            int: The x size of the interpolation matrix.
        """
        return self.__sizex

    def set_sizex(self, n: int) -> None:
        """Set the x size of the interpolation matrix for corr centroider.

        Args:
            n (int): The x size of the interpolation matrix.
        """
        self.__sizex = csu.enforce_int(n)

    sizex: int = property(get_sizex, set_sizex)

    def get_sizey(self) -> int:
        """Get the y size of the interpolation matrix for corr centroider.

        Returns:
            int: The y size of the interpolation matrix.
        """
        return self.__sizey

    def set_sizey(self, n: int) -> None:
        """Set the y size of the interpolation matrix for corr centroider.

        Args:
            n (int): The y size of the interpolation matrix.
        """
        self.__sizey = csu.enforce_int(n)

    sizey: int = property(get_sizey, set_sizey)

    def get_interpmat(self) -> np.ndarray:
        """Get the interpolation matrix for corr centroider.

        Returns:
            np.ndarray: The interpolation matrix.
        """
        return self.__interpmat

    def set_interpmat(self, imap: np.ndarray) -> None:
        """Set the interpolation matrix for corr centroider.

        Args:
            imap (np.ndarray): The interpolation matrix.
        """
        self.__interpmat = csu.enforce_arrayMultiDim(imap, imap.shape, dtype=np.float32)

    interpmat: np.ndarray = property(get_interpmat, set_interpmat)

    def get_method(self) -> int:
        """Get the method used by a pyr centroider.

        Returns:
            int: The method used by a pyr centroider.
                0: nosinus global
                1: sinus global
                2: nosinus local
                3: sinus local
        """
        return self.__method

    def set_method(self, n: int) -> None:
        """Set the method used by a pyr centroider.

        Args:
            n (int): The method used by a pyr centroider.
                0: nosinus global
                1: sinus global
                2: nosinus local
                3: sinus local
        """
        self.__method = csu.enforce_int(n)

    method: int = property(get_method, set_method)

    def get_pyrscale(self) -> float:
        """Get the pyrscale value.

        Returns:
            float: The pyrscale value.
        """
        return self.__pyrscale

    def set_pyrscale(self, t: float) -> None:
        """Set the pyrscale value.

        Args:
            t (float): The pyrscale value.
        """
        self.__pyrscale = csu.enforce_float(t)

    pyrscale: float = property(get_pyrscale, set_pyrscale)

    def get_filter_TT(self) -> bool:
        """Get the filter_TT flag.

        Returns:
            bool: The filter_TT flag.
        """
        return self.__filter_TT

    def set_filter_TT(self, t: bool) -> None:
        """Set the filter_TT flag.

        Args:
            t (bool): The filter_TT flag.
        """
        self.__filter_TT = csu.enforce_or_cast_bool(t)

    filter_TT: bool = property(get_filter_TT, set_filter_TT)
