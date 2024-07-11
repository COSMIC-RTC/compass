## @package   shesha.config.pCoronagraph
## @brief     Class for the coronagraph configuration
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


import shesha.config.config_setter_utils as csu
import numpy as np
from typing import List

class ParamCoronagraph:
    """Class for the coronagraph configuration

    Attributes:
        _type (str): coronograph type
        _asterix_parfile (str): asterix parfile path
        _asterix_datadir (str): asterix datadir path
        _wavelength_0 (float): central wavelength
        _delta_wav (float): spectral bandwidth
        _nb_wav (int): number of wavelengths
        _wav_vec (np.ndarray): wavelengths array
        _apodizer (np.ndarray): apodizer
        _apodizer_name (str): apodizer keyword or path
        _focal_plane_mask (list of np.ndarray): focal plane mask
        _focal_plane_mask_name (str): focal plane mask keyword or path
        _fpm_sampling (float): focal plane mask sampling
        _lyot_fpm_radius (float): classical Lyot fpm radius
        _dim_fpm (int): fpm support size in pixel
        _babinet_trick (bool): Babinet's trick flag
        _lyot_stop (np.ndarray): Lyot stop pupil
        _lyot_stop_name (str): Lyot stop keyword or path
        _dim_image (int): image size in pixel
        _image_sampling (float): image sampling
    """

    def __init__(self):
        self.__type = None  # Type of coronograph

        self.__asterix_parfile = None  # ASTERIX parameter file path
        self.__asterix_datadir = None  # ASTERIX data directory path

        self.__wavelength_0 = None  # Central wavelength in the coronagraph
        self.__delta_wav = 0  # optional : spectral bandwidth
        self.__nb_wav = 1  # optional : number of simulated wavelength in the spectral bandwidth
        self.__wav_vec = None  # unsettable : array of simulated wavelengths

        self.__apodizer = None  # Apodizer pupil
        self.__apodizer_name = None  # optional : apodizer string name or user path

        self.__focal_plane_mask = None  # Focal plane mask complex amplitudes
        self.__focal_plane_mask_name = None  # Focal plane mask string name or user path
        self.__fpm_sampling = None  # optional : size of lambda / D in the fpm plane, in pixel unit
        self.__lyot_fpm_radius = (
            None  # Focal plane mask radius in lamda / D unit, for a classical Lyot fpm only
        )
        self.__dim_fpm = None  # Size of the focal plane mask in pixel
        self.__babinet_trick = False  # Flag for using Babinet's trick

        self.__lyot_stop = None  # Lyot pupil
        self.__lyot_stop_name = None  # optional : Lyot stop string name or user path

        self.__dim_image = None  # Size of the science image in pixel
        self.__image_sampling = None  # Size of lambda / D in the image plane, in pixel unit

    def get_type(self) -> str:
        """Get the coronograph type.

        Returns:
            str: The coronograph type.
        """
        return self.__type

    def set_type(self, t: str) -> None:
        """Set the coronograph type.

        Args:
            t (str): The coronograph type.
        """
        self.__type = t

    _type: str = property(get_type, set_type)

    def get_asterix_parfile(self) -> str:
        """Get the path of asterix parfile.

        Returns:
            str: The asterix parfile path.
        """
        return self.__asterix_parfile

    def set_asterix_parfile(self, f: str) -> None:
        """Set the path of asterix parfile.

        Args:
            f (str): The asterix parfile path.
        """
        self.__asterix_parfile = f

    _asterix_parfile: str = property(get_asterix_parfile, set_asterix_parfile)

    def get_asterix_datadir(self) -> str:
        """Get the path of asterix datadir.

        Returns:
            str: The asterix datadir path.
        """
        return self.__asterix_datadir

    def set_asterix_datadir(self, f: str) -> None:
        """Set the path of asterix datadir.

        Args:
            f (str): The asterix datadir path.
        """
        self.__asterix_datadir = f

    _asterix_datadir: str = property(get_asterix_datadir, set_asterix_datadir)

    def get_wavelength_0(self) -> float:
        """Get the central wavelength in the coronagraph.

        Returns:
            float: The central wavelength.
        """
        return self.__wavelength_0

    def set_wavelength_0(self, w: float) -> None:
        """Set the central wavelength in the coronagraph.

        Args:
            w (float): The central wavelength.
        """
        self.__wavelength_0 = w

    _wavelength_0: float = property(get_wavelength_0, set_wavelength_0)

    def get_delta_wav(self) -> float:
        """Get the spectral bandwidth.

        Returns:
            float: The spectral bandwidth.
        """
        return self.__delta_wav

    def set_delta_wav(self, w: float) -> None:
        """Set the spectral bandwidth.

        Args:
            w (float): The spectral bandwidth.
        """
        self.__delta_wav = w

    _delta_wav: float = property(get_delta_wav, set_delta_wav)

    def get_nb_wav(self) -> int:
        """Get the number of simulated wavelength in the spectral bandwidth.

        Returns:
            int: The number of wavelengths.
        """
        return self.__nb_wav

    def set_nb_wav(self, n: int) -> None:
        """Set the number of simulated wavelength in the spectral bandwidth.

        Args:
            n (int): The number of wavelengths.
        """
        self.__nb_wav = csu.enforce_int(n)

    _nb_wav: int = property(get_nb_wav, set_nb_wav)

    def get_wav_vec(self) -> np.ndarray:
        """Get the wavelengths array.

        Returns:
            np.ndarray: The wavelengths array.
        """
        return self.__wav_vec

    def set_wav_vec(self, w: np.ndarray) -> None:
        """Set the wavelengths array.

        Args:
            w (np.ndarray): The wavelengths array.
        """
        self.__wav_vec = w

    _wav_vec: np.ndarray = property(get_wav_vec, set_wav_vec)

    def get_apodizer(self) -> np.ndarray:
        """Get the apodizer pupil.

        Returns:
            np.ndarray: The apodizer pupil.
        """
        return self.__apodizer

    def set_apodizer(self, apod: np.ndarray) -> None:
        """Set the apodizer pupil.

        Args:
            apod (np.ndarray): The apodizer pupil.
        """
        self.__apodizer = apod

    _apodizer: np.ndarray = property(get_apodizer, set_apodizer)

    def get_apodizer_name(self) -> str:
        """Get the apodizer keyword or user path.

        Returns:
            str: The apodizer keyword or path.
        """
        return self.__apodizer_name

    def set_apodizer_name(self, apod: str) -> None:
        """Set the apodizer keyword or user path.

        Args:
            apod (str): The apodizer keyword or path.
        """
        self.__apodizer_name = apod

    _apodizer_name: str = property(get_apodizer_name, set_apodizer_name)

    def get_focal_plane_mask(self) -> List[np.ndarray]:
        """Get the focal plane mask complex amplitudes.

        Returns:
            List[np.ndarray]: The focal plane mask.
        """
        return self.__focal_plane_mask

    def set_focal_plane_mask(self, fpm: List[np.ndarray]) -> None:
        """Set the focal plane mask complex amplitudes.

        Args:
            fpm (List[np.ndarray]): The focal plane mask.
        """
        self.__focal_plane_mask = fpm

    _focal_plane_mask: List[np.ndarray] = property(get_focal_plane_mask, set_focal_plane_mask)

    def get_focal_plane_mask_name(self) -> str:
        """Get the focal plane mask keyword or user path.

        Returns:
            str: The focal plane mask keyword or path.
        """
        return self.__focal_plane_mask_name

    def set_focal_plane_mask_name(self, fpm: str) -> None:
        """Set the focal plane mask keyword or user path.

        Args:
            fpm (str): The focal plane mask keyword or path.
        """
        self.__focal_plane_mask_name = fpm

    _focal_plane_mask_name: str = property(get_focal_plane_mask_name, set_focal_plane_mask_name)

    def get_fpm_sampling(self) -> float:
        """Get the sampling in the focal plane mask.

        Returns:
            float: The sampling in the focal plane mask (size of lambda / D in pixel units).
        """
        return self.__fpm_sampling

    def set_fpm_sampling(self, sp: float) -> None:
        """Set the sampling in the focal plane mask.

        Args:
            sp (float): The sampling in the focal plane mask (size of lambda / D in pixel units).
        """
        self.__fpm_sampling = sp

    _fpm_sampling: float = property(get_fpm_sampling, set_fpm_sampling)

    def get_lyot_fpm_radius(self) -> float:
        """Get the radius of the classical Lyot focal plane mask
        in lambda / D units

        Returns:
            float: The radius of the classical Lyot focal plane mask
        """
        return self.__lyot_fpm_radius

    def set_lyot_fpm_radius(self, r: float) -> None:
        """Set the radius of the classical Lyot focal plane mask
        in lambda / D units

        Args:
            r (float): The radius of the classical Lyot focal plane mask
        """
        self.__lyot_fpm_radius = r

    _lyot_fpm_radius: float = property(get_lyot_fpm_radius, set_lyot_fpm_radius)

    def get_dim_fpm(self) -> int:
        """Get the size of the focal plane mask support in pixel units

        Returns:
            int: fpm support size in pixel
        """
        return self.__dim_fpm

    def set_dim_fpm(self, n: int) -> None:
        """Set the size of the focal plane mask support in pixel units

        Args:
            n (int): fpm support size in pixel
        """
        self.__dim_fpm = csu.enforce_int(n)

    _dim_fpm: int = property(get_dim_fpm, set_dim_fpm)

    def get_babinet_trick(self) -> bool:
        """Get the Babinet's trick flag

        Returns:
            bool: Babinet's trick flag
        """
        return self.__babinet_trick

    def set_babinet_trick(self, b: bool) -> None:
        """Set the Babinet's trick flag

        Args:
            b (bool): Babinet's trick flag
        """
        self.__babinet_trick = csu.enforce_or_cast_bool(b)

    _babinet_trick: bool = property(get_babinet_trick, set_babinet_trick)

    def get_lyot_stop(self) -> np.ndarray:
        """Get the Lyot stop pupil

        Returns:
            np.ndarray: Lyot stop pupil
        """
        return self.__lyot_stop

    def set_lyot_stop(self, ls: np.ndarray) -> None:
        """Set the Lyot stop pupil

        Args:
            ls (np.ndarray): Lyot stop pupil
        """
        self.__lyot_stop = ls

    _lyot_stop: np.ndarray = property(get_lyot_stop, set_lyot_stop)

    def get_lyot_stop_name(self) -> str:
        """Get the Lyot stop keyword or user path

        Returns:
            str: Lyot stop keyword or path
        """
        return self.__lyot_stop_name

    def set_lyot_stop_name(self, ls: str) -> None:
        """Set the Lyot stop keyword or user path

        Args:
            ls (str): Lyot stop keyword or path
        """
        self.__lyot_stop_name = ls

    _lyot_stop_name: str = property(get_lyot_stop_name, set_lyot_stop_name)

    def get_dim_image(self) -> int:
        """Get the size of the science image in pixel

        Returns:
            int: image size in pixel
        """
        return self.__dim_image

    def set_dim_image(self, n: int) -> None:
        """Set the size of the science image in pixel

        Args:
            n (int): image size in pixel
        """
        self.__dim_image = csu.enforce_int(n)

    _dim_image: int = property(get_dim_image, set_dim_image)

    def get_image_sampling(self) -> float:
        """Get the sampling in the image
        sampling = size of lambda / D in pixel units

        Returns:
            float: image sampling
        """
        return self.__image_sampling

    def set_image_sampling(self, sp: float) -> None:
        """Set the sampling in the image
        sampling = size of lambda / D in pixel units

        Args:
            sp (float): image sampling
        """
        self.__image_sampling = sp

    _image_sampling: float = property(get_image_sampling, set_image_sampling)
