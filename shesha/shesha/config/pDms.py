## @package   shesha.config.pDms
## @brief     Class definition for the DM parameters
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
import shesha.constants as scons
import shesha.config.config_setter_utils as csu


class ParamDm:
    """Class definition for the DM parameters

    This class contains the parameters for the DM.

    Attributes:
        ap (float): The value of ap.
        nfunc (int): The value of nfunc.
        pzt_extent (int): The extent of pzt dm.
        segmented_mirror (bool): The flag indicating whether the mirror is segmented.
        influ_type (str): The influence function type for the pzt DM.
        type_kl (str): The type of KL used for computation.
        type (str): Type of dm.
        type_pattern (str): Type of pattern.
        file_influ_fits (str): Filename for influ hdf5 file.
        center_name (str): Center name in hdf5.
        cube_name (str): Influence function cube name in hdf5.
        x_name (str): x coord name in hdf5.
        y_name (str): y coord name in hdf5.
        influ_res (str): Influence resolution name in hdf5.
        diam_dm (float): Diameter of the tel pupil projected on the dm plane.
        diam_dm_proj (float): Diameter of the dm pupil projected on the tel pupil plane.
        nact (int): Number of actuators in the dm.
        margin_out (float): Margin for outside actuator select.
        margin_in (float): Margin for inside actuator select (central obstruction).
        alt (float): Conjugaison altitude.
        thresh (float): Threshold on response for selection.
        keep_all_actu (bool): Keep all actuator flag.
        coupling (float): Actuators coupling.
        unitpervolt (float): Influence function sensitivity in unit/volt.
        push4imat (float): Nominal voltage for imat.
        dim_screen (int): Phase screen dimension.
        nkl (int): Number of KL modes used for computation of covmat in case of minimum variance controller.
        outscl (float): Outer scale in units of telescope diam for Karman KL.
        nr (int): Number of radial points.
        npp (int): Number of elements.
        ord (np.ndarray): The radial orders of the basis.
        rabas (np.ndarray): The radial array of the basis.
        azbas (np.ndarray): The azimuthal array of the basis.
        ncp (int): Dim of grid.
        cr (np.ndarray): Radial coord in cartesien grid.
        cp (np.ndarray): Phi coord in cartesien grid.
        ap (float): The value of ap.
    """

    def __init__(self):
        # DM properties
        self.__nact = 0  # linear number of actuators across the pupil diameter
        self.__alt = 0.0  # DM conjugation altitude
        self.__thresh = 0.0  # Threshold on response for selection
        self.__keep_all_actu = False  # if True, don't mask actu by pupil
        self.__coupling = 0.2  # Actuator coupling (< .3)
        self.__gain = 1.0  # Actuator gains
        self.__pupoffset = np.array([0, 0])
        self.__dim_screen = 0  # Phase screen dimension
        # Global offset in pupil (x,y) of the whole actuator pattern

        self.__unitpervolt = 0.01
        # Influence function sensitivity in unit/volt. Optional [0.01]
        # Stackarray: mic/volt, Tip-tilt: arcsec/volt.
        self.__push4imat = 1.0  # nominal voltage for imat

        # Margins for actuator selection
        self.__margin_out = None  # outer margin (pitches) from pupil diameter
        # inner margin (pitches) from central obstruction
        self.__margin_in = 0.0
        self.__pzt_extent = 5.0  # Extent of pzt DM (pitches)
        self.__segmented_mirror = False  # Crop influence functions where spiders are.

        # KL DM
        self.__nfunc = 0
        self.__nkl = 0  # Number of KL for KL dm
        self.__outscl = None  # Outer scale in units of telescope diam for Karman KL
        self.__nr = None  # number of radial points
        self.__npp = None  # number of elements
        self.__ord = None  # the radial orders of the basis
        self.__rabas = None  # the radial array of the basis
        self.__azbas = None  # the azimuthal array of the basis
        self.__ncp = None  # dim of grid
        self.__cr = None  # radial coord in cartesien grid
        self.__cp = None  # phi coord in cartesien grid
        self.__ap = None
        self.__nfunc = 0

        # Hidden variable safe-typed in shesha_constants
        self.__type = None  # Private storage of type
        self.__type_pattern = None  # Private storage of type_pattern
        self.__influ_type = scons.InfluType.DEFAULT  # Private storage of influ_type
        self.__type_kl = scons.KLType.KOLMO  # Private storage for KL type

        # HDF5 storage management
        self.__file_influ_fits = None  # Filename for influ hdf5 file
        self.__center_name = None  # Center name in hdf5
        self.__cube_name = None  # Influence function cube name in hdf5
        self.__x_name = None  # x coord name in hdf5
        self.__y_name = None  # y coord name in hdf5
        self.__influ_res = None  # influence resolution name in hdf5
        self.__diam_dm = None  # diam of the tel pupil projected on the dm plane
        self.__diam_dm_proj = None  # diam of the dm pupil projected on the tel pupil plane

        # PXD cleanup
        # internal kwrd
        self.__pitch = None  # inter-actuator space in pixels
        self.__ntotact = None  # total number of actuators over the full area of the pupil
        self.__influsize = None  # influ function support size
        self.__n1 = None  # position of leftmost pixel in largest support
        self.__n2 = None  # position of rightmost pixel in largest support
        self.__puppixoffset = None
        self.__influ = None  # influence functions
        self.__xpos = None  # x positions of influ functions
        self.__ypos = None  # y positions of influ functions
        self.__i1 = None  # position of bottom left pixel in the largest support
        self.__j1 = None  # position of bottom left pixel in the largest support
        self.__influpos = None  # pixels that contributes to each DM pixel
        self.__ninflu = None  # number of pixels in each influence function
        self.__influstart = None  # np.ndarray - Influence function handling

        # Registration
        self.__G = 1.0  # Magnifying factor
        self.__theta = 0.0  # Rotation angle in the pupil
        self.__dx = 0.0  # X axis misalignment in meters
        self.__dy = 0.0  # Y axis misalignment in meters

    def get_ap(self) -> float:
        """Get ap.

        Returns:
            float: The value of ap.
        """
        return self.__ap

    def set_ap(self, ap: float) -> None:
        """Set ap.

        Args:
            ap (float): The value of ap.
        """
        self.__ap = csu.enforce_arrayMultiDim(ap, (ap.shape[0], ap.shape[1]), dtype=np.float32)

    ap: float = property(get_ap, set_ap)

    def get_nfunc(self) -> int:
        """Get nfunc.

        Returns:
            int: The value of nfunc.
        """
        return self.__nfunc

    def set_nfunc(self, nfunc: int) -> None:
        """Set nfunc.

        Args:
            nfunc (int): The value of nfunc.
        """
        self.__nfunc = csu.enforce_int(nfunc)

    nfunc: int = property(get_nfunc, set_nfunc)

    def get_pzt_extent(self) -> int:
        """Get extent of pzt dm in pitch unit.

        Returns:
            int: The extent of pzt dm.
        """
        return self.__pzt_extent

    def set_pzt_extent(self, p: int) -> None:
        """Set extent of pzt dm in pitch unit.

        Args:
            p (int): The extent of pzt dm.
        """
        self.__pzt_extent = csu.enforce_int(p)

    pzt_extent: int = property(get_pzt_extent, set_pzt_extent)

    def get_segmented_mirror(self) -> bool:
        """Get the flag indicating whether the mirror is segmented.

        Returns:
            bool: The flag indicating whether the mirror is segmented.
        """
        return self.__segmented_mirror

    def set_segmented_mirror(self, b: bool) -> None:
        """Set the flag indicating whether the mirror is segmented.

        Args:
            b (bool): The flag indicating whether the mirror is segmented.
        """
        self.__segmented_mirror = csu.enforce_or_cast_bool(b)

    segmented_mirror: bool = property(get_segmented_mirror, set_segmented_mirror)

    def get_influ_type(self) -> str:
        """Get the influence function type for the pzt DM.

        Returns:
            str: The influence function type.
        """
        return self.__influ_type

    def set_influ_type(self, t: str) -> None:
        """Set the influence function type for the pzt DM.

        Args:
            t (str): The influence function type.
        """
        self.__influ_type = scons.check_enum(scons.InfluType, t)

    influ_type: str = property(get_influ_type, set_influ_type)

    def get_influpos(self) -> np.ndarray:
        """Get the influence function pixels that contribute to each DM pixel.

        Returns:
            np.ndarray: The influence function pixels.
        """
        return self.__influpos

    def set_influpos(self, ip: np.ndarray) -> None:
        """Set the influence function pixels that contribute to each DM pixel.

        Args:
            ip (np.ndarray): The influence function pixels.
        """
        self.__influpos = csu.enforce_array(ip, ip.size, dtype=np.int32)

    _influpos: np.ndarray = property(get_influpos, set_influpos)

    def get_ninflu(self) -> np.ndarray:
        """Get the number of influence functions pixels that contribute
        to each DM pixel.

        Returns:
            np.ndarray: The number of influence functions pixels.
        """
        return self.__ninflu

    def set_ninflu(self, n: np.ndarray) -> None:
        """Set the number of influence functions pixels that contribute
        to each DM pixel.

        Args:
            n (np.ndarray): The number of influence functions pixels.
        """
        self.__ninflu = csu.enforce_array(n, n.size, dtype=np.int32)

    _ninflu: np.ndarray = property(get_ninflu, set_ninflu)

    def get_influstart(self) -> np.ndarray:
        """Get the index where to start a new DM pixel shape in the array influpos
        to each DM pixel.

        Returns:
            np.ndarray: The index where to start a new DM pixel shape.
        """
        return self.__influstart

    def set_influstart(self, n: np.ndarray) -> None:
        """Set the index where to start a new DM pixel shape in the array influpos
        to each DM pixel.

        Args:
            n (np.ndarray): The index where to start a new DM pixel shape.
        """
        self.__influstart = csu.enforce_array(n, n.size, dtype=np.int32)

    _influstart: np.ndarray = property(get_influstart, set_influstart)

    def get_gain(self) -> float:
        """Get the gain to apply to the actuators of the DM.

        Returns:
            float: The gain.
        """
        return self.__gain

    def set_gain(self, g: float) -> None:
        """Set the gain to apply to the actuators of the DM.

        Args:
            g (float): The gain.
        """
        self.__gain = csu.enforce_float(g)

    gain: float = property(get_gain, set_gain)

    def _get_dim_screen(self) -> int:
        """Get the phase screen dimension.

        Returns:
            int: The phase screen dimension.
        """
        return self.__dim_screen

    def _set_dim_screen(self, n: int) -> None:
        """Set the phase screen dimension.

        Args:
            n (int): The phase screen dimension.
        """
        self.__dim_screen = csu.enforce_int(n)

    _dim_screen: int = property(_get_dim_screen, _set_dim_screen)

    def get_nkl(self) -> int:
        """Get the number of KL modes used for computation of covmat in case of minimum variance controller.

        Returns:
            int: The number of KL modes.
        """
        return self.__nkl

    def set_nkl(self, n: int) -> None:
        """Set the number of KL modes used for computation of covmat in case of minimum variance controller.

        Args:
            n (int): The number of KL modes.
        """
        self.__nkl = csu.enforce_int(n)

    nkl: int = property(get_nkl, set_nkl)

    def get_type_kl(self) -> str:
        """Get the type of KL used for computation.

        Returns:
            str: The KL types: kolmo or karman.
        """
        return self.__type_kl

    def set_type_kl(self, t: str) -> None:
        """Set the type of KL used for computation.

        Args:
            t (str): The KL types: kolmo or karman.
        """
        self.__type_kl = scons.check_enum(scons.KLType, t)

    type_kl: str = property(get_type_kl, set_type_kl)

    def get_type(self) -> str:
        """Get the dm type

        Returns:
            str: Type of dm
        """
        return self.__type

    def set_type(self, t: str) -> None:
        """Set the dm type

        Args:
            t (str): Type of dm
        """
        self.__type = scons.check_enum(scons.DmType, t)

    type: str = property(get_type, set_type)

    def get_type_pattern(self) -> str:
        """Get the pattern type

        Returns:
            str: Type of pattern
        """
        return self.__type_pattern

    def set_type_pattern(self, t: str) -> None:
        """Set the pattern type

        Args:
            t (str): Type of pattern
        """
        self.__type_pattern = scons.check_enum(scons.PatternType, t)

    type_pattern: str = property(get_type_pattern, set_type_pattern)

    def get_file_influ_fits(self) -> str:
        """Get the name of hdf5 influence file

        Returns:
            str: Hdf5 file influence name
        """
        return self.__file_influ_fits

    def set_file_influ_fits(self, f: str) -> None:
        """Set the name of hdf5 influence file

        Args:
            f (str): Hdf5 file influence name
        """
        self.__file_influ_fits = f

    file_influ_fits: str = property(get_file_influ_fits, set_file_influ_fits)

    def get_center_name(self) -> str:
        """Get the name of hdf5 influence file

        Returns:
            str: Hdf5 file influence name
        """
        return self.__center_name

    def set_center_name(self, f: str) -> None:
        """Set the name of hdf5 influence file

        Args:
            f (str): Hdf5 file influence name
        """
        self.__center_name = f

    center_name: str = property(get_center_name, set_center_name)

    def get_cube_name(self) -> str:
        """Get the name of influence cube in hdf5

        Returns:
            str: Name of influence cube
        """
        return self.__cube_name

    def set_cube_name(self, cubename: str) -> None:
        """Set the name of influence cube in hdf5

        Args:
            cubename (str): Name of influence cube
        """
        self.__cube_name = cubename

    cube_name: str = property(get_cube_name, set_cube_name)

    def get_x_name(self) -> str:
        """Get the name of x coord of influence function in file

        Returns:
            str: Name of x coord of influence
        """
        return self.__x_name

    def set_x_name(self, xname: str) -> None:
        """Set the name of x coord of influence function in file

        Args:
            xname (str): Name of x coord of influence
        """
        self.__x_name = xname

    x_name: str = property(get_x_name, set_x_name)

    def get_y_name(self) -> str:
        """Get the name of y coord of influence function in file

        Returns:
            str: Name of y coord of influence
        """
        return self.__y_name

    def set_y_name(self, yname: str) -> None:
        """Set the name of y coord of influence function in file

        Args:
            yname (str): Name of y coord of influence
        """
        self.__y_name = yname

    y_name: str = property(get_y_name, set_y_name)

    def get_influ_res(self) -> str:
        """Get the name of influence function resolution in file

        Returns:
            str: Name of resolution (meter/pixel) of influence
        """
        return self.__influ_res

    def set_influ_res(self, res: str) -> None:
        """Set the name of influence function resolution in file

        Args:
            res (str): Name of resolution (meter/pixel) of influence
        """
        self.__influ_res = res

    influ_res: str = property(get_influ_res, set_influ_res)

    def get_diam_dm(self) -> float:
        """Get the diameter of the tel pupil projected on the dm plane

        Returns:
            float: Diameter (meters) of the tel pupil projected on the dm plane
        """
        return self.__diam_dm

    def set_diam_dm(self, di: float) -> None:
        """Set the diameter of the tel pupil projected on the dm plane

        Args:
            di (float): Diameter (meters) of the tel pupil projected on the dm plane
        """
        self.__diam_dm = di

    diam_dm: float = property(get_diam_dm, set_diam_dm)

    def get_diam_dm_proj(self) -> float:
        """Get the diameter of the dm pupil projected on the tel pupil plane

        Returns:
            float: Diameter (meters) of the dm pupil projected on the tel pupil plane
        """
        return self.__diam_dm_proj

    def set_diam_dm_proj(self, dp: float) -> None:
        """Set the diameter of the dm pupil projected on the tel pupil plane

        Args:
            dp (float): Diameter (meters) of the dm pupil projected on the tel pupil plane
        """
        self.__diam_dm_proj = dp

    diam_dm_proj: float = property(get_diam_dm_proj, set_diam_dm_proj)

    def get_nact(self) -> int:
        """Get the number of actuator

        Returns:
            int: Number of actuators in the dm
        """
        return self.__nact

    def set_nact(self, n: int) -> None:
        """Set the number of actuator

        Args:
            n (int): Number of actuators in the dm
        """
        self.__nact = csu.enforce_int(n)

    nact: int = property(get_nact, set_nact)

    def get_margin_out(self) -> float:
        """Get the margin for outside actuator select

        Returns:
            float: Unit is actuator pitch (+) for extra (-) for intra
        """
        return self.__margin_out

    def set_margin_out(self, n: float) -> None:
        """Set the margin for outside actuator select

        Args:
            n (float): Unit is actuator pitch (+) for extra (-) for intra
        """
        self.__margin_out = csu.enforce_float(n)

    margin_out: float = property(get_margin_out, set_margin_out)

    def get_margin_in(self) -> float:
        """Get the margin for inside actuator select (central obstruction)

        Returns:
            float: Unit is actuator pitch (+) for extra (-) for intra
        """
        return self.__margin_in

    def set_margin_in(self, n: float) -> None:
        """Set the margin for inside actuator select (central obstruction)

        Args:
            n (float): Unit is actuator pitch (+) for extra (-) for intra
        """
        self.__margin_in = csu.enforce_float(n)

    margin_in: float = property(get_margin_in, set_margin_in)

    def get_alt(self) -> float:
        """Get the conjugaison altitude

        Returns:
            float: Conjugaison altitude (in m)
        """
        return self.__alt

    def set_alt(self, a: float) -> None:
        """Set the conjugaison altitude

        Args:
            a (float): Conjugaison altitude (in m)
        """
        self.__alt = csu.enforce_float(a)

    alt: float = property(get_alt, set_alt)

    def get_thresh(self) -> float:
        """Get the threshold on response for selection

        Returns:
            float: Threshold on response for selection (<1)
        """
        return self.__thresh

    def set_thresh(self, t: float) -> None:
        """Set the threshold on response for selection

        Args:
            t (float): Threshold on response for selection (<1)
        """
        self.__thresh = csu.enforce_float(t)

    thresh: float = property(get_thresh, set_thresh)

    def get_keep_all_actu(self) -> bool:
        """Get the flag for keeping all actuators

        Returns:
            bool: keep all actuator flag
        """
        return self.__keep_all_actu

    def set_keep_all_actu(self, k: bool) -> None:
        """Set the flag for keeping all actuators

        Args:
            k (bool): keep all actuator flag
        """
        self.__keep_all_actu = csu.enforce_or_cast_bool(k)

    keep_all_actu: bool = property(get_keep_all_actu, set_keep_all_actu)

    def get_coupling(self) -> float:
        """Get the actuators coupling

        Returns:
            float: actuators coupling
        """
        return self.__coupling

    def set_coupling(self, c: float) -> None:
        """Set the actuators coupling

        Args:
            c (float): actuators coupling
        """
        self.__coupling = csu.enforce_float(c)

    coupling: float = property(get_coupling, set_coupling)

    def get_unitpervolt(self) -> float:
        """Get the Influence function sensitivity

        Returns:
            float: Influence function sensitivity in unit/volt
        """
        return self.__unitpervolt

    def set_unitpervolt(self, u: float) -> None:
        """Set the Influence function sensitivity

        Args:
            u (float): Influence function sensitivity in unit/volt
        """
        self.__unitpervolt = csu.enforce_float(u)

    unitpervolt: float = property(get_unitpervolt, set_unitpervolt)

    def get_push4imat(self) -> float:
        """Get the nominal voltage for imat

        Returns:
            float: nominal voltage for imat
        """
        return self.__push4imat

    def set_push4imat(self, p: float) -> None:
        """Set the nominal voltage for imat

        Args:
            p (float): nominal voltage for imat
        """
        self.__push4imat = csu.enforce_float(p)

    push4imat: float = property(get_push4imat, set_push4imat)

    def get_ntotact(self) -> int:
        """Get the total number of actuators

        Returns:
            int: Total number of actuators
        """
        return self.__ntotact

    def set_ntotact(self, n: int) -> None:
        """Set the total number of actuators

        Args:
            n (int): Total number of actuators
        """
        self.__ntotact = csu.enforce_int(n)

    _ntotact: int = property(get_ntotact, set_ntotact)

    def get_pitch(self) -> float:
        """Get the actuators pitch [pixels]

        Returns:
            float: Actuators pitch [pixels]
        """
        return self.__pitch

    def set_pitch(self, p: float) -> None:
        """Set the actuators pitch [pixels]

        Args:
            p (float): Actuators pitch [pixels]
        """
        self.__pitch = csu.enforce_float(p)

    _pitch: float = property(get_pitch, set_pitch)

    def get_influsize(self) -> int:
        """Get the actuators influsize [pixels]

        Returns:
            int: Actuators influsize [pixels]
        """
        return self.__influsize

    def set_influsize(self, s: int) -> None:
        """Set the actuators influsize [pixels]

        Args:
            s (int): Actuators influsize [pixels]
        """
        self.__influsize = csu.enforce_int(s)

    _influsize: int = property(get_influsize, set_influsize)

    def get_n1(self) -> int:
        """Get the position of bottom left pixel in the largest support

        Returns:
            int: Position of bottom left pixel in the largest support
        """
        return self.__n1

    def set_n1(self, n: int) -> None:
        """Set the position of bottom left pixel in the largest support

        Args:
            n (int): Position of bottom left pixel in the largest support
        """
        self.__n1 = csu.enforce_int(n)

    _n1: int = property(get_n1, set_n1)

    def get_n2(self) -> int:
        """Get the position of bottom right pixel in the largest support

        Returns:
            int: Position of bottom right pixel in the largest support
        """
        return self.__n2

    def set_n2(self, n: int) -> None:
        """Set the position of bottom right pixel in the largest support

        Args:
            n (int): Position of bottom right pixel in the largest support
        """
        self.__n2 = csu.enforce_int(n)

    _n2: int = property(get_n2, set_n2)

    def get_xpos(self) -> np.ndarray:
        """Get the x positions of influ functions (lower left corner)

        Returns:
            np.ndarray[ndim=1,dtype=np.float32_t]: x positions of influ functions
        """
        return self.__xpos

    def set_xpos(self, xpos: np.ndarray) -> None:
        """Set the x positions of influ functions (lower left corner)

        Args:
            xpos (np.ndarray[ndim=1,dtype=np.float32_t]): x positions of influ functions
        """
        self.__xpos = csu.enforce_array(xpos, self.__ntotact, dtype=np.float32)

    _xpos: np.ndarray = property(get_xpos, set_xpos)

    def get_ypos(self) -> np.ndarray:
        """Get the y positions of influ functions (lower left corner)

        Returns:
            np.ndarray[ndim=1,dtype=np.float32_t]: y positions of influ functions
        """
        return self.__ypos

    def set_ypos(self, ypos: np.ndarray) -> None:
        """Set the y positions of influ functions (lower left corner)

        Args:
            ypos (np.ndarray[ndim=1,dtype=np.float32_t]): y positions of influ functions
        """
        self.__ypos = csu.enforce_array(ypos, self.__ntotact, dtype=np.float32)

    _ypos: np.ndarray = property(get_ypos, set_ypos)

    def get_i1(self) -> np.ndarray:
        """Get the X-position of the bottom left corner of each influence function

        Returns:
            np.ndarray[ndim=1,dtype=np.int32_t]:
        """
        return self.__i1

    def set_i1(self, i1: np.ndarray) -> None:
        """Set the X-position of the bottom left corner of each influence function

        Args:
            i1 (np.ndarray[ndim=1,dtype=np.int32_t]):
        """
        self.__i1 = csu.enforce_array(i1, self.__ntotact, dtype=np.int32)

    _i1: np.ndarray = property(get_i1, set_i1)

    def get_j1(self) -> np.ndarray:
        """Get the Y-position of the bottom left corner of each influence function

        Returns:
            np.ndarray[ndim=1,dtype=np.int32_t]:
        """
        return self.__j1

    def set_j1(self, j1: np.ndarray) -> None:
        """Set the Y-position of the bottom left corner of each influence function

        Args:
            j1 (np.ndarray[ndim=1,dtype=np.int32_t]):
        """
        self.__j1 = csu.enforce_array(j1, self.__ntotact, dtype=np.int32)

    _j1: np.ndarray = property(get_j1, set_j1)

    def get_influ(self) -> np.ndarray:
        """Get the influence function

        Returns:
            np.ndarray[ndim=3,dtype=np.float32_t]: influence function
        """
        return self.__influ

    def set_influ(self, influ: np.ndarray) -> None:
        """Set the influence function

        Args:
            influ (np.ndarray[ndim=3,dtype=np.float32_t]): influence function
        """
        self.__influ = csu.enforce_arrayMultiDim(
            influ, (self.__influsize, self.__influsize, self._ntotact), dtype=np.float32
        )

    _influ: np.ndarray = property(get_influ, set_influ)

    def get_pupoffset(self) -> np.ndarray:
        """Get the pupil offset in meters

        Returns:
            np.ndarray: offsets [m]
        """
        return self.__pupoffset

    def set_pupoffset(self, off: np.ndarray) -> None:
        """Set the pupil offset in meters

        Args:
            off (np.ndarray): offsets [m]
        """
        self.__pupoffset = csu.enforce_array(off, 2, dtype=np.float32)

    pupoffset: np.ndarray = property(get_pupoffset, set_pupoffset)

    def get_puppixoffset(self) -> np.ndarray:
        """Get the pupil offset in pixels

        Returns:
            np.ndarray: offsets [pixels]
        """
        return self.__puppixoffset

    def set_puppixoffset(self, off: np.ndarray) -> None:
        """Set the pupil offset in pixels

        Args:
            off (np.ndarray): offsets [pixels]
        """
        self.__puppixoffset = csu.enforce_array(off, 2, dtype=np.float32)

    _puppixoffset: np.ndarray = property(get_puppixoffset, set_puppixoffset)

    def get_outscl(self) -> float:
        """Get the outer scale for KL with Von Karman spectrum

        Returns:
            float: outer scale [m]
        """
        return self.__outscl

    def set_outscl(self, L0: float) -> None:
        """Set the outer scale for KL with Von Karman spectrum

        Args:
            L0 (float): outer scale [m]
        """
        self.__outscl = csu.enforce_float(L0)

    outscl: float = property(get_outscl, set_outscl)

    def get_nr(self) -> int:
        """Get the number of radial points for KL

        Returns:
            int: Number of radial points
        """
        return self.__nr

    def set_nr(self, n: int) -> None:
        """Set the number of radial points for KL

        Args:
            n (int): Number of radial points
        """
        self.__nr = csu.enforce_int(n)

    _nr: int = property(get_nr, set_nr)

    def get_npp(self) -> int:
        """Get the number of elements (?) for KL

        Returns:
            int: Number of elements
        """
        return self.__npp

    def set_npp(self, n: int) -> None:
        """Set the number of elements (?) for KL

        Args:
            n (int): Number of elements
        """
        self.__npp = csu.enforce_int(n)

    _npp: int = property(get_npp, set_npp)

    def get_ncp(self) -> int:
        """Get the dimension of grid (?)

        Returns:
            int: Dimension
        """
        return self.__ncp

    def set_ncp(self, n: int) -> None:
        """Set the dimension of grid (?)

        Args:
            n (int): Dimension
        """
        self.__ncp = csu.enforce_int(n)

    _ncp: int = property(get_ncp, set_ncp)

    def get_ord(self) -> np.ndarray:
        """Get the radial orders of the basis

        Returns:
            np.ndarray: Radial order of the basis
        """
        return self.__ord

    def set_ord(self, n: np.ndarray) -> None:
        """Set the radial orders of the basis

        Args:
            n (np.ndarray): Radial order of the basis
        """
        self.__ord = csu.enforce_array(n, n.size, dtype=np.int32)

    _ord: np.ndarray = property(get_ord, set_ord)

    def get_rabas(self) -> np.ndarray:
        """Get the radial array of the KL basis

        Returns:
            np.ndarray: Radial array
        """
        return self.__rabas

    def set_rabas(self, r: np.ndarray) -> None:
        """Set the radial array of the KL basis

        Args:
            r (np.ndarray): Radial array
        """
        self.__rabas = csu.enforce_arrayMultiDim(r, r.shape, dtype=np.float32)

    _rabas: np.ndarray = property(get_rabas, set_rabas)

    def get_azbas(self) -> np.ndarray:
        """Get the azimuthal array of the KL basis

        Returns:
            np.ndarray: Azimuthal array
        """
        return self.__azbas

    def set_azbas(self, r: np.ndarray) -> None:
        """Set the azimuthal array of the KL basis

        Args:
            r (np.ndarray): Azimuthal array
        """
        self.__azbas = csu.enforce_arrayMultiDim(r, r.shape, dtype=np.float32)

    _azbas: np.ndarray = property(get_azbas, set_azbas)

    def get_cr(self) -> np.ndarray:
        """Get the radial coordinates in Cartesian grid

        Returns:
            np.ndarray: Radial coordinates in Cartesian grid
        """
        return self.__cr

    def set_cr(self, r: np.ndarray) -> None:
        """Set the radial coordinates in Cartesian grid

        Args:
            r (np.ndarray): Radial coordinates in Cartesian grid
        """
        self.__cr = csu.enforce_arrayMultiDim(r, r.shape, dtype=np.float32)

    _cr: np.ndarray = property(get_cr, set_cr)

    def get_cp(self) -> np.ndarray:
        """Get the phi coordinates in Cartesian grid

        Returns:
            np.ndarray: Phi coordinates in Cartesian grid
        """
        return self.__cp

    def set_cp(self, r: np.ndarray) -> None:
        """Set the phi coordinates in Cartesian grid

        Args:
            r (np.ndarray): Phi coordinates in Cartesian grid
        """
        self.__cp = csu.enforce_arrayMultiDim(r, r.shape, dtype=np.float32)

    _cp: np.ndarray = property(get_cp, set_cp)

    def get_G(self) -> float:
        """Get the magnifying factor

        Returns:
            float: Magnifying factor
        """
        return self.__G

    def set_G(self, G: float) -> None:
        """Set the magnifying factor

        Args:
            G (float): Magnifying factor
        """
        self.__G = csu.enforce_float(G)

    G: float = property(get_G, set_G)

    def get_theta(self) -> float:
        """Get the rotation angle in the pupil.

        Returns:
            float: The rotation angle in radians.
        """
        return self.__theta

    def set_theta(self, theta: float) -> None:
        """Set the rotation angle in the pupil.

        Args:
            theta (float): The rotation angle in radians.
        """
        self.__theta = csu.enforce_float(theta)

    theta: float = property(get_theta, set_theta)

    def get_dx(self) -> float:
        """Get the X axis misalignment.

        Returns:
            float: The X axis misalignment in pixels.
        """
        return self.__dx

    def set_dx(self, dx: float) -> None:
        """Set the X axis misalignment.

        Args:
            dx (float): The X axis misalignment in pixels.
        """
        self.__dx = csu.enforce_float(dx)

    dx: float = property(get_dx, set_dx)

    def get_dy(self) -> float:
        """Get the Y axis misalignment.

        Returns:
            float: The Y axis misalignment in pixels.
        """
        return self.__dy

    def set_dy(self, dy: float) -> None:
        """Set the Y axis misalignment.

        Args:
            dy (float): The Y axis misalignment in pixels.
        """
        self.__dy = csu.enforce_float(dy)

    dy: float = property(get_dy, set_dy)
