## @package   shesha.supervisor.abstractSupervisor
## @brief     Abstract layer for initialization and execution of a COMPASS supervisor
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   5.0.0
## @date      2020/05/18
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
#  All rights reserved.
#  Distributed under GNU - LGPL
#
#  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
#  General Public License as published by the Free Software Foundation, either version 3 of the License,
#  or any later version.
#
#  COMPASS: End-to-end AO simulation tool using GPU acceleration
#  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
#
#  The final product includes a software package for simulating all the critical subcomponents of AO,
#  particularly in the context of the ELT and a real-time core based on several control approaches,
#  with performances consistent with its integration into an instrument. Taking advantage of the specific
#  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
#  conduct large simulation campaigns called to the ELT.
#
#  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
#  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
#  various systems configurations such as multi-conjugate AO.
#
#  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
#  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
#  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.

from abc import ABC, abstractmethod
import numpy as np
from tqdm import trange


class AbstractSupervisor(ABC):

    @abstractmethod
    def get_config(self):
        """ Returns the configuration in use, in a supervisor specific format ?
        """
        ...

    @abstractmethod
    def load_config(self):
        """ Load the configuration for the supervisor
        """
        ...

    @abstractmethod
    def init_config(self):
        """ Init the configuration for the supervisor
        """
        ...

    @abstractmethod
    def set_command(self, command: np.ndarray) -> None:
        """ Immediately sets provided command to DMs - does not affect integrator

        Parameters:
            command : (np.ndarray) : DM Command vector
        """
        ...

    @abstractmethod
    def set_perturbation_voltage(self, nControl: int, name: str,
                                 command: np.ndarray) -> None:
        """ Add circular buffer of offset values to integrator (will be applied at the end of next iteration)

        Parameters:
            nControl : (int) : Controller index

            name : (str) : Buffer name

            command : (np.ndarray) : perturbation voltage circular buffer
        """
        ...

    @abstractmethod
    def get_slopes(self) -> np.ndarray:
        """ Immediately gets one slope vector for all WFS at the current state of the system

        Return:
            slopes : (np.ndarray) : Current WFS slopes
        """
        ...

    @abstractmethod
    def single_next(self, move_atmos: bool = True, show_atmos: bool = True,
                   get_tar_image: bool = False, get_residual: bool = False) -> None:
        """ Performs a single loop iteration 

        Parameters:
            move_atmos : (bool, optional) : Move the atmosphere layers. Default is True

            show_atmos : (bool, optional) : WFS and targets see the atmosphere layers. Default is True

            get_tar_image : (bool, optional) : 

            get_residual : (bool, optional) :
        """
        ...

    def next(self, nbiters, see_atmos=True):
        """ Performs nbiters loop iterations

        Parameters:
            nbiters : (int) : Number of iterations to perform

            see_atmos : (bool, optional) : Flag to enable atmosphere. Default is True
        """
        for i in range(nbiters):
            self.single_next(show_atmos=see_atmos)

    @abstractmethod
    def close_loop(self) -> None:
        """ DM receives controller output + pertuVoltage
        """
        ...

    @abstractmethod
    def open_loop(self) -> None:
        """ Integrator computation goes to /dev/null but pertuVoltage still applied
        """
        ...

    @abstractmethod
    def set_ref_slopes(self, ref_slopes: np.ndarray) -> None:
        """ Set given ref slopes in controller
        """
        ...

    @abstractmethod
    def get_ref_slopes(self) -> np.ndarray:
        """ Get the currently used reference slopes
        """
        ...

    @abstractmethod
    def set_gain(self, gain: float):
        """ Set the scalar gain of feedback controller loop
        """
        ...

    @abstractmethod
    def set_command_matrix(self, cMat: np.ndarray):
        """ Set the cmat for the controller to use
        """
        ...

    @abstractmethod
    def get_tar_image(self, tar_index):
        """ Get an image from a target 
        """
        ...

    @abstractmethod
    def get_wfs_image(self, numWFS: int = 0, cal_pix: bool = False) -> np.ndarray:
        """ Get an image from the WFS. If calpix = True returns the calibrated image (background + flat + pixels selection )
        """
        ...

    @abstractmethod
    def get_intensities(self):
        """ Return sum of intensities in subaps. Size nSubaps, same order as slopes
        """
        ...

    @abstractmethod
    def get_frame_counter(self):
        """ Return the current frame counter
        """
        ...
