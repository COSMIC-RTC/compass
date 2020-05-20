## @package   shesha.supervisor.aoSupervisor
## @brief     Abstract layer for initialization and execution of a AO supervisor
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

from abc import abstractmethod
import numpy as np
from shesha.util.utilities import load_config_from_file



class GenericSupervisor(object):
    """ This class defines generic methods and behavior of a supervisor
    It is not intended to be instantiated as it is : prefer to build
    a supervisor class inheriting from it. This approach allows to build multiple
    supervisors with various components and less effort

    Attributes:
        context : (CarmaContext) : a CarmaContext instance

        config : (config) : Parameters structure

        telescope : (TelescopeComponent) : a TelescopeComponent instance

        atmos : (AtmosComponent) : An AtmosComponent instance

        target : (targetComponent) : A TargetComponent instance

        wfs : (WfsComponent) : A WfsComponent instance

        dms : (DmComponent) : A DmComponent instance

        rtc : (RtcComponent) : A Rtc component instance

        is_loaded : (bool) : Flag equals to True if a config has already been loaded

        is_init : (bool) : Flag equals to True if the supervisor has already been initialized

        iter : (int) : Frame counter     
    """
    def __init__(self, config_file : str = None):
        """ Init the a supervisor

        Parameters:
            config_file: (str, optional) : Path to a configuration file.
                                           If provided, loads it directly
        """
        self.config = None
        self.telescope = None
        self.atmos = None
        self.target = None
        self.wfs = None
        self.dms = None
        self.rtc = None
        self.is_loaded = False
        self.is_init = False
        if config_file is not None:
            self.load_config(config_file)
    
    def load_config(self, config_file: str) -> None:
        """ Load the parameters configuration from config_file

        Parameters:
            config_file : (str) : path to the configuration file
        """
        self.is_loaded = False
        self.is_init = False
        self.config = load_config_from_file(config_file)
        self.is_loaded = True
    
    def get_config(self):
        """ Returns the configuration in use, in a supervisor specific format ?

        Return:
            config : (config module) : Current supervisor configuration
        """
        return self.config

    def get_frame_counter(self) -> int:
        """Return the current iteration number of the loop

        Return:
            framecounter : (int) : Number of iteration already performed
        """
        return self.iter
    
    def init(self) -> None:
        """ Initialize all the components
        """
        if self.is_loaded:
            if self.config.p_tel is not None:
                self.init_tel()
            if self.config.p_atmos is not None:
                self.init_atmos()
            if self.config.p_dms is not None:
                self.init_dms()
            if self.config.p_targets is not None:
                self.init_target()
            if self.config.p_wfss is not None:
                self.init_wfs()
            if self.config.p_controllers is not None or self.config.p_centroiders is not None:
                self.init_rtc()
        self.is_init = True
    
    @abstractmethod
    def init_tel(self):
        """ Initialize the telescope component of the supervisor
        """
        pass

    @abstractmethod
    def init_atmos(self):
        """ Initialize the atmos component of the supervisor
        """
        pass

    @abstractmethod
    def init_dms(self):
        """ Initialize the dms component of the supervisor
        """
        pass

    @abstractmethod
    def init_target(self):
        """ Initialize the target component of the supervisor
        """
        pass

    @abstractmethod
    def init_wfs(self):
        """ Initialize the wfs component of the supervisor
        """
        pass

    @abstractmethod
    def init_rtc(self):
        """ Initialize the rtc component of the supervisor
        """
        pass

    def next(self, *, move_atmos: bool = True, nControl: int = 0,
             tar_trace: Iterable[int] = None, wfs_trace: Iterable[int] = None,
             do_control: bool = True, apply_control: bool = True,
             compute_tar_psf: bool = True) -> None:
        """Iterates the AO loop, with optional parameters

        Parameters :
             move_atmos: (bool): move the atmosphere for this iteration, default: True

             nControl: (int): Controller number to use, default 0 (single control configurations)

             tar_trace: (None or list[int]): list of targets to trace. None equivalent to all.

             wfs_trace: (None or list[int]): list of WFS to trace. None equivalent to all.

             apply_control: (bool): (optional) if True (default), apply control on DMs
        """
        if tar_trace is None and self.target is not None:
            tar_trace = range(len(self.config.p_targets))
        if wfs_trace is None and self.wfs is not None:
            wfs_trace = range(len(self.config.p_wfss))

        if move_atmos and self.atmos is not None:
            self.atmos.move_atmos()

        if tar_trace is not None:
            for t in tar_trace:
                if self.atmos.is_enable:
                    self.target.raytrace(t, tel=self.tel, atm=self.atmos, dm=self.dms)
                else:
                    self.target.raytrace(t, tel=self.tel, dm=self.dms)

        if wfs_trace is not None:
            for w in wfs_trace:
                if self.atmos.is_enable:
                    self.wfs.raytrace(w, tel=self.tel, atm=self.atmos)
                else:
                    self.wfs.raytrace(w, tel=self.tel)

                if not self.config.p_wfss[w].open_loop and self.dms is not None:
                    self.wfs.raytrace(w, dm=self.dms, reset=False)
                self.wfs.compute_wfs_image(w)
        if do_control and self.rtc is not None:
            for ncontrol in range(len(self.config.p_controllers)):
                self.rtc.do_centroids(ncontrol)
                self.rtc.do_control(ncontrol)
                self.rtc.do_clipping(ncontrol)

        if apply_control:
            self.rtc.apply_control(ncontrol)

        if compute_tar_psf:
            for tar_index in tar_trace:
                self.target.comp_tar_image(tar_index)
                self.target.comp_strehl(tar_index)

        self.iter += 1
