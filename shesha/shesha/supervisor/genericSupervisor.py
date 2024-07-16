## @package   shesha.supervisor.aoSupervisor
## @brief     Abstract layer for initialization and execution of a AO supervisor
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


from abc import ABC, abstractmethod
from shesha.sutra_wrap import carma_context


class GenericSupervisor(ABC):
    """ This class defines generic methods and behavior of a supervisor
    It is not intended to be instantiated as it is : prefer to build
    a supervisor class inheriting from it. This approach allows to build multiple
    supervisors with various components and less effort

    Attributes:
        context : (CarmaContext) : a CarmaContext instance

        config : (config) : Parameters structure

        is_init : (bool) : Flag equals to True if the supervisor has already been initialized

        iter : (int) : Frame counter
    """

    def __init__(self, config):
        """ Init the a supervisor

        Args:
            config : (config module) : Configuration module

        """
        self.context = None
        self.config = config
        self.is_init = False
        self.iter = 0

        if (self.config.p_loop.devices.size > 1):
            self.context = carma_context.get_instance_ngpu(
                    self.config.p_loop.devices.size, self.config.p_loop.devices)
        else:
            self.context = carma_context.get_instance_1gpu(
                    self.config.p_loop.devices[0])
        self.force_context()

        self._init_components()

    def get_config(self):
        """ Returns the configuration in use, in a supervisor specific format ?

        Returns:
            config : (config module) : Current supervisor configuration
        """
        return self.config

    def get_frame_counter(self) -> int:
        """Return the current iteration number of the loop

        Returns:
            framecounter : (int) : Number of iteration already performed
        """
        return self.iter

    def force_context(self) -> None:
        """ Active all the GPU devices specified in the parameters file
        """
        if self.context is not None:
            current_device = self.context.active_device
            for device in range(len(self.config.p_loop.devices)):
                self.context.set_active_device_force(device)
            self.context.set_active_device(current_device)

    @abstractmethod
    def _init_components(self) -> None:
        """ Initialize all the components
        """
        self.is_init = True
