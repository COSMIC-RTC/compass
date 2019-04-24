"""

Class SimulatorBrahma: Brahma overloaded simulator

"""

import sys
import os

from .simulator import Simulator
from shesha.init.rtc_init import rtc_init
from shesha.init.target_init import target_init


class SimulatorBrahma(Simulator):
    """
        Class SimulatorBrahma: Brahma overloaded simulator
        _tar_init and _rtc_init to instantiate Brahma classes instead of regular classes
        next() to call rtc/tar.publish()
    """

    def _tar_init(self) -> None:
        '''
        Initializes a Target_Brahma object
        '''
        if self.config.p_targets is not None:
            print("->target")
            self.tar = target_init(self.c, self.tel, self.config.p_targets,
                                   self.config.p_atmos, self.config.p_tel,
                                   self.config.p_geom, self.config.p_dms, brahma=True)
        else:
            self.tar = None

    def _rtc_init(self, ittime) -> None:
        '''
        Initializes a Rtc_Brahma object
        '''
        if self.config.p_controllers is not None or self.config.p_centroiders is not None:
            print("->rtc")
            #   rtc
            self.rtc = rtc_init(self.c, self.tel, self.wfs, self.dms, self.atm,
                                self.config.p_wfss, self.config.p_tel,
                                self.config.p_geom, self.config.p_atmos, ittime,
                                self.config.p_centroiders, self.config.p_controllers,
                                self.config.p_dms, brahma=True, tar=None)
        else:
            self.rtc = None

    def next(self, **kwargs) -> None:
        """
        Overload of the Simulator.next() function with BRAHMA publications
        """
        Simulator.next(self, **kwargs)
        if self.rtc is not None:
            self.rtc.publish()
        if self.tar is not None:
            self.tar.publish()
