"""
Class SimulatorCacao: Cacao overloaded simulator
"""
import sys
import os

from .simulator import Simulator
from shesha.init.rtc_init import rtc_init


class SimulatorCacao(Simulator):
    """
        Class SimulatorCacao: Cacao overloaded simulator
        _tar_init and _rtc_init to instantiate Cacao classes instead of regular classes
        next() to call rtc/tar.publish()
    """

    def _rtc_init(self, ittime) -> None:
        '''
        Initializes a Rtc_Cacao object
        '''
        if self.config.p_controllers is not None or self.config.p_centroiders is not None:
            print("->rtc")
            #   rtc
            self.rtc = rtc_init(self.c, self.tel, self.wfs, self.dms, self.atm,
                                self.config.p_wfss, self.config.p_tel,
                                self.config.p_geom, self.config.p_atmos, ittime,
                                self.config.p_centroiders, self.config.p_controllers,
                                self.config.p_dms, cacao=True, tar=None)
        else:
            self.rtc = None

    def next(self, **kwargs) -> None:
        """
        Overload of the Simulator.next() function with cacao publications
        """
        self.rtc.d_centro[0].load_img(self.wfs.d_wfs[0].d_binimg)
        Simulator.next(self, **kwargs)
        if self.rtc is not None:
            self.rtc.publish()
