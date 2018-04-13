"""
Class SimulatorRTC: COMPASS simulation linked to real RTC with Octopus
"""
import sys
import os

from .simulator import Simulator, init, Iterable

from CacaoInterfaceWrap import CacaoInterfaceFloat as GFInterface

import numpy as np


class SimulatorRTC(Simulator):
    """
        Class SimulatorRTC: COMPASS simulation linked to real RTC with Octopus
    """

    def init_sim(self) -> None:
        super().init_sim()

        self.fakewfs = GFInterface("compass_wfs", self.wfs.get_binimg(0).shape, 1)
        self.nact = self.rtc.get_voltage(0).size
        self.fakedms = GFInterface("compass_dms", (1, self.rtc.get_voltage(0).size), 1)
        tmp_cmat = self.rtc.get_cmat(0)
        self.cmat = GFInterface("compass_cmat", tmp_cmat.shape, 1)
        self.cmat.send(tmp_cmat)
        tmp_valid = self.config.p_wfss[0].get_validsub()
        self.valid = GFInterface("compass_valid", tmp_valid.shape, 1)
        self.valid.send(tmp_valid * self.config.p_wfss[0].npix)

    def next(self, *, move_atmos: bool=True, see_atmos: bool=True, nControl: int=0,
             tar_trace: Iterable[int]=None, wfs_trace: Iterable[int]=None,
             do_control: bool=True, apply_control: bool=True) -> None:
        """
        Overload of the Simulator.next() function
        """
        Simulator.next(self, move_atmos=move_atmos, see_atmos=see_atmos,
                       nControl=nControl, tar_trace=[0], wfs_trace=[0], do_control=False)

        # print("Send a frame")
        self.fakewfs.send(self.wfs.get_binimg(0))
        if apply_control:
            comp = np.zeros(self.nact, dtype=np.float32)
            # print("Wait a command...")
            self.fakedms.recv(comp, 0)
            self.dms.set_full_comm(comp)
