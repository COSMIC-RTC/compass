"""
Class SimulatorRTC: COMPASS simulation linked to real RTC with Octopus
"""
import sys
import os

from .simulator import Simulator, init, Iterable, load_config_from_file
from shesha_constants import WFSType, CentroiderType

from CacaoInterfaceWrap import CacaoInterfaceFloat as GFInterface

import numpy as np


class SimulatorRTC(Simulator):
    """
        Class SimulatorRTC: COMPASS simulation linked to real RTC with Octopus
    """

    def __init__(self, filepath: str=None, rtcfilepath: str=None,
                 use_DB: bool=False) -> None:
        """
        Initializes a Simulator instance

        :parameters:
            filepath: (str): (optional) path to the parameters file
            rtcfilepath: (str): (optional) path to the rtc parameters file
            use_DB: (bool): (optional) flag to use dataBase system
        """
        Simulator.__init__(self, filepath, use_DB)
        self.rtcconf = lambda x: None

        if filepath is not None:
            self.load_from_file(filepath, rtcfilepath)

    def load_from_file(self, filepath: str, rtcfilepath: str=None) -> None:
        """
        Load the parameters from the parameters file

        :parameters:
            filepath: (str): path to the parameters file

        """
        Simulator.load_from_file(self, filepath)

        if rtcfilepath is not None:
            load_config_from_file(self.rtcconf, rtcfilepath)

    def init_sim(self) -> None:
        Simulator.init_sim(self)

        wfsNb = len(self.rtcconf.config.p_wfss)
        if wfsNb > 1:
            raise RuntimeError("multi WFS not supported")
        p_wfs = self.rtcconf.config.p_wfss[0]

        nact = self.rtcconf.config.p_dms[0].nact
        framesizex = p_wfs._framesizex
        framesizey = p_wfs._framesizey
        nvalid = p_wfs._nvalid
        if p_wfs.type == WFSType.SH:
            frame = self.wfs.get_binimg(0)
        elif p_wfs.type == WFSType.PYRHR:
            frame = self.wfs.get_pyrimg(0)
        else:
            raise RuntimeError("WFS Type not usable")

        if frame.shape != (framesizex, framesizey):
            raise RuntimeError("framesize not match with the simulation")

        if self.rtc.get_voltage(0).size != nact:
            raise RuntimeError("nact not match with the simulation")

        if self.rtc.get_cmat(0).shape != (nact, nvalid * 2):
            raise RuntimeError("cmat not match with the simulation")

        self.fakewfs = GFInterface(p_wfs._frameShmName)

        self.comp = np.zeros(nact, dtype=np.float32)
        self.fakedms = GFInterface(self.rtcconf.config.p_dms[0]._actuShmName)

        tmp_cmat = self.rtc.get_cmat(0)
        self.cmat = GFInterface(self.rtcconf.config.p_controllers[0]._cmatShmName)
        self.cmat.send(tmp_cmat)

        tmp_valid = self.config.p_wfss[0].get_validsub()
        self.valid = GFInterface(p_wfs._validsubsShmName)
        if p_wfs.type == WFSType.SH:
            self.valid.send(tmp_valid * self.config.p_wfss[0].npix)
        elif p_wfs.type == WFSType.PYRHR:
            self.valid.send(tmp_valid)
        else:
            raise RuntimeError("WFS Type not usable")

    def next(self, *, move_atmos: bool=True, see_atmos: bool=True, nControl: int=0,
             tar_trace: Iterable[int]=None, wfs_trace: Iterable[int]=None,
             do_control: bool=True, apply_control: bool=True) -> None:
        """
        Overload of the Simulator.next() function
        """
        Simulator.next(self, move_atmos=move_atmos, see_atmos=see_atmos,
                       nControl=nControl, tar_trace=[0], wfs_trace=[0], do_control=False)

        # print("Send a frame")
        p_wfs = self.rtcconf.config.p_wfss[0]
        if p_wfs.type == WFSType.SH:
            frame = self.wfs.get_binimg(0)
        elif p_wfs.type == WFSType.PYRHR:
            frame = self.wfs.get_pyrimg(0)
        else:
            raise RuntimeError("WFS Type not usable")
        self.fakewfs.send(frame)
        if apply_control:
            # print("Wait a command...")
            self.fakedms.recv(self.comp, 0)
            self.dms.set_full_comm(self.comp)
