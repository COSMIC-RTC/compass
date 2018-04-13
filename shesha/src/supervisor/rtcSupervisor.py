from .benchSupervisor import BenchSupervisor, naga_context, rtc_standalone
from shesha_sim.simulator import load_config_from_file

import numpy as np

from CacaoInterfaceWrap import CacaoInterfaceFloat as GFInterface

from shesha_constants import WFSType, CentroiderType


class RTCSupervisor(BenchSupervisor):

    def __init__(self, configFile: str=None, BRAHMA: bool=False):
        '''
        Init the COMPASS wih the configFile
        '''
        self._sim = lambda: None
        self.cam = None
        self.rtc = None
        self.npix = 0
        self.BRAHMA = BRAHMA

        self._sim.is_init = False  # type: bool
        self._sim.loaded = False  # type: bool
        self._sim.config = None  # type: Any # types.ModuleType ?

        if configFile is not None:
            self.loadConfig(configFile)

    def clearInitSim(self) -> None:
        """
        Delete objects initialized in a previous simulation
        """
        if self._sim.loaded and self._sim.is_init:
            self._sim.iter = 0

            del self._sim.atm
            self._sim.atm = None
            del self._sim.tel
            self._sim.tel = None
            del self._sim.tar
            self._sim.tar = None
            del self._sim.rtc
            self._sim.rtc = None
            del self._sim.wfs
            self._sim.wfs = None
            del self._sim.dms
            self._sim.dms = None

            # del self._sim.c  # What is this supposed to do ... ?
            # self._sim.c = None

        self._sim.is_init = False

    def forceContext(self) -> None:
        '''
        Clear the initialization of the simulation
        '''
        ...

    def singleNext(self, moveAtmos: bool=True, showAtmos: bool=True, getPSF: bool=False,
                   getResidual: bool=False) -> None:
        '''
        Move atmos -> getSlope -> applyControl ; One integrator step
        '''
        frame_size = int(np.sqrt(self.fakewfs.size))
        frame = np.zeros((frame_size, frame_size), dtype=np.float32)
        # print("Wait a frame...")
        self.fakewfs.recv(frame, 0)
        self.rtc.load_rtc_img(0, frame.astype(np.float32))
        if self._sim.config.p_wfss[0].type == WFSType.SH:
            #for SH
            self.rtc.fill_rtc_bincube(0, self.npix)
        self.rtc.do_centroids(0)
        self.rtc.do_control(0)
        self.rtc.save_com(0)
        # print("Send a command")
        comms = self.rtc.get_com(0)
        self.fakedms.send(comms)

    def loadConfig(self, configFile: str) -> None:
        '''
        Init the COMPASS wih the configFile
        '''
        load_config_from_file(self._sim, configFile)

    def initConfig(self) -> None:
        '''
        Initialize the simulation
        '''
        print("->cam")
        self.fakewfs = GFInterface("compass_wfs")
        self.fakedms = GFInterface("compass_dms")
        nact = self.fakedms.size

        self.cmat = GFInterface("compass_cmat")
        cMat_data = np.zeros(self.cmat.size, dtype=np.float32)
        self.cmat.recv(cMat_data)
        nvalid = self.cmat.size // nact // 2

        self.valid = GFInterface("compass_valid")
        tmp_valid = np.zeros(self.valid.size, dtype=np.float32)
        self.valid.recv(tmp_valid)
        self._sim.config.p_nvalid = np.reshape(tmp_valid, (2, nvalid))

        print("->RTC")
        wfsNb = len(self._sim.config.p_wfss)
        if wfsNb > 1:
            raise RuntimeError("multi WFS not supported")

        if self._sim.config.p_wfss[0].type == WFSType.SH:
            self.npix = self._sim.config.p_wfss[0].npix

            if "p_nvalid" not in self._sim.config.__dict__.keys(
            ) or self._sim.config.p_nvalid is None:
                import rtcData.DataInit as di
                dataS = di.makeSH(wfsNb=wfsNb, frameSize=self.fakewfs.data.md.size[0],
                                  roiSize=self._sim.config.p_wfss[0].nxsub,
                                  subSize=self.npix)
                xvalid = dataS.data["roiTab"].data[0, :] / self.npix
                yvalid = dataS.data["roiTab"].data[1, :] / self.npix
            else:
                xvalid = self._sim.config.p_nvalid[0, :] / self.npix
                yvalid = self._sim.config.p_nvalid[1, :] / self.npix

            self.context = naga_context(devices=np.array([0], dtype=np.int32))
            print("nvalid : %d" % nvalid)
            offset = self.npix // 2 + 0.5
            scale = 0.29005988378497927
            gain = .4
            print("nact : %d" % nact)

            self.rtc = rtc_standalone(self.context, wfsNb, [nvalid], nact,
                                      self._sim.config.p_centroiders[0].type, 1, offset,
                                      scale)
            self.rtc.set_decayFactor(0, np.ones(nact, dtype=np.float32))
            self.rtc.set_matE(0, np.identity(nact, dtype=np.float32))
            self.rtc.set_mgain(0, np.ones(nact, dtype=np.float32) * gain)
            self.rtc.set_cmat(0, np.reshape(cMat_data, (nact, cMat_data.size // nact)))
            self.rtc.load_rtc_validpos(0,
                                       xvalid.astype(np.int32), yvalid.astype(np.int32))

        elif self._sim.config.p_wfss[0].type == WFSType.PYRHR:
            ...
        self._sim.is_init = True

    def getRawWFSImage(self, numWFS: int=0) -> np.ndarray:
        '''
        Get an image from the WFS
        '''
        return np.array(self.fakewfs.data)
