"""Widget to simulate a closed loop

Usage:
  canapassSupervisor.py [<parameters_filename>]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
"""

import os, sys
import numpy as np
import time
from collections import OrderedDict

from tqdm import trange
import astropy.io.fits as pfits
from threading import Thread
from subprocess import Popen, PIPE

import shesha.ao as ao
import shesha.constants as scons

from typing import Any, Dict, Tuple, Callable, List
from .compassSupervisor import CompassSupervisor

from naga.obj import obj_Double2D
from naga.magma import syevd_Double, svd_host_Double
from naga.context import context as naga_context

from naga.host_obj import host_obj_Double1D, host_obj_Double2D


class CanapassSupervisor(CompassSupervisor):

    def __init__(self, configFile: str=None, BRAHMA: bool=True) -> None:
        CompassSupervisor.__init__(self, configFile, True)

        #############################################################
        #                 CONNECTED BUTTONS                         #
        #############################################################
        # Default path for config files

        self.ph2modes = None
        self.KL2V = None
        self.P = None
        self.currentBuffer = 1
        #############################################################
        #                       METHODS                             #
        #############################################################

    """
          ____ ___  __  __ ____   _    ____ ____
         / ___/ _ \|  \/  |  _ \ / \  / ___/ ___|
        | |  | | | | |\/| | |_) / _ \ \___ \___ \
        | |__| |_| | |  | |  __/ ___ \ ___) |__) |
         \____\___/|_|  |_|_| /_/   \_\____/____/
         ____  _   _ ____  _____ ______     _____ ____   ___  ____
        / ___|| | | |  _ \| ____|  _ \ \   / /_ _/ ___| / _ \|  _ \
        \___ \| | | | |_) |  _| | |_) \ \ / / | |\___ \| | | | |_) |
         ___) | |_| |  __/| |___|  _ < \ V /  | | ___) | |_| |  _ <
        |____/ \___/|_|   |_____|_| \_\ \_/  |___|____/ \___/|_| \_\
         __  __ _____ _____ _   _  ___  ____  ____
        |  \/  | ____|_   _| | | |/ _ \|  _ \/ ___|
        | |\/| |  _|   | | | |_| | | | | | | \___ \
        | |  | | |___  | | |  _  | |_| | |_| |___) |
        |_|  |_|_____| |_| |_| |_|\___/|____/|____/
    """

    def getConfig(self, path=None):
        ''' Returns the configuration in use, in a supervisor specific format '''
        if path:
            self.writeConfigOnFile(path)
            return
        return self._sim.config

    def loadConfig(self, configFile: str, BRAMA: bool=True) -> None:
        ''' Load the configuration for the compass supervisor'''
        CompassSupervisor.loadConfig(self, configFile, BRAMA)
        print("switching to a generic controller")
        self._sim.config.p_controllers[0].type = scons.ControllerType.GENERIC

    """
          ____    _    _   _    _    ____   _    ____ ____
         / ___|  / \  | \ | |  / \  |  _ \ / \  / ___/ ___|
        | |     / _ \ |  \| | / _ \ | |_) / _ \ \___ \___ \
        | |___ / ___ \| |\  |/ ___ \|  __/ ___ \ ___) |__) |
         \____/_/   \_\_| \_/_/   \_\_| /_/   \_\____/____/
         __  __ _____ _____ _   _  ___  ____  ____
        |  \/  | ____|_   _| | | |/ _ \|  _ \/ ___|
        | |\/| |  _|   | | | |_| | | | | | | \___ \
        | |  | | |___  | | |  _  | |_| | |_| |___) |
        |_|  |_|_____| |_| |_| |_|\___/|____/|____/
    """

    def computePh2Modes(self):
        oldnoise = self._sim.config.p_wfs0.noise
        self.setNoise(-1)

        if (self.modalBasis is None):
            self.modalBasis, self.P = self.getModes2VBasis("Btt")
        nbmode = self.modalBasis.shape[1]
        pup = self._sim.config.p_geom._spupil
        ph = np.array(self._sim.tar.d_targets[0].d_phase)
        ph2KL = np.zeros((nbmode, ph.shape[0], ph.shape[1]))
        S = np.sum(pup)
        for i in trange(nbmode):
            self._sim.tar.d_targets[0].reset_phase()
            self._sim.dms.set_full_comm(
                    (self.modalBasis[:, i]).astype(np.float32).copy())
            self._sim.next(see_atmos=False)
            ph = np.array(self._sim.tar.d_targets[0].d_phase) * pup
            # Normalisation pour les unites rms en microns !!!
            norm = np.sqrt(np.sum((ph)**2) / S)
            ph2KL[i] = ph / norm
        self.ph2modes = ph2KL
        self.setNoise(oldnoise)
        return ph2KL

    def next(self, nbiters, see_atmos=True):
        for i in trange(nbiters):
            self._sim.next(see_atmos=see_atmos)

    def computePh2ModesFits(self, fullpath):
        data = self.computePh2Modes()
        self.writeDataInFits(data, fullpath)

    def getModes2VBasis(self, ModalBasisType):
        if (ModalBasisType == "KL2V"):
            print("Computing KL2V basis...")
            self.modalBasis, _ = self.returnkl2V()
            return self.modalBasis, 0
        elif (ModalBasisType == "Btt"):
            print("Computing Btt basis...")
            self.modalBasis, self.P = self.compute_Btt2(inv_method="cpu_svd")
            return self.modalBasis, self.P

    def returnkl2V(self):
        """
        KL2V = ao.compute_KL2V(self._sim.config.p_controllers[
                               0], self._sim.dms, self._sim.config.p_dms, self._sim.config.p_geom, self._sim.config.p_atmos, self._sim.config.p_tel)
        """
        if (self.KL2V is None):
            print("Computing KL2V...")
            KL2V = ao.compute_KL2V(self._sim.config.p_controllers[0], self._sim.dms,
                                   self._sim.config.p_dms, self._sim.config.p_geom,
                                   self._sim.config.p_atmos, self._sim.config.p_tel)
            print("KL2V Done!")
            self.KL2V = KL2V
            return KL2V, 0
        else:
            return self.KL2V, 0

    def setGain(self, gain: float) -> None:
        CompassSupervisor.setGain(self, gain)

    # TODO : use ao.computeBtt instead, use compute_DMbasis to get IF sparse if pzt or IF dense if TT
    def compute_Btt2(self, inv_method: str="gpu_evd"):
        IF = self._sim.rtc.d_control[1].d_IFsparse.get_csr()
        n = IF.shape[0]
        N = IF.shape[1]
        T = np.array(self._sim.rtc.d_control[1].d_TT)

        delta = IF.dot(IF.T).toarray() / N

        # Tip-tilt + piston
        Tp = np.ones((T.shape[0], T.shape[1] + 1))
        Tp[:, :2] = T.copy()  #.toarray()
        deltaT = IF.dot(Tp) / N
        # Tip tilt projection on the pzt dm
        tau = np.linalg.inv(delta).dot(deltaT)

        # Famille generatrice sans tip tilt
        G = np.identity(n)
        tdt = tau.T.dot(delta).dot(tau)
        subTT = tau.dot(np.linalg.inv(tdt)).dot(tau.T).dot(delta)
        G -= subTT

        # Base orthonormee sans TT
        gdg = G.T.dot(delta).dot(G)

        startTimer = time.time()
        if inv_method == "cpu_svd":
            print("Doing SVD on CPU of a matrix...")
            U, s, _ = np.linalg.svd(gdg)
        elif inv_method == "gpu_svd":
            print("Doing SVD on CPU of a matrix...")
            m = gdg.shape[0]
            h_mat = host_obj_Double2D(data=gdg, mallocType="pagelock")
            h_eig = host_obj_Double1D(data=np.zeros([m], dtype=np.float64),
                                      mallocType="pagelock")
            h_U = host_obj_Double2D(data=np.zeros((m, m), dtype=np.float64),
                                    mallocType="pagelock")
            h_VT = host_obj_Double2D(data=np.zeros((m, m), dtype=np.float64),
                                     mallocType="pagelock")
            svd_host_Double(h_mat, h_eig, h_U, h_VT)
            U = h_U.getData().T.copy()
            s = h_eig.getData()[::-1].copy()
        elif inv_method == "gpu_evd":
            print("Doing EVD on GPU of a matrix...")
            c = naga_context()
            m = gdg.shape[0]
            d_mat = obj_Double2D(c, data=gdg)
            d_U = obj_Double2D(c, data=np.zeros([m, m], dtype=np.float64))
            h_s = np.zeros(m, dtype=np.float64)
            syevd_Double(d_mat, h_s, d_U)
            U = d_U.device2host().T.copy()
            s = h_s[::-1].copy()
        else:
            raise "ERROR cannot recognize inv_method"
        print("Done in %fs" % (time.time() - startTimer))
        U = U[:, :U.shape[1] - 3]
        s = s[:s.size - 3]
        L = np.identity(s.size) / np.sqrt(s)
        B = G.dot(U).dot(L)

        # Rajout du TT
        TT = T.T.dot(T) / N  #.toarray()/N
        Btt = np.zeros((n + 2, n - 1))
        Btt[:B.shape[0], :B.shape[1]] = B
        mini = 1. / np.sqrt(np.abs(TT))
        mini[0, 1] = 0
        mini[1, 0] = 0
        Btt[n:, n - 3:] = mini

        # Calcul du projecteur actus-->modes
        delta = np.zeros((n + T.shape[1], n + T.shape[1]))
        delta[:-2, :-2] = IF.dot(IF.T).toarray() / N
        delta[-2:, -2:] = T.T.dot(T) / N
        P = Btt.T.dot(delta)

        return Btt.astype(np.float32), P.astype(np.float32)

    def doImatModal(self, ampliVec, KL2V, Nslopes, noise=False, nmodesMax=0,
                    withTurbu=False, pushPull=False):
        iMatKL = np.zeros((KL2V.shape[1], Nslopes))
        #currentVolts = self._sim.rtc.get_voltage(0).copy()[None,:]

        if (nmodesMax):
            KLMax = nmodesMax
        else:
            KLMax = KL2V.shape[1]
        for kl in trange(KLMax, desc="Modal IM"):
            v = ampliVec[kl] * KL2V[:, kl:kl + 1].T.copy()
            if ((pushPull is True) or
                (withTurbu is True)):  # with turbulence/aberrations => push/pull
                self._sim.rtc.d_control[0].set_perturbcom(
                        v)  # Adding Perturbation voltage on current iteration
                devpos = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                self._sim.rtc.d_control[0].set_perturbcom(-v)
                devmin = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                iMatKL[kl, :] = (devpos - devmin) / (2. * ampliVec[kl])
                #imat[:-2, :] /= pushDMMic
                #if(nmodesMax == 0):# i.e we measured all modes including TT
                #imat[-2:, :] /= pushTTArcsec
            else:  # No turbulence => push only
                self._sim.rtc.d_control[0].set_openloop(1)  # openLoop
                self._sim.rtc.d_control[0].set_perturbcom(v)
                iMatKL[kl, :] = self.applyVoltGetSlopes(noise=noise) / ampliVec[kl]

        # print("Modal interaction matrix done in %3.0f seconds" % (time.time() - st))

        return iMatKL

    def doImatPhase(self, cubePhase, Nslopes, noise=False, nmodesMax=0, withTurbu=False,
                    pushPull=False, wfsnum=0):
        iMatPhase = np.zeros((cubePhase.shape[0], Nslopes))
        for nphase in trange(cubePhase.shape[0], desc="Phase IM"):
            if ((pushPull is True) or
                (withTurbu is True)):  # with turbulence/aberrations => push/pull
                self.setNcpaWfs(cubePhase[nphase, :, :], wfsnum=wfsnum)
                devpos = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                self.setNcpaWfs(-cubePhase[nphase, :, :], wfsnum=wfsnum)
                devmin = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
                iMatPhase[nphase, :] = (devpos - devmin) / 2
            else:  # No turbulence => push only
                self._sim.rtc.d_control[0].set_openloop(1)  # openLoop
                self.setNcpaWfs(cubePhase[nphase, :, :], wfsnum=wfsnum)
                iMatPhase[nphase, :] = self.applyVoltGetSlopes(noise=noise)
        self.setNcpaWfs(cubePhase[nphase, :, :] * 0.,
                        wfsnum=wfsnum)  # Remove the Phase on WFS
        _ = self.applyVoltGetSlopes(turbu=withTurbu, noise=noise)
        # print("Phase interaction matrix done in %3.0f seconds" % (time.time() - st))

        return iMatPhase

    def applyVoltGetSlopes(self, noise=False, turbu=False, reset=1):
        self._sim.rtc.apply_control(0, self._sim.dms)
        for w in range(len(self._sim.config.p_wfss)):
            if (turbu):
                self._sim.wfs.d_wfs[w].d_gs.raytrace(self._sim.tel, self._sim.atm,
                                                     self._sim.dms)
            else:
                self._sim.wfs.d_wfs[w].d_gs.raytrace(self._sim.tel, rst=reset)
                self._sim.wfs.d_wfs[w].d_gs.raytrace(self._sim.dms)

            self._sim.wfs.d_wfs[w].comp_image(noise=noise)
        self._sim.rtc.do_centroids(0)
        c = np.array(self._sim.rtc.d_control[0].d_centroids)
        return c

    def computeModalResiduals(self):
        self._sim.rtc.d_control[1].comp_dphi(self._sim.tar.d_targets[0], False)
        self._sim.rtc.do_control(1)
        v = np.array(
                self._sim.rtc.d_control[1].d_com
        )  # We compute here the residual phase on the DM modes. Gives the Equivalent volts to apply/
        if (self.P is None):
            self.modalBasis, self.P = self.getModes2VBasis("Btt")
        ai = self.P.dot(v) * 1000.  # np rms units
        return ai

    def writeConfigOnFile(self,
                          filepath=os.environ["SHESHA_ROOT"] + "/widgets/canapass.conf"):
        aodict = OrderedDict()

        aodict.update({"Fe": 1 / self._sim.config.p_loop.ittime})
        aodict.update({"teldiam": self._sim.config.p_tel.diam})
        aodict.update({"telobs": self._sim.config.p_tel.cobs})

        # WFS
        aodict.update({"nbWfs": len(self._sim.config.p_wfss)})
        aodict.update({"nbTargets": len(self._sim.config.p_targets)})
        aodict.update({"nbCam": aodict["nbWfs"]})
        aodict.update({"nbOffaxis": 0})
        aodict.update({"nbNgsWFS": 1})
        aodict.update({"nbLgsWFS": 0})
        aodict.update({"nbFigSensor": 0})
        aodict.update({"nbSkyWfs": aodict["nbWfs"]})
        aodict.update({"nbOffNgs": 0})

        # DMS
        aodict.update({"nbDms": len(self._sim.config.p_dms)})
        aodict.update({"Nactu": self._sim.rtc.d_control[0].nactu})

        # List of things
        aodict.update({"list_NgsOffAxis": []})
        aodict.update({"list_Fig": []})
        aodict.update({"list_Cam": [0]})
        aodict.update({"list_SkyWfs": [0]})
        aodict.update({"list_ITS": []})
        aodict.update({"list_Woofer": []})
        aodict.update({"list_Tweeter": []})
        aodict.update({"list_Steering": []})

        # fct of Nb of wfss
        NslopesList = []
        NsubapList = []
        listWfsType = []
        pyrModulationList = []
        pyr_npts = []
        pyr_pupsep = []
        pixsize = []
        xPosList = []
        yPosList = []
        fstopsize = []
        fstoptype = []
        npixPerSub = []
        nxsubList = []
        nysubList = []
        lambdaList = []
        dms_seen = []
        colTmpList = []
        noise = []
        new_hduwfsl = pfits.HDUList()
        new_hduwfsSubapXY = pfits.HDUList()
        for i in range(aodict["nbWfs"]):
            new_hduwfsl.append(pfits.ImageHDU(
                    self._sim.config.p_wfss[i]._isvalid))  # Valid subap array
            new_hduwfsl[i].header["DATATYPE"] = "valid_wfs%d" % i

            xytab = np.zeros((2, self._sim.config.p_wfss[i]._validsubsx.shape[0]))
            xytab[0, :] = self._sim.config.p_wfss[i]._validsubsx
            xytab[1, :] = self._sim.config.p_wfss[i]._validsubsy

            new_hduwfsSubapXY.append(
                    pfits.ImageHDU(xytab))  # Valid subap array inXx Y on the detector
            new_hduwfsSubapXY[i].header["DATATYPE"] = "validXY_wfs%d" % i

            pixsize.append(self._sim.config.p_wfss[i].pixsize)
            NslopesList.append(self._sim.config.p_wfss[i]._nvalid * 2)  # slopes per wfs
            NsubapList.append(self._sim.config.p_wfss[i]._nvalid)  # subap per wfs
            listWfsType.append(self._sim.config.p_wfss[i].type)
            xPosList.append(self._sim.config.p_wfss[i].xpos)
            yPosList.append(self._sim.config.p_wfss[i].ypos)
            fstopsize.append(self._sim.config.p_wfss[i].fssize)
            fstoptype.append(self._sim.config.p_wfss[i].fstop)
            nxsubList.append(self._sim.config.p_wfss[i].nxsub)
            nysubList.append(self._sim.config.p_wfss[i].nxsub)
            lambdaList.append(self._sim.config.p_wfss[i].Lambda)
            dms_seen.append(list(self._sim.config.p_wfss[i].dms_seen))
            noise.append(self._sim.config.p_wfss[i].noise)

            if (self._sim.config.p_wfss[i].type == b"pyrhr"):
                pyrModulationList.append(self._sim.config.p_wfss[i].pyr_ampl)
                pyr_npts.append(self._sim.config.p_wfss[i].pyr_npts)
                pyr_pupsep.append(self._sim.config.p_wfss[i].pyr_pup_sep)
                npixPerSub.append(1)
            else:
                pyrModulationList.append(0)
                pyr_npts.append(0)
                pyr_pupsep.append(0)
                npixPerSub.append(self._sim.config.p_wfss[i].npix)
        confname = filepath.split("/")[-1].split('.conf')[0]
        new_hduwfsl.writeto(
                filepath.split(".conf")[0] + '_wfsConfig.fits', overwrite=True)
        new_hduwfsSubapXY.writeto(
                filepath.split(".conf")[0] + '_wfsValidXYConfig.fits', overwrite=True)
        aodict.update({"listWFS_NslopesList": NslopesList})
        aodict.update({"listWFS_NsubapList": NsubapList})
        aodict.update({"listWFS_WfsType": listWfsType})
        aodict.update({"listWFS_pixarc": pixsize})
        aodict.update({"listWFS_pyrModRadius": pyrModulationList})
        aodict.update({"listWFS_pyrModNPts": pyr_npts})
        aodict.update({"listWFS_pyrPupSep": pyr_pupsep})
        aodict.update({"listWFS_fstopsize": fstopsize})
        aodict.update({"listWFS_fstoptype": fstoptype})
        aodict.update({"listWFS_dms_seen": dms_seen})
        aodict.update({"listWFS_NsubX": nxsubList})
        aodict.update({"listWFS_NsubY": nysubList})
        aodict.update({"listWFS_Nsub": nysubList})
        aodict.update({"listWFS_NpixPerSub": npixPerSub})
        aodict.update({"listWFS_Lambda": lambdaList})
        aodict.update({"listWFS_noise": noise})

        listDmsType = []
        NactuX = []
        unitPerVolt = []
        push4imat = []
        coupling = []
        push4iMatArcSec = []
        new_hdudmsl = pfits.HDUList()

        for j in range(aodict["nbDms"]):
            listDmsType.append(self._sim.config.p_dms[j].type)
            NactuX.append(self._sim.config.p_dms[j].nact)
            unitPerVolt.append(self._sim.config.p_dms[j].unitpervolt)
            push4imat.append(self._sim.config.p_dms[j].push4imat)
            coupling.append(self._sim.config.p_dms[j].coupling)
            tmp = []
            if (self._sim.config.p_dms[j].type != 'tt'):
                tmpdata = np.zeros((2, len(self._sim.config.p_dm0._i1)))
                tmpdata[0, :] = self._sim.config.p_dm0._j1
                tmpdata[1, :] = self._sim.config.p_dm0._i1
                new_hdudmsl.append(pfits.ImageHDU(tmpdata))  # Valid subap array
                new_hdudmsl[j].header["DATATYPE"] = "valid_dm%d" % j
            #for k in range(aodict["nbWfs"]):
            #    tmp.append(self._sim.computeDMrange(j, k))

            push4iMatArcSec.append(tmp)
        new_hdudmsl.writeto(
                filepath.split(".conf")[0] + '_dmsConfig.fits', overwrite=True)
        aodict.update({"listDMS_push4iMatArcSec": push4iMatArcSec})
        aodict.update({"listDMS_push4iMat": push4imat})
        aodict.update({"listDMS_unitPerVolt": unitPerVolt})
        aodict.update({"listDMS_Nxactu": NactuX})
        aodict.update({"listDMS_Nyactu": NactuX})
        aodict.update({"listDMS_type": listDmsType})
        aodict.update({"listDMS_coupling": coupling})

        listTargetsLambda = []
        listTargetsXpos = []
        listTargetsYpos = []
        listTargetsDmsSeen = []
        listTargetsMag = []
        for k in range(aodict["nbTargets"]):
            listTargetsLambda.append(self._sim.config.p_targets[k].Lambda)
            listTargetsXpos.append(self._sim.config.p_targets[k].xpos)
            listTargetsYpos.append(self._sim.config.p_targets[k].ypos)
            listTargetsMag.append(self._sim.config.p_targets[k].mag)
            listTargetsDmsSeen.append(self._sim.config.p_targets[k].dms_seen)

        aodict.update({"listTARGETS_Lambda": listTargetsLambda})
        aodict.update({"listTARGETS_Xpos": listTargetsXpos})
        aodict.update({"listTARGETS_Ypos": listTargetsYpos})
        aodict.update({"listTARGETS_Mag": listTargetsMag})
        aodict.update({"listTARGETS_DmsSeen": listTargetsDmsSeen})

        listDmsType = []
        Nslopes = sum(NslopesList)
        Nsubap = sum(NsubapList)
        aodict.update({"Nslopes": Nslopes})
        aodict.update({"Nsubap": Nsubap})
        f = open(filepath, 'w+')
        for dictval in aodict:
            f.write(dictval + ":" + str(aodict[dictval]) + "\n")
        f.close()
        print("OK: Config File wrote in:" + filepath)
        #return aodict

    def setIntegratorLaw(self):
        self._sim.rtc.d_control[0].set_commandlaw("integrator")

    def setDecayFactor(self, decay):
        self._sim.rtc.d_control[0].set_decayFactor(decay.astype(np.float32).copy())

    def setEMatrix(self, eMat):
        self._sim.rtc.d_control[0].set_matE(eMat.astype(np.float32).copy())

    def doRefslopes(self):
        print("Doing refslopes...")
        self._sim.rtc.do_centroids_ref(0)
        print("refslopes done")

    def resetRefslopes(self):
        self._sim.rtc.d_control[0].d_centroids_ref.reset()

    def setPyrModulation(self, pyrmod):
        CompassSupervisor.setPyrModulation(self, pyrmod)
        self._sim.rtc.do_centroids(0)  # To be ready for the next getSlopes

    def setPyrMethod(self, pyrmethod):
        CompassSupervisor.setPyrMethod(self, pyrmethod)
        self._sim.rtc.do_centroids(0)  # To be ready for the next getSlopes

    def setNoise(self, noise, numwfs=0):
        CompassSupervisor.setNoise(self, noise, numwfs)

    def setNcpaWfs(self, ncpa, wfsnum):
        self._sim.wfs.set_ncpa_phase(wfsnum, ncpa.astype(np.float32).copy())

    def setNcpaTar(self, ncpa, tarnum):
        self._sim.tar.d_targets[tarnum].set_ncpa(tarnum, ncpa.astype(np.float32).copy())

    def set_phaseWFS(self, numwfs, phase):
        pph = phase.astype(np.float32)
        self._sim.wfs.d_wfs[numwfs].set_phase(pph)
        self.computeSlopes()

    def setMpupil(self, mpupil, numwfs=0):
        oldmpup = self.getMpupil()
        dimx = oldmpup.shape[0]
        dimy = oldmpup.shape[1]
        if ((mpupil.shape[0] != dimx) or (mpupil.shape[1] != dimy)):
            print("Error mpupil shape on wfs %d must be: (%d,%d)" % (numwfs, dimx, dimy))
        else:
            self._sim.wfs.d_wfs[numwfs].set_pupil(mpupil.copy())

    def getIpupil(self):
        return self._sim.config.p_geom._ipupil

    def getSpupil(self):
        return self._sim.config.p_geom._spupil

    def getMpupil(self):
        return self._sim.config.p_geom._mpupil

    def getAmplipup(self, tarnum):
        return self._sim.config.tar.get_amplipup(tarnum)

    def getPhase(self, tarnum):
        return np.array(self._sim.tar.d_targets[tarnum].d_phase)

    def getWFSPhase(self, wfsnum):
        return np.array(self._sim.wfs.d_wfs[wfsnum].d_gs.d_phase)

    def getTargetPhase(self, tarnum):
        pup = self.getSpupil()
        ph = self.getPhase(tarnum) * pup
        return ph

    def getNcpaWfs(self, wfsnum):
        return np.array(self._sim.wfs.d_wfs[wfsnum].d_gs.d_ncpa_phase)

    def getNcpaTar(self, tarnum):
        return np.array(self._sim.tar.d_targets[tarum].d_ncpa_phase)

    #def getVolts(self):
    #    return self._sim.rtc.get_voltage(0)
    #
    #def getSlopes(self):
    #     return self._sim.rtc.get_centroids(0)

    def writeDataInFits(self, data, fullpath):
        pfits.writeto(fullpath, data, overwrite=True)

    def getFrameCounter(self):
        return self._sim.iter

    def recordCB(self, CBcount, subSample=1, tarnum=0, seeAtmos=True,
                 tarPhaseFilePath=""):
        slopesdata = None
        voltsdata = None
        tarPhaseData = None
        aiData = None
        k = 0
        # Resets the target so that the PSF LE is synchro with the data
        for i in range(len(self._sim.config.p_targets)):
            self._sim.tar.d_targets[i].reset_strehlmeter()

        # Starting CB loop...
        for j in trange(CBcount, desc="recording"):
            self._sim.next(see_atmos=seeAtmos)
            for t in range(len(self._sim.config.p_targets)):
                self._sim.tar.d_targets[t].comp_image()

            if (j % subSample == 0):
                aiVector = self.computeModalResiduals()
                if (aiData is None):
                    aiData = np.zeros((len(aiVector), int(CBcount / subSample)))
                aiData[:, k] = aiVector

                slopesVector = np.array(self._sim.rtc.d_control[0].d_centroids)
                if (slopesdata is None):
                    slopesdata = np.zeros((len(slopesVector), int(CBcount / subSample)))
                slopesdata[:, k] = slopesVector

                voltsVector = np.array(self._sim.rtc.d_control[0].d_com)
                if (voltsdata is None):
                    voltsdata = np.zeros((len(voltsVector), int(CBcount / subSample)))
                voltsdata[:, k] = voltsVector

                if (tarPhaseFilePath != ""):
                    tarPhaseArray = self.getTargetPhase(tarnum)
                    if (tarPhaseData is None):
                        tarPhaseData = np.zeros((*tarPhaseArray.shape,
                                                 int(CBcount / subSample)))
                    tarPhaseData[:, :, k] = tarPhaseArray
                k += 1
        if (tarPhaseFilePath != ""):
            print("Saving tarPhase cube at: ", tarPhaseFilePath)
            pfits.writeto(tarPhaseFilePath, tarPhaseData, overwrite=True)
        psfLE = self._sim.tar.get_image(tarnum, b"le")
        return slopesdata, voltsdata, aiData, psfLE

        #wao.sim.config.p_geom._ipupil
        """
        plt.matshow(wao.config.p_geom._ipupil, fignum=1)

        # DM positions in iPupil:
        dmposx = wao.config.p_dm0._xpos -0.5
        dmposy = wao.config.p_dm0._ypos -0.5
        plt.scatter(dmposy, dmposx, color="blue")

        #WFS position in ipupil
        ipup = wao.config.p_geom._ipupil
        spup = wao.config.p_geom._spupil
        s2ipup = (ipup.shape[0] - spup.shape[0]) / 2.
        posx = wao.config.p_wfss[0]._istart + s2ipup
        posx = np.tile(posx,(posx.size,1))
        posy = posx.T.copy()
        posx = posx *  wao.config.p_wfss[0]._isvalid
        posx = posx[np.where(posx > 0)] - ipup.shape[0] / 2 -wao.config.p_wfs0.npix
        posy = posy * wao.config.p_wfss[0]._isvalid
        posy = posy[np.where(posy > 0)] - ipup.shape[0] / 2 -wao.config.p_wfs0.npix

        #center of ssp position in ipupil
        demissp = (posx[1]-posx[0])/2
        sspx = posx+ipup.shape[0]/2+demissp
        sspy = posy+ipup.shape[0]/2+demissp
        plt.scatter(sspy, sspx, color="red")



        imat = wao.rtc.get_imat(0)

        influDM = wao.dms.get_influ(b"pzt", 0)
        influTT = wao.dms.get_influ(b"tt", 0)

        """


if __name__ == '__main__':
    from docopt import docopt
    arguments = docopt(__doc__)
    supervisor = CanapassSupervisor(arguments["<parameters_filename>"], True)
    supervisor.initConfig()

    try:
        from subprocess import Popen, PIPE
        from hraa.server.pyroServer import PyroServer

        p = Popen("whoami", shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if (err != b''):
            print(err)
            raise ValueError("ERROR CANNOT RECOGNIZE USER")
        else:
            user = out.split(b"\n")[0].decode("utf-8")
            print("User is " + user)
        server = PyroServer()
        server.add_device(supervisor, "waoconfig_" + user)
        server.start()
    except:
        raise EnvironmentError(
                "Missing dependencies (code HRAA or Pyro4 or Dill Serializer)")
