## @package   shesha.supervisor.canapassSupervisor
## @brief     Initialization and execution of a CANAPASS supervisor
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   4.4.1
## @date      2011/01/28
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
"""
Initialization and execution of a CANAPASS supervisor

Usage:
  canapassSupervisor.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h, --help          Show this help message and exit
  -f, --freq freq       change the frequency of the loop
  -d, --delay delay     change the delay of the loop
"""

import os, sys
import numpy as np
import time
from collections import OrderedDict

from tqdm import tqdm
import astropy.io.fits as pfits
from threading import Thread
from subprocess import Popen, PIPE

import shesha.ao as ao
import shesha.constants as scons
from shesha.constants import CentroiderType, WFSType

from typing import Any, Dict, Tuple, Callable, List
from shesha.supervisor.compassSupervisor import CompassSupervisor

# from carmaWrap.obj import obj_Double2D
# from carmaWrap.magma import syevd_Double, svd_host_Double
# from carmaWrap.context import context as carmaWrap_context

# from carmaWrap.host_obj import host_obj_Double1D, host_obj_Double2D


class CanapassSupervisor(CompassSupervisor):

    def __init__(self, config_file: str = None, cacao: bool = True) -> None:
        CompassSupervisor.__init__(self, config_file=config_file, cacao=cacao)

    def get_config(self):
        ''' Returns the configuration in use, in a supervisor specific format '''
        return CompassSupervisor.get_config(self)

    def get_configFab(self):
        aodict = OrderedDict()
        dataDict = {}
        if (root is None):
            root = self

        if (root.config.p_tel is not None):
            aodict.update({"teldiam": root.config.p_tel.diam})
            aodict.update({"telobs": root.config.p_tel.cobs})
            aodict.update({"pixsize": root.config.p_geom._pixsize})
            # TURBU
            aodict.update({"r0": root.config.p_atmos.r0})
            aodict.update({"Fe": 1 / root.config.p_loop.ittime})
            aodict.update({"nbTargets": len(root.config.p_targets)})
        else:
            aodict.update({"nbTargets": 1})

        # WFS
        aodict.update({"nbWfs": len(root.config.p_wfss)})
        aodict.update({"nbCam": aodict["nbWfs"]})
        aodict.update({"nbOffaxis": 0})
        aodict.update({"nbNgsWFS": 1})
        aodict.update({"nbLgsWFS": 0})
        aodict.update({"nbFigSensor": 0})
        aodict.update({"nbSkyWfs": aodict["nbWfs"]})
        aodict.update({"nbOffNgs": 0})

        # DMS
        aodict.update({"nbDms": len(root.config.p_dms)})
        aodict.update({"Nactu": root.rtc.d_control[0].nactu})
        # List of things
        aodict.update({"list_NgsOffAxis": []})
        aodict.update({"list_Fig": []})
        aodict.update({"list_Cam": [0]})
        aodict.update({"list_SkyWfs": [0]})
        aodict.update({"list_ITS": []})
        aodict.update({"list_Woofer": []})
        aodict.update({"list_Tweeter": []})
        aodict.update({"list_Steering": []})

        listOfNstatesPerController = []
        listOfcontrolLawTypePerController = []
        for control in self.config.p_controllers:
            listOfNstatesPerController.append(control.nstates)
            listOfcontrolLawTypePerController.append(control.type)
        aodict.update({"list_nstatesPerController": listOfNstatesPerController})
        aodict.update({"list_controllerType": listOfcontrolLawTypePerController})

        # fct of Nb of wfss
        NslopesList = []
        NsubapList = []
        listWfsType = []
        listCentroType = []

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
        #new_hduwfsl = pfits.HDUList()
        #new_hduwfsSubapXY = pfits.HDUList()
        for i in range(aodict["nbWfs"]):
            #new_hduwfsl.append(pfits.ImageHDU(root.config.p_wfss[i]._isvalid))  # Valid subap array
            #new_hduwfsl[i].header["DATATYPE"] = "valid_wfs%d" % i
            dataDict["wfsValid_" + str(i)] = root.config.p_wfss[i]._isvalid

            xytab = np.zeros((2, root.config.p_wfss[i]._validsubsx.shape[0]))
            xytab[0, :] = root.config.p_wfss[i]._validsubsx
            xytab[1, :] = root.config.p_wfss[i]._validsubsy
            dataDict["wfsValidXY_" + str(i)] = xytab

            #new_hduwfsSubapXY.append(pfits.ImageHDU(xytab))  # Valid subap array inXx Y on the detector
            #new_hduwfsSubapXY[i].header["DATATYPE"] = "validXY_wfs%d" % i
            pixsize.append(root.config.p_wfss[i].pixsize)
            """
            if (root.config.p_centroiders[i].type == "maskedpix"):
                factor = 4
            else:
                factor = 2
            NslopesList.append(
                    root.config.p_wfss[i]._nvalid * factor)  # slopes per wfs
            """
            listCentroType.append(
                    root.config.p_centroiders[i].
                    type)  # assumes that there is the same number of centroiders and wfs
            NsubapList.append(root.config.p_wfss[i]._nvalid)  # subap per wfs
            listWfsType.append(root.config.p_wfss[i].type)
            xPosList.append(root.config.p_wfss[i].xpos)
            yPosList.append(root.config.p_wfss[i].ypos)
            fstopsize.append(root.config.p_wfss[i].fssize)
            fstoptype.append(root.config.p_wfss[i].fstop)
            nxsubList.append(root.config.p_wfss[i].nxsub)
            nysubList.append(root.config.p_wfss[i].nxsub)
            lambdaList.append(root.config.p_wfss[i].Lambda)
            if (root.config.p_wfss[i].dms_seen is not None):
                dms_seen.append(list(root.config.p_wfss[i].dms_seen))
                noise.append(root.config.p_wfss[i].noise)

            if (root.config.p_centroiders[i].type == CentroiderType.MASKEDPIX):
                NslopesList.append(root.config.p_wfss[i]._nvalid * 4)  # slopes per wfs
            else:
                NslopesList.append(root.config.p_wfss[i]._nvalid * 2)  # slopes per wfs

            if (root.config.p_wfss[i].type == "pyrhr"):
                pyrModulationList.append(root.config.p_wfss[i].pyr_ampl)
                pyr_npts.append(root.config.p_wfss[i].pyr_npts)
                pyr_pupsep.append(root.config.p_wfss[i].pyr_pup_sep)
                npixPerSub.append(1)
            else:
                pyrModulationList.append(0)
                pyr_npts.append(0)
                pyr_pupsep.append(0)
                npixPerSub.append(root.config.p_wfss[i].npix)
        """
        confname = filepath.split("/")[-1].split('.conf')[0]
        print(filepath.split(".conf")[0] + '_wfsConfig.fits')
        new_hduwfsl.writeto(
                filepath.split(".conf")[0] + '_wfsConfig.fits', overwrite=True)
        new_hduwfsSubapXY.writeto(
                filepath.split(".conf")[0] + '_wfsValidXYConfig.fits', overwrite=True)
        """
        if (len(dms_seen) != 0):
            aodict.update({"listWFS_dms_seen": dms_seen})

        aodict.update({"listWFS_NslopesList": NslopesList})
        aodict.update({"listWFS_NsubapList": NsubapList})
        aodict.update({"listWFS_CentroType": listCentroType})
        aodict.update({"listWFS_WfsType": listWfsType})
        aodict.update({"listWFS_pixarc": pixsize})
        aodict.update({"listWFS_pyrModRadius": pyrModulationList})
        aodict.update({"listWFS_pyrModNPts": pyr_npts})
        aodict.update({"listWFS_pyrPupSep": pyr_pupsep})
        aodict.update({"listWFS_fstopsize": fstopsize})
        aodict.update({"listWFS_fstoptype": fstoptype})
        aodict.update({"listWFS_NsubX": nxsubList})
        aodict.update({"listWFS_NsubY": nysubList})
        aodict.update({"listWFS_Nsub": nysubList})
        aodict.update({"listWFS_NpixPerSub": npixPerSub})
        aodict.update({"listWFS_Lambda": lambdaList})
        if (len(noise) != 0):
            aodict.update({"listWFS_noise": noise})

        listDmsType = []
        NactuX = []
        Nactu = []
        unitPerVolt = []
        push4imat = []
        coupling = []
        push4iMatArcSec = []
        #new_hdudmsl = pfits.HDUList()

        for j in range(aodict["nbDms"]):
            listDmsType.append(root.config.p_dms[j].type)
            NactuX.append(
                    root.config.p_dms[j].nact)  # nb of actuators across the diameter !!
            Nactu.append(root.config.p_dms[j]._ntotact)  # nb of actuators in total
            unitPerVolt.append(root.config.p_dms[j].unitpervolt)
            push4imat.append(root.config.p_dms[j].push4imat)
            coupling.append(root.config.p_dms[j].coupling)
            tmp = []
            if (root.config.p_dms[j]._i1 is
                        not None):  # Simu Case where i1 j1 is known (simulated)
                if (root.config.p_dms[j].type != 'tt'):
                    tmpdata = np.zeros((4, len(root.config.p_dms[j]._i1)))
                    tmpdata[0, :] = root.config.p_dms[j]._j1
                    tmpdata[1, :] = root.config.p_dms[j]._i1
                    tmpdata[2, :] = root.config.p_dms[j]._xpos
                    tmpdata[3, :] = root.config.p_dms[j]._ypos
                else:
                    tmpdata = np.zeros((4, 2))

                dataDict["dmData" + str(j)] = tmpdata
                """
                new_hdudmsl.append(pfits.ImageHDU(tmpdata))  # Valid subap array
                new_hdudmsl[j].header["DATATYPE"] = "valid_dm%d" % j
                """
                #for k in range(aodict["nbWfs"]):
                #    tmp.append(root.computeDMrange(j, k))

                push4iMatArcSec.append(tmp)

        # new_hdudmsl.writeto(filepath.split(".conf")[0] + '_dmsConfig.fits', overwrite=True)
        if (len(push4iMatArcSec) != 0):
            aodict.update({"listDMS_push4iMat": push4imat})
            aodict.update({"listDMS_unitPerVolt": unitPerVolt})
        aodict.update({"listDMS_Nxactu": NactuX})
        aodict.update({"listDMS_Nyactu": NactuX})
        aodict.update({"listDMS_Nactu": Nactu})

        aodict.update({"listDMS_type": listDmsType})
        aodict.update({"listDMS_coupling": coupling})

        if (root.config.p_targets is not None):  # simu case
            listTargetsLambda = []
            listTargetsXpos = []
            listTargetsYpos = []
            listTargetsDmsSeen = []
            listTargetsMag = []
            listTARGETS_pixsize = []
            for k in range(aodict["nbTargets"]):
                listTargetsLambda.append(root.config.p_targets[k].Lambda)
                listTargetsXpos.append(root.config.p_targets[k].xpos)
                listTargetsYpos.append(root.config.p_targets[k].ypos)
                listTargetsMag.append(root.config.p_targets[k].mag)
                listTargetsDmsSeen.append(list(root.config.p_targets[k].dms_seen))
                PSFPixsize = (root.config.p_targets[k].Lambda * 1e-6) / (
                        root.config.p_geom._pixsize *
                        root.config.p_geom.get_ipupil().shape[0]) * 206265.
                listTARGETS_pixsize.append(PSFPixsize)

            aodict.update({"listTARGETS_Lambda": listTargetsLambda})
            aodict.update({"listTARGETS_Xpos": listTargetsXpos})
            aodict.update({"listTARGETS_Ypos": listTargetsYpos})
            aodict.update({"listTARGETS_Mag": listTargetsMag})
            aodict.update({"listTARGETS_DmsSeen": listTargetsDmsSeen})
            aodict.update({"listTARGETS_pixsize": listTARGETS_pixsize})

        listDmsType = []
        Nslopes = sum(NslopesList)
        Nsubap = sum(NsubapList)
        aodict.update({"Nslopes": Nslopes})
        aodict.update({"Nsubap": Nsubap})
        return aodict, dataDict

    def load_config(self, config_file: str = None, sim=None) -> None:
        ''' Load the configuration for the compass supervisor'''
        CompassSupervisor.load_config(self, config_file=config_file, sim=sim)
        print("switching to a generic controller")
        self._sim.config.p_controllers[0].type = scons.ControllerType.GENERIC

    def write_config_on_file(self, root=None):
        """
        TODO : je sais pas
        """




########################## PROTO #############################

    def initModalGain(self, gain, cmatModal, modal_basis, control=0, reset_gain=True):
        """
        Given a gain, cmat and btt2v initialise the modal gain mode
        """
        print("TODO: A RECODER !!!!")
        nmode_total = modal_basis.shape[1]
        nactu_total = modal_basis.shape[0]
        nfilt = nmode_total - cmatModal.shape[0]
        ctrl = self._sim.rtc.d_control[control]
        ctrl.set_commandlaw('modal_integrator')
        cmat = np.zeros((nactu_total, cmatModal.shape[1]))
        dec = cmat.shape[0] - cmatModal.shape[0]
        cmat[:-dec, :] += cmatModal  # Fill the full Modal with all non-filtered modes
        modes2V = np.zeros((nactu_total, nactu_total))
        dec2 = modes2V.shape[1] - modal_basis.shape[1]
        modes2V[:, :-dec2] += modal_basis
        mgain = np.ones(len(modes2V)) * gain  # Initialize the gain
        ctrl.set_matE(modes2V)
        ctrl.set_cmat(cmat)
        if reset_gain:
            ctrl.set_modal_gains(mgain)

    def leaveModalGain(self, control=0):
        ctrl = self._sim.rtc.d_control[control]
        ctrl.set_commandlaw('integrator')




if __name__ == '__main__':
    from docopt import docopt
    arguments = docopt(__doc__)
    supervisor = CanapassSupervisor(arguments["<parameters_filename>"], cacao=True)
    if (arguments["--freq"]):
        print("Warning changed frequency loop to: ", arguments["--freq"])
        supervisor.config.p_loop.set_ittime(1 / float(arguments["--freq"]))
    if (arguments["--delay"]):
        print("Warning changed delay loop to: ", arguments["--delay"])
        supervisor.config.p_controllers[0].set_delay(float(arguments["--delay"]))
    supervisor.init_config()

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
