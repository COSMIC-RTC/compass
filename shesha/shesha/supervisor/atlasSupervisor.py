"""Widget to simulate a closed loop

Usage:
  AtlasSupervisor.py [<parameters_filename>] [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -d, --GPUdevices GPUdevices      Specify the GPU devices

"""

import os
import numpy as np
from collections import OrderedDict

import astropy.io.fits as pfits

import shesha.constants as scons

from .benchSupervisor import BenchSupervisor

# from naga.obj import obj_Double2D
# from naga.magma import syevd_Double, svd_host_Double
# from naga.context import context as naga_context

# from naga.host_obj import host_obj_Double1D, host_obj_Double2D


class AtlasSupervisor(BenchSupervisor):

    def __init__(self, configFile: str = None, BRAHMA: bool = True) -> None:
        BenchSupervisor.__init__(self, configFile=configFile, BRAHMA=BRAHMA)

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

    def initConfig(self) -> None:
        BenchSupervisor.initConfig(self)

        from hraa.devices.dm.kakou import Kakou
        self.dm = Kakou(reset=True)

    def getConfig(self, path=None):
        ''' Returns the configuration in use, in a supervisor specific format '''
        if path:
            self.writeConfigOnFile(path)
            return
        return BenchSupervisor.getConfig(self)

    def loadConfig(self, configFile: str = None, sim=None) -> None:
        ''' Load the configuration for the compass supervisor'''
        BenchSupervisor.loadConfig(self, configFile=configFile, sim=sim)
        print("switching to a generic controller")
        self.config.p_controllers[0].type = scons.ControllerType.GENERIC

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

    def writeConfigOnFile(self,
                          filepath=os.environ["SHESHA_ROOT"] + "/widgets/Atlas.conf"):
        aodict = OrderedDict()

        #aodict.update({"Fe": 1 / self.config.p_loop.ittime})
        #aodict.update({"teldiam": self.config.p_tel.diam})
        #aodict.update({"telobs": self.config.p_tel.cobs})

        # WFS
        aodict.update({"nbWfs": len(self.config.p_wfss)})
        aodict.update({"nbTargets": 1})
        aodict.update({"nbCam": aodict["nbWfs"]})
        aodict.update({"nbOffaxis": 0})
        aodict.update({"nbNgsWFS": 1})
        aodict.update({"nbLgsWFS": 0})
        aodict.update({"nbFigSensor": 0})
        aodict.update({"nbSkyWfs": aodict["nbWfs"]})
        aodict.update({"nbOffNgs": 0})

        # DMS
        aodict.update({"nbDms": len(self.config.p_dms)})
        aodict.update({"Nactu": self.rtc.d_control[0].nactu})

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
                    self.config.p_wfss[i]._isvalid))  # Valid subap array
            new_hduwfsl[i].header["DATATYPE"] = "valid_wfs%d" % i

            xytab = np.zeros((2, self.config.p_wfss[i]._validsubsx.shape[0]))
            xytab[0, :] = self.config.p_wfss[i]._validsubsx
            xytab[1, :] = self.config.p_wfss[i]._validsubsy

            new_hduwfsSubapXY.append(
                    pfits.ImageHDU(xytab))  # Valid subap array inXx Y on the detector
            new_hduwfsSubapXY[i].header["DATATYPE"] = "validXY_wfs%d" % i

            pixsize.append(self.config.p_wfss[i].pixsize)
            NslopesList.append(self.config.p_wfss[i]._nvalid * 2)  # slopes per wfs
            NsubapList.append(self.config.p_wfss[i]._nvalid)  # subap per wfs
            listWfsType.append(self.config.p_wfss[i].type)
            xPosList.append(self.config.p_wfss[i].xpos)
            yPosList.append(self.config.p_wfss[i].ypos)
            fstopsize.append(self.config.p_wfss[i].fssize)
            fstoptype.append(self.config.p_wfss[i].fstop)
            nxsubList.append(self.config.p_wfss[i].nxsub)
            nysubList.append(self.config.p_wfss[i].nxsub)
            lambdaList.append(self.config.p_wfss[i].Lambda)
            #dms_seen.append(list(self.config.p_wfss[i].dms_seen))
            #noise.append(self.config.p_wfss[i].noise)

            if (self.config.p_wfss[i].type == "pyrhr"):
                pyrModulationList.append(self.config.p_wfss[i].pyr_ampl)
                pyr_npts.append(self.config.p_wfss[i].pyr_npts)
                pyr_pupsep.append(self.config.p_wfss[i].pyr_pup_sep)
                npixPerSub.append(1)
            else:
                pyrModulationList.append(0)
                pyr_npts.append(0)
                pyr_pupsep.append(0)
                npixPerSub.append(self.config.p_wfss[i].npix)
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
        #aodict.update({"listWFS_dms_seen": dms_seen})
        aodict.update({"listWFS_NsubX": nxsubList})
        aodict.update({"listWFS_NsubY": nysubList})
        aodict.update({"listWFS_Nsub": nysubList})
        aodict.update({"listWFS_NpixPerSub": npixPerSub})
        aodict.update({"listWFS_Lambda": lambdaList})
        #aodict.update({"listWFS_noise": noise})

        listDmsType = []
        NactuX = []
        unitPerVolt = []
        push4imat = []
        coupling = []
        push4iMatArcSec = []
        new_hdudmsl = pfits.HDUList()

        for j in range(aodict["nbDms"]):
            listDmsType.append(self.config.p_dms[j].type)
            NactuX.append(self.config.p_dms[j].nact)
            unitPerVolt.append(self.config.p_dms[j].unitpervolt)
            push4imat.append(self.config.p_dms[j].push4imat)
            coupling.append(self.config.p_dms[j].coupling)
            tmp = []
            """
            if (self.config.p_dms[j].type != 'tt'):
                tmpdata = np.zeros((2, len(self.config.p_dm0._i1)))
                tmpdata[0, :] = self.config.p_dm0._j1
                tmpdata[1, :] = self.config.p_dm0._i1
                new_hdudmsl.append(pfits.ImageHDU(tmpdata))  # Valid subap array
                new_hdudmsl[j].header["DATATYPE"] = "valid_dm%d" % j
            """
            #for k in range(aodict["nbWfs"]):
            #    tmp.append(self.computeDMrange(j, k))

            #push4iMatArcSec.append(tmp)
        new_hdudmsl.writeto(
                filepath.split(".conf")[0] + '_dmsConfig.fits', overwrite=True)
        #aodict.update({"listDMS_push4iMatArcSec": push4iMatArcSec})
        #aodict.update({"listDMS_push4iMat": push4imat})
        aodict.update({"listDMS_unitPerVolt": unitPerVolt})
        aodict.update({"listDMS_Nxactu": NactuX})
        aodict.update({"listDMS_Nyactu": NactuX})
        aodict.update({"listDMS_type": listDmsType})
        #aodict.update({"listDMS_coupling": coupling})

        listTargetsLambda = []
        listTargetsXpos = []
        listTargetsYpos = []
        listTargetsDmsSeen = []
        listTargetsMag = []
        """
        for k in range(aodict["nbTargets"]):
            listTargetsLambda.append(self.config.p_targets[k].Lambda)
            listTargetsXpos.append(self.config.p_targets[k].xpos)
            listTargetsYpos.append(self.config.p_targets[k].ypos)
            listTargetsMag.append(self.config.p_targets[k].mag)
            listTargetsDmsSeen.append(list(self.config.p_targets[k].dms_seen))

        aodict.update({"listTARGETS_Lambda": listTargetsLambda})
        aodict.update({"listTARGETS_Xpos": listTargetsXpos})
        aodict.update({"listTARGETS_Ypos": listTargetsYpos})
        aodict.update({"listTARGETS_Mag": listTargetsMag})
        aodict.update({"listTARGETS_DmsSeen": listTargetsDmsSeen})
        """
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

    def writeDataInFits(self, data, fullpath):
        pfits.writeto(fullpath, data, overwrite=True)


from threading import Thread, Event


class MyThread(Thread):

    def __init__(self, event, func):
        Thread.__init__(self)
        self.stopped = event
        self.func = func

    def run(self):
        while not self.stopped.wait(0.01):
            self.func()
        print("THREAD STOPPED")


if __name__ == '__main__':
    from docopt import docopt
    from hraa.server.pyroServer import PyroServer

    arguments = docopt(__doc__)
    supervisor = AtlasSupervisor(arguments["<parameters_filename>"], True)
    USEPYRO = True
    if arguments["--GPUdevices"]:
        supervisor.config.p_loop.set_devices([
                int(device) for device in arguments["--GPUdevices"].split(",")
        ])

    supervisor.initConfig()
    # supervisor.loop(supervisor.config.p_loop.niter)

    listDevices = [supervisor]
    listNames = ['ATLAS_SUPERVISOR']

    srv = None
    if USEPYRO:
        ATLAS_PYRO_PORT = 6667
        srv = PyroServer(
                bindTo=('wsmiao2.obspm.fr',
                        ATLAS_PYRO_PORT), nsAddress=('wsmiao2.obspm.fr', 6666),
                listDevices=listDevices, listNames=listNames)
        srv.start()

    stop = Event()
    thread = MyThread(stop, supervisor.singleNext)
    thread.start()
