from collections import OrderedDict
import os
import astropy.io.fits as pfits
import numpy as np
"""
This file parse the compass parameters ((wao object)) and writes the configuration in a text file to be read by any software.
Arrays are written in .fits files

"""


def returnConfigfromWao(wao,
                        filepath=os.environ["SHESHA_ROOT"] + "/widgets/canapass.conf"):
    aodict = OrderedDict()

    aodict.update({"Fe": 1 / wao.sim.config.p_loop.ittime})
    aodict.update({"teldiam": wao.sim.config.p_tel.diam})
    aodict.update({"telobs": wao.sim.config.p_tel.cobs})

    # WFS
    aodict.update({"nbWfs": len(wao.sim.config.p_wfss)})
    aodict.update({"nbCam": aodict["nbWfs"]})
    aodict.update({"nbOffaxis": 0})
    aodict.update({"nbNgsWFS": 1})
    aodict.update({"nbLgsWFS": 0})
    aodict.update({"nbFigSensor": 0})
    aodict.update({"nbSkyWfs": aodict["nbWfs"]})
    aodict.update({"nbOffNgs": 0})

    # DMS
    aodict.update({"Ndm": len(wao.sim.config.p_dms)})
    aodict.update({"Nactu": sum(wao.sim.config.p_controllers[0].nactu)})

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
    new_hduwfsl = pfits.HDUList()
    new_hduwfsSubapXY = pfits.HDUList()
    for i in range(aodict["nbWfs"]):
        new_hduwfsl.append(
                pfits.ImageHDU(wao.sim.config.p_wfss[i]._isvalid))  # Valid subap array
        new_hduwfsl[i].header["DATATYPE"] = "valid_wfs%d" % i

        xytab = np.zeros((wao.sim.config.p_wfss[i]._validsubsx.shape[0],
                          wao.sim.config.p_wfss[i]._validsubsy.shape[0]))
        xytab[0] = wao.sim.config.p_wfss[i]._validsubsx
        xytab[1] = wao.sim.config.p_wfss[i]._validsubsy

        new_hduwfsSubapXY.append(
                pfits.ImageHDU(xytab))  # Valid subap array inXx Y on the detector
        new_hduwfsSubapXY[i].header["DATATYPE"] = "validXY_wfs%d" % i

        pixsize.append(wao.sim.config.p_wfss[i].pixsize)
        NslopesList.append(wao.sim.config.p_wfss[i]._nvalid * 2)  # slopes per wfs
        NsubapList.append(wao.sim.config.p_wfss[i]._nvalid)  # subap per wfs
        listWfsType.append(wao.sim.config.p_wfss[i].type)
        xPosList.append(wao.sim.config.p_wfss[i].xpos)
        yPosList.append(wao.sim.config.p_wfss[i].ypos)
        fstopsize.append(wao.sim.config.p_wfss[i].fssize)
        fstoptype.append(wao.sim.config.p_wfss[i].fstop)
        nxsubList.append(wao.sim.config.p_wfss[i].nxsub)
        nysubList.append(wao.sim.config.p_wfss[i].nxsub)
        lambdaList.append(wao.sim.config.p_wfss[i].Lambda)
        dms_seen.append(list(wao.sim.config.p_wfss[i].dms_seen))

        if (wao.sim.config.p_wfss[i].type == "pyrhr"):
            pyrModulationList.append(wao.sim.config.p_wfss[i].pyr_ampl)
            pyr_npts.append(wao.sim.config.p_wfss[i].pyr_npts)
            pyr_pupsep.append(wao.sim.config.p_wfss[i].pyr_pup_sep)
            npixPerSub.append(1)
        else:
            pyrModulationList.append(0)
            pyr_npts.append(0)
            pyr_pupsep.append(0)
            npixPerSub.append(wao.sim.config.p_wfss[i].npix)
    confname = filepath.split("/")[-1].split('.conf')[0]
    new_hduwfsl.writeto(filepath.split(".conf")[0] + '_wfsConfig.fits', clobber=True)
    new_hduwfsSubapXY.writeto(
            filepath.split(".conf")[0] + '_wfsValidXYConfig.fits', clobber=True)
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

    listDmsType = []
    NactuX = []
    unitPerVolt = []
    push4imat = []
    coupling = []
    push4iMatArcSec = []
    new_hdudmsl = pfits.HDUList()

    for j in range(aodict["Ndm"]):
        listDmsType.append(wao.sim.config.p_dms[j].type)
        NactuX.append(wao.sim.config.p_dms[j].nact)
        unitPerVolt.append(wao.sim.config.p_dms[j].unitpervolt)
        push4imat.append(wao.sim.config.p_dms[j].push4imat)
        coupling.append(wao.sim.config.p_dms[j].coupling)
        tmp = []
        if (wao.sim.config.p_dms[j].type != 'tt'):
            tmpdata = np.zeros((2, len(wao.sim.config.p_dm0._i1)))
            tmpdata[0, :] = wao.sim.config.p_dm0._j1
            tmpdata[1, :] = wao.sim.config.p_dm0._i1
            new_hdudmsl.append(pfits.ImageHDU(tmpdata))  # Valid subap array
            new_hdudmsl[j].header["DATATYPE"] = "valid_dm%d" % j
        #for k in range(aodict["nbWfs"]):
        #    tmp.append(wao.sim.computeDMrange(j, k))

        push4iMatArcSec.append(tmp)
    new_hdudmsl.writeto(filepath.split(".conf")[0] + '_dmsConfig.fits', clobber=True)
    aodict.update({"listDMS_push4iMatArcSec": push4iMatArcSec})
    aodict.update({"listDMS_push4iMat": push4imat})
    aodict.update({"listDMS_unitPerVolt": unitPerVolt})
    aodict.update({"listDMS_Nxactu": NactuX})
    aodict.update({"listDMS_Nyactu": NactuX})
    aodict.update({"listDMS_type": listDmsType})
    aodict.update({"listDMS_coupling": coupling})

    listDmsType = []
    Nslopes = sum(NslopesList)
    Nsubap = sum(NsubapList)
    aodict.update({"Nslopes": Nslopes})
    aodict.update({"Nsubap": Nsubap})
    f = open(filepath, 'w+')
    for dictval in aodict:
        f.write(dictval + ":" + str(aodict[dictval]) + "\n")
    f.close()
    #return aodict