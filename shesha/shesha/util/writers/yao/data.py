from astropy.io import fits
import numpy as np

def get_yao_subapPos_single(sup,wfsId):
    """Return the coordinates of the valid subapertures of a given WFS

    this coordinates are given in meters and centered

    sup   : compass supervisor  :
    wfsId : int                 : index of the WFS
    """

    config=sup.config
    wfs=config.p_wfss[wfsId]
    geom=config.p_geom
    total=geom.pupdiam/wfs.nxsub*(wfs.nxsub-1)
    validX=wfs._validpuppixx-2
    validY=wfs._validpuppixy-2
    toMeter=(config.p_tel.diam/wfs.nxsub/wfs._pdiam)
    validX=(validX-total/2)*toMeter
    validY=(validY-total/2)*toMeter
    return validX,validY


def get_yao_subapPos(sup,nwfs=-1):

    """return the number of valid subapertures for all WFS as well as their coordinates
    
    the coordinates are given in meters and centered

    sup   : compass supervisor  :
    """

    config=sup.config
    if(nwfs<0):
        nwfs=len(config.p_wfss)
    nvalid=np.array([w._nvalid for w in config.p_wfss[:nwfs]])
    validSubap=np.zeros((2,np.sum(nvalid)))
    ind=0
    for w in range(nwfs):
        validSubap[0,ind:ind+nvalid[w]],validSubap[1,ind:ind+nvalid[w]]=get_yao_subapPos_single(sup,w)
        ind+=nvalid[w]
    return nvalid.astype(np.int32),validSubap.astype(np.float64)


def get_yao_actuPos_single(sup,dmId):
    """return the coordinates of a given DM actuators for YAO

    sup : compass supervisor
    dmId: int : index of the DM
    """

    dm=sup.config.p_dms[dmId]
    return dm._xpos+1, dm._ypos+1


def get_yao_actuPos(sup):
    """return the coordinates of all  DM actuators for YAO

    sup : compass supervisor
    """
    config=sup.config
    nactu=np.array([dm._ntotact for dm in config.p_dms])
    nactuPos=np.zeros((2,np.sum(nactu)))
    ind=0
    for dm in range(len(config.p_dms)):
        if(sup.config.p_dms[dm].type !="tt"):
            nactuPos[0,ind:ind+nactu[dm]],nactuPos[1,ind:ind+nactu[dm]]=get_yao_actuPos_single(sup,dm)
        ind+=nactu[dm]
    return nactu.astype(np.int32),nactuPos.astype(np.float64)


def write_data(fileName,sup,nwfs=-1,controllerId=0):
    """ Write data for yao compatibility

    sup : compass supervisor 

    write into a single fits:
        * number of valide subapertures
        * number of actuators
        * subapertures position (2-dim array x,y) in meters centered
        * actuator position (2-dim array x,y) in pixels starting from 0
        * interaction matrix (2*nSubap , nactu)
        * command matrix (nacy , 2*nSubap)
    """

    print("writing data to"+fileName)
    hdu=fits.PrimaryHDU(np.zeros(1,dtype=np.int32))
    config=sup.config
    nactu=config.p_controllers[controllerId].nactu

    #get nb of subap and their position
    nvalid,subapPos=get_yao_subapPos(sup,nwfs)
    nTotalValid=np.sum(nvalid)
    hdu.header["NTOTSUB"]=nTotalValid
    hdu_nsubap=fits.ImageHDU(nvalid,name="NSUBAP")
    hdu_subapPos=fits.ImageHDU(subapPos,name="SUBAPPOS")

    #get nb of actu and their position
    nactu,actuPos=get_yao_actuPos(sup)
    nTotalActu=np.sum(nactu)
    hdu.header["NTOTACTU"]=nTotalActu
    hdu_nactu=fits.ImageHDU(nactu,name="NACTU")
    hdu_actuPos=fits.ImageHDU(actuPos,name="ACTUPOS")

    #IMAT
    hdu_imat=fits.ImageHDU(sup.getImat(),name="IMAT")

    #CMAT
    hdu_cmat=fits.ImageHDU(sup.getCmat(),name="CMAT")

    print("\t* number of subaperture per WFS")
    print("\t* subapertures position")
    print("\t* number of actuator per DM")
    print("\t* actuators position")
    print("\t* Imat")
    print("\t* Cmat")

    hdul=fits.HDUList([hdu,hdu_nsubap,hdu_subapPos,hdu_nactu,hdu_actuPos,hdu_imat,hdu_cmat])
    hdul.writeto(fileName,overwrite=1)


