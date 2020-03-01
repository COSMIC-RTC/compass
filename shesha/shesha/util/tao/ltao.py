from shesha.ao import imats 
from shesha.ao import cmats 
from astropy.io import fits 
import numpy as np 
from shesha.util import write_sysParam 
import os
from shesha.util import fits_io

def init(VARS,sup,nfilt=10):
    """Initialize the LTAO mode

    compute meta matrix of interaction / command and write parameter files

    VARS
    """
 
    #compute meta imat 
    metaD=imats.get_metaD(sup,0,0) 
    #get svd of (D.T*D) 
    SVD=cmats.svd_for_cmat(metaD) 
    #plt.plot(SVD[1]) 
    nfilt=10 
    metaDx=cmats.get_cmat(metaD,nfilt=nfilt,svd=SVD) 

    #write MOAO pipeline inputs 
    dataPath=VARS["INPUTPATH"]
    write_sysParam.generate_files(sup,dataPath,singleFile=True,dm_tt=False) 
    write_sysParam.write_metaDx(metaDx,nTS=sup.config.NTS,path=dataPath) 


def reconstructor(VARS,applyLog="./log"):
    flags=VARS["STARPU_FLAGS"]
    taoPath=VARS["TAOPATH"]
    dataPath=VARS["INPUTPATH"]
    gpus=VARS["GPUIDS"]
    ts=str(VARS["TILESIZE"])
    applyCmd=flags+" "+taoPath+"/ltao_reconstructor --sys_path="+dataPath+" --atm_path="+dataPath+" --ncores=1 --gpuIds="+gpus+" --ts="+ts+" --sync=1 --warmup=0  >"+applyLog+" 2>&1"
    os.system(applyCmd)
    return fits_io.fitsread("M_ltao_0.fits").T
