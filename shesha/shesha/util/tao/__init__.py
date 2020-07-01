from shesha.ao import imats 
from shesha.ao import cmats 
from astropy.io import fits 
import numpy as np 
import os
from shutil import copyfile
from shesha.util import write_sysParam 
from shesha.util import fits_io


from . import ltao
from . import mcao
from importlib import reload
reload(ltao)
reload(mcao)

TILESIZE="1000"

STARPU_FLAGS=""

#variable necessary to run TAO
VARS={"SCHED":"dmdas",
      "STARPU_FLAGS":"",
      "GPUIDS":0,
      "TILESIZE":1000,
      "INPUTPATH":0,
      "TAOPATH":0
      }


def check():
    """Checks that variable are initialized
    """
    stop=0
    try :
        if (not isinstance(VARS["SCHED"], str)):
            print("you must select a scheduler (dmda,dmdas,dmdar...)\n\tex: VARS[\"SCHED\"]=\"dmdas\"")
            stop=1
    except:
        print("you must select a scheduler (dmda,dmdas,dmdar...)\n\tex: VARS[\"SCHED\"]=\"dmdas\"")
        stop=1
    try :
        if( not isinstance(VARS["GPUIDS"], str)):
            print("you must define the GPUs to use as a string \n\tex:VARS[\"GPUIDS\"]=\"1,2\"")
            stop=1
    except:
        print("you must define the GPUs to use as a string \n\tex:VARS[\"GPUIDS\"]=\"1,2\"")
        stop=1
    try :
        if( not isinstance(VARS["INPUTPATH"], str)):
            print("you must define the location of the system parameters \n\tex: VARS[\"INPUTPATH\"]=\"~/workspace/compass/params\"")
            stop=1
    except:
        print("you must define the location of the system parameters \n\tex: VARS[\"INPUTPATH\"]=\"~/workspace/compass/params\"")
        stop=1
    try :
        if( not isinstance(VARS["TAOPATH"], str)):
            print("you must define the location of the tao executables \n\tex: VARS[\"TAOPATH\"]=\"~/workspace/tao/install/bin\"")
            stop=1
    except:
        print("you must define the location of the tao executables \n\tex: VARS[\"TAOPATH\"]=\"~/workspace/tao/install/bin\"")
        stop=1
    try :
        STARPU_FLAGS
    except:
        STARPU_FLAGS=""

    return stop


def init(sup,mod,WFS="all",DM_TT=False,nfilt=None):
    """ Set up the compass loop

    set the interaction matrix, loop gain and write parameter files for TAO

    sup : CompassSupervisor :
    mod : module            : AO mode requested (among: ltao , mcao)
    """

    #setting open loop
    sup.rtc._rtc.d_control[0].set_polc(True)

    if nfilt is None:
        mod.init(VARS,sup,DM_TT=DM_TT,WFS=WFS)
    else:
        mod.init(VARS,sup,DM_TT=DM_TT,WFS=WFS,nfilt=nfilt)

def reconstructor(mod):
    """ Compute the TAO reconstructor for a given AO mode

    mod : module    : AO mode requested (among: ltao , mcao)
    """
    return mod.reconstructor(VARS)

def updateCmat(sup, cmatFile):
    """ Update the compass command matrix from an input fits file
    
    sup         : CompassSupervisor :
    cmatFile    : str               : name of the cmat fits file
    """
    M=fits_io.fitsread(cmatFile).T
    sup.setCommandMatrix(M)
    return M


def run(sup,mod,nIter=1000,initialisation=0,reset=1,WFS="all",DM_TT=False,nfilt=None):
    check()

    #setting open loop
    sup.rtc._rtc.d_control[0].set_polc(True)

    #if generic: need to update imat in controller
    if(np.abs(np.array(sup.rtc._rtc.d_control[0].d_imat)).max()==0):
        #imat not set yet for controller
        sup.rtc._rtc.d_control[0].set_imat(sup.config.p_controllers[0]._imat)
    #update gain
    sup.rtc._rtc.set_gain(0,sup.config.p_controllers[0].gain)

    if(initialisation):
        init(sup,mod,WFS=WFS,DM_TT=DM_TT,nfilt=nfilt)
    M=reconstructor(mod)
    if(reset):
        sup.reset()
    cmat_shape=sup.rtc.get_command_matrix(0).shape
    if(M.shape[0] != cmat_shape[0] or M.shape[1] != cmat_shape[1]):
        print("ToR shape is not valid:\n\twaiting for:",cmat_shape,"\n\tgot        :",M.shape)
    else:
        sup.rtc.set_command_matrix(0,M)
        if(nIter>0):
            sup.loop(nIter)
