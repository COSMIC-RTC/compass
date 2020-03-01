
from .general import *
from .wfs     import *
from .dm      import *
from .targets import *
from .atmos   import *
from .loop    import *
from .gs      import *
from .data    import *

def write_parfiles(sup, paramfile="./yao.par",
                        fitsfile="./yao.fits",
                        screenfile="./yao_screen",
                        nwfs=-1):
    """Write parameter files for YAO simulations

    sup         : CompassSupervisor :
    paramfile   : str               : name of the yao parameter file
    fitsfile    : str               : name of fits file containing sub-apertures and actuator position
    screenfile  : str               : path to the yao turbulent screen files
    nwfs        : int               : (optional) number of WFS
    """
    conf=sup.config
    zerop=conf.p_wfss[0].zerop
    lgsreturnperwatt=max([w.lgsreturnperwatt for w in conf.p_wfss])

    print("writing parameter file to "+paramfile)
    write_general(paramfile,conf.p_geom,conf.p_controllers[0],conf.p_tel,conf.simul_name)
    write_wfss(paramfile,conf.p_wfss,nwfs=nwfs)
    write_dms(paramfile,conf.p_dms)
    write_targets(paramfile,conf.p_targets)
    write_atm(paramfile,conf.p_atmos,screenfile)
    write_loop(paramfile,conf.p_loop,conf.p_controllers[0])
    write_gs(paramfile,zerop,lgsreturnperwatt,conf.p_geom.zenithangle)
    write_data(fitsfile,sup,nwfs=nwfs)
