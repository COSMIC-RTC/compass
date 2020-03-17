
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
    write_general(paramfile,conf.p_geom,conf.p_controllers,conf.p_tel,conf.simul_name)
    wfs_offset=0
    dm_offset=0
    ndm=init_dm(paramfile)
    for subSyst, c in enumerate(conf.p_controllers):
        dms =[ conf.p_dms[i]  for i in c.get_ndm() ]
        ndm+=write_dms (paramfile,dms ,subSystem=subSyst+1,offset=dm_offset)
        dm_offset=dm_offset+len(dms)
    finish_dm(paramfile,ndm)
    gs=init_wfs(paramfile)
    for subSyst, c in enumerate(conf.p_controllers):
        wfss=[ conf.p_wfss[i] for i in c.get_nwfs()]
        nngs,nlgs=write_wfss(paramfile,wfss,subSystem=subSyst+1,nwfs=nwfs,offset=wfs_offset)
        gs=(gs[0]+nngs,gs[1]+nlgs)
        wfs_offset=wfs_offset+len(wfss)
    finish_wfs(paramfile,gs[0],gs[1])
    write_targets(paramfile,conf.p_targets)
    write_gs(paramfile,zerop,lgsreturnperwatt,conf.p_geom.zenithangle)
    write_atm(paramfile,conf.p_atmos,screenfile)
    write_loop(paramfile,conf.p_loop,conf.p_controllers[0])
    write_data(fitsfile,sup,nwfs=nwfs)
