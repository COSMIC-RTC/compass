
from shesha.util.writers.yao.general import *
from shesha.util.writers.yao.wfs     import *
from shesha.util.writers.yao.dm      import *
from shesha.util.writers.yao.targets import *
from shesha.util.writers.yao.atmos   import *
from shesha.util.writers.yao.loop    import *
from shesha.util.writers.yao.gs      import *
from shesha.util.writers.yao.data    import *

def write_parfiles(sup, *,  param_file_name="./yao.par",
                            fits_file_name="./yao.fits",
                            screen_dir="./yao_screen",
                            n_wfs=-1,
                            imat_type="controller"):
    """Write parameter files for YAO simulations

    Parameters:
        sup : (CompassSupervisor) : supervisor

        param_file_name : (str) : (optional), default "./yao.par" name of the yao parameter file

        fits_file_name : (str) : (optional), default "./yao.fits" name of fits file containing sub-apertures and actuator position

        screen_dir : (str) : (optional), default "./yao_screen" path to the yao turbulent screen files

        n_wfs : (int) : (optional), default -1 number of WFS (-1: all wfs)

        imat_type : (str) : (optional), default "controller" use of regular controller or split tomography (among "controller", "splitTomo")
    """
    conf = sup.config
    zerop = conf.p_wfss[0].zerop
    lgs_return_per_watt = max([w.lgsreturnperwatt for w in conf.p_wfss])

    print("writing parameter file to " + param_file_name)
    write_general(param_file_name, conf.p_geom, conf.p_controllers,
                  conf.p_tel, conf.simul_name)
    wfs_offset = 0
    dm_offset = 0
    ndm = init_dm(param_file_name)
    for sub_system, c in enumerate(conf.p_controllers):
        dms = [ conf.p_dms[i]  for i in c.get_ndm() ]
        ndm += write_dms (param_file_name, dms ,sub_system=sub_system + 1,
                        offset=dm_offset)
        dm_offset = dm_offset+len(dms)
    finish_dm(param_file_name, ndm)
    gs = init_wfs(param_file_name)
    for sub_system, c in enumerate(conf.p_controllers):
        wfss = [ conf.p_wfss[i] for i in c.get_nwfs()]
        n_ngs, n_lgs = write_wfss(param_file_name, wfss, sub_system=sub_system + 1,
                                n_wfs=n_wfs, offset=wfs_offset)
        gs = (gs[0] + n_ngs, gs[1] + n_lgs)
        wfs_offset = wfs_offset + len(wfss)
    finish_wfs(param_file_name, gs[0], gs[1])
    write_targets(param_file_name, conf.p_targets)
    write_gs(param_file_name, zerop, lgs_return_per_watt,
             conf.p_geom.zenithangle)
    write_atm(param_file_name, conf.p_atmos, screen_dir)
    write_loop(param_file_name, conf.p_loop, conf.p_controllers[0])
    write_data(fits_file_name, sup, n_wfs=n_wfs, compose_type=imat_type)
