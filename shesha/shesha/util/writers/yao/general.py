import numpy as np

def write_general(filename,geom,controller,tel,simulname):
    """Write (append) cgeneral simulation parameter to file for YAO use

    filename    : str              : name of the file to append the parameter to
    geom        : Param_geom       : compass AO geometry parameters  
    controller  : Param_controller : compass controller parameters  
    tel         : Param_tel        : compass telescope parameters  
    simulname   : str              : simulation name
    """
    f=open(filename,"w")
    f.write("\n\n//------------------------------")
    f.write("\n//general parameters")
    f.write("\n//------------------------------")
    f.write("\nsim.name        = \""+simulname+"\";")
    f.write("\nsim.pupildiam   = "+str(geom.pupdiam)+";")
    f.write("\nsim.debug       = 0;")
    f.write("\nsim.verbose     = 1;")

    f.write("\nmat.file            = \"\";")
    f.write("\nmat.condition = &(["+str(np.sqrt(controller.maxcond))+"]);")
    f.write("\nmat.method = \"none\";")
    #f.write("\nhfield = 15")
    f.write("\nYAO_SAVEPATH = \"\"; // where to save the output to the simulations")

    f.write("\ntel.diam = "+str(tel.diam)+";")
    f.write("\ntel.cobs = "+str(tel.cobs)+";")
