YAO_DMTYPE={"pzt":"\"stackarray\""}

def write_dm(filename,dm,index,subSystem=1):
    """Write (append) dm parameter to file for YAO use for a single dm

    filename : str      : name of the file to append the parameter to
    dm       : Param_dm : compass dm  parameters
    index    : int      : YAO index for dm
    """
    obj="dm("+str(index)+")"
    f=open(filename,"a+")
    f.write("\n"+obj+".type          = "+YAO_DMTYPE[dm.type]+";")
    f.write("\n"+obj+".subsystem     = "+str(subSystem)+";")
    f.write("\n"+obj+".iffile        = \"\"; // not set by compass")
    f.write("\n"+obj+".nxact         = "+str(dm.nact)+";")
    f.write("\n"+obj+".pitch         = "+str(int(dm._pitch))+";")#+sim.pupildiam/wfs(n).shnxsub)
    f.write("\n"+obj+".alt           = "+str(dm.alt)+";")
    f.write("\n"+obj+".thresholdresp = "+str(dm.thresh)+";")
    f.write("\n"+obj+".unitpervolt   = "+str(dm.unitpervolt)+";")
    f.write("\n"+obj+".push4imat     = "+str(dm.push4imat)+";")
    f.write("\n"+obj+".pitchMargin   = "+str(2.2)+"; // not set by compass")
    f.write("\n"+obj+".elt           = "+str(1)+"; // not set by compass")
    f.write("\n"+obj+".coupling    = "+str(dm.coupling)+";")
    f.close()

def write_dms(filename,dms,subSystem=1):
    """Write (append) dm parameter to file for YAO

    filename : str       : name of the file to append the parameter to
    dms       : list[Param_dm] : compass dm  parameters list
    """
    ndm=len(dms)
    f=open(filename,"a+")
    f.write("\n\n//------------------------------")
    f.write("\n//DM  parameters")
    f.write("\n//------------------------------")
    f.write("\nndm = "+str(len(dms))+";")
    f.write("\ndm = array(dms,ndm);")

    i=1
    for d in dms:
        f.write("\n\n//DM "+str(i))
        f.flush()
        write_dm(filename,d,i)
        i+=1

    f.close()
