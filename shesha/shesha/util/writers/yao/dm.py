YAO_DMTYPE={"pzt":"\"stackarray\"",
            "tt" :"\"tiptilt\""}

def init_dm(filename):
    f=open(filename,"a+")
    f.write("\n\n//------------------------------")
    f.write("\n//DM  parameters")
    f.write("\n//------------------------------")
    f.close()
    return 0
    
def write_dm(filename,dm,index,subSystem=1,offset=0):
    """Write (append) dm parameter to file for YAO use for a single dm

    filename : str      : name of the file to append the parameter to
    dm       : Param_dm : compass dm  parameters
    index    : int      : YAO index for dm
    """
    obj="dm("+str(index)+")"
    f=open(filename,"a+")
    f.write("\ngrow,dm,dms;")
    f.write("\n"+obj+".type          = "+YAO_DMTYPE[dm.type]+";")
    f.write("\n"+obj+".subsystem     = "+str(subSystem)+";")
    f.write("\n"+obj+".iffile        = \"\"; // not set by compass")
    f.write("\n"+obj+".alt           = "+str(dm.alt)+";")
    f.write("\n"+obj+".unitpervolt   = "+str(dm.unitpervolt)+";")
    f.write("\n"+obj+".push4imat     = "+str(dm.push4imat)+";")

    if(dm.type != "tt"):
        f.write("\n"+obj+".nxact         = "+str(dm.nact)+";")
        f.write("\n"+obj+".pitch         = "+str(dm._pitch)+";")
        f.write("\n"+obj+".thresholdresp = "+str(dm.thresh)+";")
        f.write("\n"+obj+".pitchMargin   = "+str(2.2)+"; // not set by compass")
        f.write("\n"+obj+".elt           = "+str(1)+"; // not set by compass")
        f.write("\n"+obj+".coupling    = "+str(dm.coupling)+";")
    f.close()

def write_dms(filename,dms,subSystem=1,offset=0):
    """Write (append) dm parameter to file for YAO

    filename : str       : name of the file to append the parameter to
    dms       : list[Param_dm] : compass dm  parameters list
    """
    f=open(filename,"a+")

    i=1
    for d in dms:
        f.write("\n\n//DM "+str(i+offset))
        f.flush()
        write_dm(filename,d,i+offset,subSystem=subSystem)
        i+=1

    f.close()
    return len(dms)

def finish_dm(filename,ndm):
    f=open(filename,"a+")
    f.write("\n\nndm = "+str(ndm)+";")
    f.close()
