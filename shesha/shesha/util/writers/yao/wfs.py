import numpy as np

YAO_WFSTYPE={"sh":"\"hartmann\"", "pyrhr":"\"pyramid\""}

def write_wfs(filename,wfs,index,subSystem=1):
    """Write (append) wfs parameter to file for YAO use for a single wfs

    filename : str       : name of the file to append the parameter to
    wfs      : Param_wfs : compass wfs parameters
    """
    obj="wfs("+str(index)+")"
    f=open(filename,"a+")
    f.write("\n"+obj+".type           = "+YAO_WFSTYPE[wfs.type]+";")
    f.write("\n"+obj+".subsystem     = "+ str(subSystem)+";")
    f.write("\n"+obj+".shmethod      = 2" +";")
    f.write("\n"+obj+".shnxsub       = "+ str(wfs.nxsub)+";")
    f.write("\n"+obj+".lambda        = "+ str(wfs.Lambda)+";")
    f.write("\n"+obj+".pixsize       = "+ str(wfs.pixsize)+";")
    f.write("\n"+obj+".npixels       = "+ str(wfs.npix)+";")
    f.write("\n"+obj+".shthreshold   = 0;   // not set by compass")
    f.write("\n"+obj+".dispzoom      = 1.0; // not set by compass")
    f.write("\n"+obj+".fracIllum     = "+str(wfs.fracsub)+";")
    f.write("\n"+obj+".gspos         = [ "+str(wfs.xpos)+" , "+str(wfs.ypos)+" ];")
    if(wfs.noise<0):
        f.write("\n"+obj+".noise         = 1;")
        f.write("\n"+obj+".ron           = 0;")
    else:
        f.write("\n"+obj+".noise         = 1;")
        f.write("\n"+obj+".ron           = "+ str(wfs.noise)+";")
    f.write("\n"+obj+".darkcurrent   = 0 ; // not set by compass ")
    if(wfs.gsalt>0):
        f.write("\n"+obj+".gsalt         = "+ str(wfs.gsalt)+";")
        f.write("\n"+obj+".gsdepth       = "+ str(1)+";") #Eventually 10000
        f.write("\n"+obj+".optthroughput = "+ str(wfs.optthroughput)+";")
        f.write("\n"+obj+".laserpower    = "+ str(wfs.laserpower)+";") #20.; // 75 ph/subap
        f.write("\n"+obj+".filtertilt    = "+ str(1)+";") 
        f.write("\n"+obj+".correctUpTT   = "+ str(1)+";") 
        f.write("\n"+obj+".uplinkgain    = "+ str(0.2)+";")
    f.close()


def write_wfss(filename,wfss,nwfs=-1,subSystem=1,offset=0):
    """Write (append) wfs parameter to file for YAO use for a wfs list

    filename : str             : name of the file to append the parameter to
    wfss     :list[ Param_wfs] : compass wfs parameters list
    """
    #counting nb of lgs and ngs
    nngs=0
    nlgs=0
    if(nwfs<0):
        nwfs=len(wfss)
    for w in wfss[:nwfs]:
        if(w.gsalt>0):
            nlgs+=1
        else:
            nngs+=1
    nwfs=nngs+nlgs
    f=open(filename,"a+")
    f.write("\n\n//------------------------------")
    f.write("\n//WFS parameters")
    f.write("\n//------------------------------")
    f.write("\nnngs = "+str(nngs)+";")
    f.write("\nnlgs = "+str(nlgs)+";")
    f.write("\nnwfs = "+str(nwfs)+";")
    f.write("\nwfs = array(wfss,nwfs);")

    i=1
    for w in wfss[:nwfs] :
        f.write("\n\n//WFS"+str(i))
        f.flush()
        write_wfs(filename,w,i+offset,subSystem=subSystem)
        i+=1

    f.close()

################################

