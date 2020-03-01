
def write_gs(filename,zeropoint,lgsreturnperwatt,zenithangle):
    """Write (append) guide stars parameters to file for YAO

    filename            : str       : name of the file to append the parameter to
    zeropoint           : float     : flux for magnitude 0 (ph/m²/s)
    lgsreturnperwatt    : float     : return per watt factor (ph/cm²/s/W)
    zenithangle         : float     : zenithal angle (degree)
    """
    f=open(filename,"a+")
    f.write("\n\n//------------------------------")
    f.write("\n//GS parameters")
    f.write("\n//------------------------------")

    f.write("\ngs.zeropoint         = "+ str(zeropoint)+"; //TODO get ")# Consider later (ngs intensity)
    f.write("\ngs.lgsreturnperwatt  = "+ str(lgsreturnperwatt)+"; //TODO check lgs case")
    f.write("\ngs.zenithangle       = "+ str(zenithangle)+";")


