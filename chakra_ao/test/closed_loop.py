
import cProfile
import pstats as ps

import sys
import numpy as np
import chakra as ch
import chakra_ao as ao
import time
import matplotlib.pyplot as pl

print "TEST CHAKRA_AO\n closed loop: call loop(int niter)"


if(len(sys.argv)!=2):
    error= 'command line should be:"python -i test.py parameters_filename"\n with "parameters_filename" the path to the parameters file'
    raise StandardError(error)

#get parameters from file
param_file=sys.argv[1]
execfile(param_file)
start=param_file.rindex("/")
end=param_file.rindex(".")
simul_name=param_file[start+1:end]
print "param_file is",param_file
print "simul name is",simul_name


#initialisation:
#   context
c=ch.chakra_context()
c.set_activeDevice(0)

#    wfs
print "->wfs"
wfs=ao.wfs_init(p_wfss,p_atmos,p_tel,p_geom,p_target,p_loop, 1,0,p_dms)

#   atmos
print "->atmos"
atm=ao.atmos_init(c,p_atmos,p_tel,p_geom,p_loop,p_wfss,p_target,rank=-1)

#   dm 
print "->dm"
dms=ao.dm_init(p_dms,p_wfs0,p_geom,p_tel)

#   target
print "->target"
tar=ao.target_init(c,p_target,p_atmos,p_geom,p_tel,p_wfss,wfs,p_dms)

print "->rtc"
#   rtc
rtc=ao.rtc_init(wfs,p_wfss,dms,p_dms,p_geom,p_rtc,p_atmos,atm,p_tel,p_loop,tar,p_target,simul_name=simul_name)

print "===================="
print "init done"
print "===================="
print "objects initialzed on GPU:"
print "--------------------------------------------------------"
print atm
print wfs
print dms
print tar
print rtc

mimg = 0.# initializing average image
print "----------------------------------------------------";
print "iter# | S.E. SR | L.E. SR | Est. Rem. | framerate";
print "----------------------------------------------------";

def loop( n):
    #fig,((turbu,image),(shak,defMir))=pl.subplots(2,2, figsize=(15,15))
    #pl.ion()
    #pl.show()
    t0=time.time()
    for i in range(n):
        atm.move_atmos()
        
        if(p_controller0.type_control == "geo"):
            for t in range(p_target.ntargets):
                tar.atmos_trace(t,atm)
                rtc.docontrol_geo(0, dms, tar, 0)
                rtc.applycontrol(0,dms)
                tar.dmtrace(0,dms)
        else:
            for t in range(p_target.ntargets):
                tar.atmos_trace(t,atm)
                tar.dmtrace(t,dms)
            for w in range(len(p_wfss)):
                wfs.sensors_trace(w,"all",atm,dms)
                wfs.sensors_compimg(w)

            rtc.docentroids(0)
            rtc.docontrol(0)
        
            rtc.applycontrol(0,dms)
        
        if((i+1)%100==0):
            """
            turbu.clear()
            image.clear()
            shak.clear()
            defMir.clear()

            screen=atm.get_screen(p_atmos.alt)
            f1=turbu.matshow(screen,cmap='Blues_r')

            im=tar.get_image(0,"se")
            im=np.roll(im,im.shape[0]/2,axis=0)
            im=np.roll(im,im.shape[1]/2,axis=1)
            f2=image.matshow(im,cmap='Blues_r')

            sh=wfs.get_binimg(0)
            f3=shak.matshow(sh,cmap='Blues_r')

            dm=dms.get_dm("pzt",0.)
            f4=defMir.matshow(dm)

            pl.draw()
            """
            strehltmp = tar.get_strehl(0)
            print i+1,"\t",strehltmp[0],"\t",strehltmp[1]
    t1=time.time()
    print " loop execution time:",t1-t0,"  (",n,"iterations), ",(t1-t0)/n,"(mean)  ", n/(t1-t0),"Hz"


#loop(p_loop.niter)
