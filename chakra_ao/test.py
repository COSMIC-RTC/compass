
import cProfile
import pstats as ps

import sys
import numpy as np
import chakra as ch
import chakra_ao as ao
import time

print "TEST CHAKRA_AO (no mpi)"


if(len(sys.argv)!=2):
    error= 'command line should be:"python -i test.py parameters_filename"\n with "parameters_filename" the path to the parameters file'
    raise StandardError(error)

#get parameters from file
param_file=sys.argv[1]
execfile(param_file)

#initialisation:
#   context
c=ch.chakra_context()
c.set_activeDevice(0)

#    wfs
print "->wfs"
wfs=ao.wfs_init(p_wfss,p_atmos,p_tel,p_geom,p_target,p_loop, 1,0,p_dms)

#   atmos
print "->atmos"
atm=p_atmos.atmos_init(c,p_tel,p_geom,p_loop)

#   dm 
print "->dm"
dms=ao.dm_init(p_dms,p_wfs0,p_geom,p_tel)

#   target
print "->target"
tar=p_target.target_init(c,p_atmos,p_geom,p_tel,p_wfss,wfs,p_dms)

print "->rtc"
#   rtc
rtc=ao.rtc_init(wfs,p_wfss,dms,p_dms,p_geom,p_rtc,p_atmos,atm,p_loop,p_target)

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


start = time.time()
#generate turbulence, raytrace through the atmos
# and display:
#   the turbulence
#   the target's image 
#ao.see_atmos_target(5,atm,tar,alt=0,n_tar=0,f=0.15)
ao.see_atmos_target_disp(1,atm,tar,wfs,alt=0,n_tar=0,f=0.15)


#for profiling purpose
#filename="Profile.prof"
#cProfile.runctx("ao.see_atmos_target(10,atm,tar,wfs,alt=0,n_tar=0,f=0.15)", globals(), locals(),filename)

end = time.time()

#for profiling purpose
#s = ps.Stats(filename)
#s.strip_dirs().sort_stats("time").print_stats()
#end = time.time()

print end - start

#shak=wfs.get_binimg(0)
#screen=atm.get_screen(0)
#img=tar.get_image(0,"se")

#np.save("wfs_binimg",shak)
#np.save("target_image",img)
#np.save("atmos_screen",screen)
#use to compare result with non mpi exec
# load numpy array with numpy.load("wfs_binimg.npy")
