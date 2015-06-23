
import cProfile
import pstats as ps

import sys
import numpy as np
import chakra as ch
import chakra_ao as ao
import time

print "TEST CHAKRA_AO (no mpi)"



#loop
p_loop = ao.Param_loop()

p_loop.set_niter(2000)
p_loop.set_ittime(0.002) #=1/500


#geom
p_geom=ao.Param_geom()

p_geom.set_zenithangle(0.)


#tel
p_tel=ao.Param_tel()

p_tel.set_diam(4.)
p_tel.set_cobs(0.30)


#atmos
p_atmos=ao.Param_atmos()

p_atmos.set_r0(0.16)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([20.0])
p_atmos.set_winddir([45])
p_atmos.set_L0([1.e3])


#target
p_target=ao.Param_target()

p_target.set_nTargets(1)
p_target.set_xpos([0])
p_target.set_ypos([0.])
p_target.set_Lambda([1.65])
p_target.set_mag([10])


#wfs
p_wfs0= ao.Param_wfs()
p_wfs1= ao.Param_wfs()
p_wfss=[p_wfs0]#,p_wfs1]

p_wfs0.set_type("sh")
p_wfs0.set_nxsub(8)#8
p_wfs0.set_npix(8)#8#20
p_wfs0.set_pixsize(0.3)#0.3#1.0
p_wfs0.set_fracsub(0.8)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.5)
p_wfs0.set_gsmag(5.)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(1.e11)
p_wfs0.set_noise(-1)
p_wfs0.set_atmos_seen(1)
p_wfs0.set_zerop(1e11)

p_wfs1.set_type("pyr")
p_wfs1.set_nxsub(8)
p_wfs1.set_npix(8)#8#20
p_wfs1.set_pixsize(0.3)#0.3#1.0
p_wfs1.set_fracsub(0.8)
p_wfs1.set_xpos(0.)
p_wfs1.set_ypos(0.)
p_wfs1.set_Lambda(0.5)
p_wfs1.set_gsmag(5.)
p_wfs1.set_optthroughput(0.5)
p_wfs1.set_zerop(1.e11)
p_wfs1.set_noise(-1)
p_wfs1.set_atmos_seen(1)
p_wfs1.set_zerop(1e11)
p_wfs1.set_fstop("round")

#dm
p_dm0=ao.Param_dm()
p_dm1=ao.Param_dm()
p_dms=[p_dm0,p_dm1]

p_dm0.set_type("pzt")
nact=p_wfs0.nxsub+1
p_dm0.set_nact(nact)
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.3)
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(0.01)
p_dm0.set_push4imat(100.)

p_dm1.set_type("tt")
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(0.0005)
p_dm1.set_push4imat(10.)



#centroiders
p_centroider0=ao.Param_centroider()
p_centroiders=[p_centroider0]

p_centroider0.set_nwfs(0)
p_centroider0.set_type("cog")
#p_centroider0.set_type_fct("gauss")

#controllers
p_controller0=ao.Param_controller()
p_controllers=[p_controller0]

p_controller0.set_type("ls")
p_controller0.set_nwfs([1])
p_controller0.set_ndm([1,2])
p_controller0.set_maxcond(150)
p_controller0.set_delay(1)
p_controller0.set_gain(0.5)


#rtc
p_rtc=ao.Param_rtc()

p_rtc.set_nwfs(1)
p_rtc.set_centroiders(p_centroiders)
p_rtc.set_controllers(p_controllers)

c=ch.chakra_context()
c.set_activeDevice(0)#0%c.get_ndevice())

#initialisation:
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
