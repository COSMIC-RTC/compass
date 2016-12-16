import shesha as ao
import numpy as np

simul_name = ""
#loop
p_loop = ao.Param_loop()

p_loop.set_niter(20000)
p_loop.set_ittime(1/500.) #=1/500


#geom
p_geom=ao.Param_geom()

p_geom.set_zenithangle(0.)


#tel
p_tel=ao.Param_tel()

p_tel.set_diam(8.0)
p_tel.set_cobs(0.)


#atmos
p_atmos=ao.Param_atmos()

p_atmos.set_r0(0.16)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([5000.])
p_atmos.set_windspeed([10.0])
p_atmos.set_winddir([90])
p_atmos.set_L0([100])


#target
p_target=ao.Param_target()

p_target.set_nTargets(1)
p_target.set_xpos([0.])
p_target.set_ypos([0.])
p_target.set_Lambda([1.65])
p_target.set_mag([10])


#wfs
#p_wfs0= ao.Param_wfs()
#p_wfs1= ao.Param_wfs()
p_wfss=[ao.Param_wfs(error_budget=True)]

for i in range(len(p_wfss)):
    p_wfss[i].set_type("sh")
    p_wfss[i].set_nxsub(40)
    p_wfss[i].set_npix(4)
    p_wfss[i].set_pixsize(0.2)
    p_wfss[i].set_fracsub(0.8)
    p_wfss[i].set_xpos(0.)
    p_wfss[i].set_ypos(0.)
    p_wfss[i].set_Lambda(0.5)
    p_wfss[i].set_gsmag(6.)
    p_wfss[i].set_optthroughput(0.5)
    p_wfss[i].set_zerop(2.5e10)
    p_wfss[i].set_noise(3)
    p_wfss[i].set_atmos_seen(1)

#lgs parameters
#p_wfss[0].set_gsalt(90*1.e3)
#p_wfss[0].set_lltx(0)
#p_wfss[0].set_llty(0)
#p_wfss[0].set_laserpower(10)
#p_wfss[0].set_lgsreturnperwatt(1.e3)
#p_wfss[0].set_proftype("Exp")
#p_wfss[0].set_beamsize(0.8)

#dm
#p_dm0=ao.Param_dm()
#p_dm1=ao.Param_dm()
p_dms=[ao.Param_dm(),ao.Param_dm()]
p_dms[0].set_type("pzt")
nact=p_wfss[0].nxsub+1
p_dms[0].set_nact(nact)
p_dms[0].set_alt(0.)
p_dms[0].set_thresh(0.3)
p_dms[0].set_coupling(0.2)
p_dms[0].set_unitpervolt(1.)
p_dms[0].set_push4imat(1.)

p_dms[1].set_type("tt")
p_dms[1].set_alt(0.)
p_dms[1].set_unitpervolt(1.)
p_dms[1].set_push4imat(1.)

#centroiders
#p_centroider0=ao.Param_centroider()
p_centroiders=[ao.Param_centroider()]

for i in range(len(p_centroiders)):

    p_centroiders[i].set_nwfs(i)
    p_centroiders[i].set_type("cog")
    #p_centroiders[i].set_nmax(8)
    p_centroiders[i].set_thresh(0)

#p_centroider0.set_type("corr")
#p_centroider0.set_type_fct("model")

#controllers
p_controller1=ao.Param_controller()
p_controllers=[p_controller1]

p_controller1.set_type("ls")
p_controller1.set_nwfs([0])
p_controller1.set_ndm([0,1])
p_controller1.set_maxcond(20)
p_controller1.set_delay(0)
p_controller1.set_gain(0.3)

#rtc
p_rtc=ao.Param_rtc()

p_rtc.set_nwfs(0)
p_rtc.set_centroiders(p_centroiders)
p_rtc.set_controllers(p_controllers)
