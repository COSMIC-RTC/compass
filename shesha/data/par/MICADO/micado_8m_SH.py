import shesha as ao

simul_name = "micado_8m_SH"

#loop
p_loop = ao.Param_loop()

p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.)  #=1/500

#geom
p_geom = ao.Param_geom()
p_geom.set_zenithangle(0.)

#tel
p_tel = ao.Param_tel()
p_tel.set_diam(8.)
p_tel.set_cobs(0.1)
"""
E_ELT PUPIL
p_tel.set_type_ap("EELT-Nominal")
p_tel.set_spiders_type("six")
p_tel.set_pupangle(0)

p_tel.set_nbrmissing(7)
p_tel.set_referr(0.01)
p_tel.set_std_tt(0.050) # microns
p_tel.set_std_piston(0.050) # microns
"""

#atmos
p_atmos = ao.Param_atmos()

p_atmos.set_r0(0.129)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.])
p_atmos.set_winddir([45.])
p_atmos.set_L0([1.e5])  # Not simulated in Yorick?

#target
p_target = ao.Param_target()
p_targets = [p_target]
p_target.set_ypos(0.)
p_target.set_Lambda(1.65)
p_target.set_mag(4.)

#wfs
#p_wfs0= ao.Param_wfs(error_budget=True)
p_wfs0 = ao.Param_wfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")
p_wfs0.set_nxsub(16)
p_wfs0.set_npix(8)
p_wfs0.set_pixsize(0.3)
p_wfs0.set_fracsub(0.8)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.5)
p_wfs0.set_gsmag(15.)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(1.e11)
p_wfs0.set_noise(-1)  # in electrons units
p_wfs0.set_atmos_seen(1)
p_wfs0.set_fstop("square")
p_wfs0.set_fssize(1.6)

#dm
p_dm0 = ao.Param_dm()
p_dm1 = ao.Param_dm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type("pzt")
nact = p_wfs0.nxsub + 1
#nact = 9

p_dm0.set_nact(nact)
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.3)  # fraction units
p_dm0.set_coupling(
        0.2)  # !!!!!!!!!!!!!!!!!!!!!!!!! attention pas autre chose que 0.2 !!!!!!!!!
p_dm0.set_unitpervolt(0.01)
p_dm0.set_push4imat(100.)

p_dm1.set_type("tt")
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(0.0005)
p_dm1.set_push4imat(10.)

#centroiders
p_centroider0 = ao.Param_centroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
p_centroider0.set_type("cog")
#sp_centroider0.set_thresh(0.)

#p_centroider0.set_type("bpcog")
#p_centroider0.set_nmax(10)

#p_centroider0.set_type("corr")
#p_centroider0.set_type_fct("model")

#controllers
p_controller0 = ao.Param_controller()
p_controllers = [p_controller0]

p_controller0.set_type("ls")
p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(150.)
p_controller0.set_delay(1)
p_controller0.set_gain(0.2)

#p_controller0.set_modopti(0)
#p_controller0.set_nrec(2048)
#p_controller0.set_nmodes(5064)
#p_controller0.set_gmin(0.001)
#p_controller0.set_gmax(0.5)
#p_controller0.set_ngain(500)

#rtc
p_rtc = ao.Param_rtc()

p_rtc.set_nwfs(1)
p_rtc.set_centroiders(p_centroiders)
p_rtc.set_controllers(p_controllers)
