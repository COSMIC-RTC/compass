import shesha.config as conf

simul_name = "micado_39m_SH"
# loop
p_loop = conf.ParamLoop()

p_loop.set_niter(30000)  # 500Hz: 1mn = 30000, 1kH = 60000
p_loop.set_ittime(1 / 500.0)  # =1/500
# p_loop.set_devices([0, 1, 2, 3])

# geom
p_geom = conf.ParamGeom()

p_geom.set_zenithangle(0.0)

# tel
p_tel = conf.ParamTel()

p_tel.set_diam(39.0)
p_tel.set_cobs(0.30)

# atmos
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

# atmos
p_atmos = conf.ParamAtmos()

p_atmos.set_r0(0.129)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.0])
p_atmos.set_winddir([45.0])
p_atmos.set_L0([1e6])  # Not simulated in Yorick?
"""
# 1 target
p_target=conf.ParamTarget()
p_target.set_nTargets(1)
p_target.set_xpos([0])
p_target.set_ypos([0.])
p_target.set_Lambda([1.65])
p_target.set_mag([4.])
"""

# 3 targets
p_targets = [conf.ParamTarget() for _ in range(3)]
Lambda = [1.2, 1.65, 2.2]
k = 0
for p_target in p_targets:
    p_target.set_xpos(0.0)
    p_target.set_ypos(0.0)
    p_target.set_Lambda(Lambda[k])
    p_target.set_mag(4.0)

# wfs

p_wfs0 = conf.ParamWfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")
p_wfs0.set_nxsub(78)
p_wfs0.set_npix(6)
p_wfs0.set_pixsize(0.3)
p_wfs0.set_fracsub(0.8)
p_wfs0.set_xpos(0.0)
p_wfs0.set_ypos(0.0)
p_wfs0.set_Lambda(0.5)
p_wfs0.set_gsmag(11)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(1.0e11)
p_wfs0.set_noise(0.1)  # in electrons units
p_wfs0.set_atmos_seen(1)
p_wfs0.set_fstop("round")
p_wfs0.set_fssize(2.4)

# dm
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type("pzt")
nact = p_wfs0.nxsub + 1
# nact = 9

p_dm0.set_nact(nact)
p_dm0.set_alt(0.0)
p_dm0.set_thresh(0.3)  # fraction units
p_dm0.set_coupling(0.2)  # !!!!!!!!!!!!!!!!!!!!!!!!! attention pas autre chose que 0.2 !!!!!!!!!
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(1.0)

p_dm1.set_type("tt")
p_dm1.set_alt(0.0)
p_dm1.set_unitpervolt(1)
p_dm1.set_push4imat(0.05)

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
# p_centroider0.set_type("cog")
# sp_centroider0.set_thresh(0.)

p_centroider0.set_type("bpcog")
p_centroider0.set_nmax(6)

# p_centroider0.set_type("corr")
# p_centroider0.set_type_fct("model")

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("ls")
p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(150.0)
p_controller0.set_delay(1)
p_controller0.set_gain(0.3)

# p_controller0.set_modopti(1)
# p_controller0.set_nrec(2048)
# p_controller0.set_nmodes(3532)
# p_controller0.set_gmin(0.001)
# p_controller0.set_gmax(0.5)
# p_controller0.set_ngain(500)
