import shesha.config as conf

simul_name = "error_80x80_8pix"

# loop
p_loop = conf.ParamLoop()

p_loop.set_niter(5000)
p_loop.set_ittime(0.002)  # =1/500

# geom
p_geom = conf.ParamGeom()

p_geom.set_zenithangle(0.0)

# tel
p_tel = conf.ParamTel()

p_tel.set_diam(40.0)
p_tel.set_cobs(0.20)

# atmos
p_atmos = conf.ParamAtmos()

p_atmos.set_r0(0.16)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([20.0])
p_atmos.set_winddir([45])
p_atmos.set_L0([1.0e5])

# target
p_target = conf.ParamTarget()

p_target.set_nTargets(1)
p_target.set_xpos([0])
p_target.set_ypos([0.0])
p_target.set_Lambda([1.65])
p_target.set_mag([10])

# wfs
p_wfs0 = conf.ParamWfs(roket=True)
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")
p_wfs0.set_nxsub(80)
p_wfs0.set_npix(8)
p_wfs0.set_pixsize(0.3)
p_wfs0.set_fracsub(0.5)
p_wfs0.set_xpos(0.0)
p_wfs0.set_ypos(0.0)
p_wfs0.set_Lambda(0.5)
p_wfs0.set_gsmag(5.0)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(1.0e11)
p_wfs0.set_noise(3)
p_wfs0.set_atmos_seen(1)

# dm
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type("pzt")
nact = p_wfs0.nxsub + 1
p_dm0.set_nact(nact)
p_dm0.set_alt(0.0)
p_dm0.set_thresh(0.3)
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(1.0)
p_dm0.set_push4imat(1.0)

p_dm1.set_type("tt")
p_dm1.set_alt(0.0)
p_dm1.set_unitpervolt(1.0)
p_dm1.set_push4imat(1.0)

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
p_centroider0.set_type("cog")
# p_centroider0.set_type("corr")
# p_centroider0.set_type_fct("model")

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("ls")
p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(1500.0)
p_controller0.set_delay(0)
p_controller0.set_gain(0.4)

p_controller0.set_modopti(0)
p_controller0.set_nrec(2048)
p_controller0.set_nmodes(5064)
p_controller0.set_gmin(0.001)
p_controller0.set_gmax(0.5)
p_controller0.set_ngain(500)

# rtc
p_rtc = conf.Param_rtc()

p_rtc.set_nwfs(1)
p_rtc.set_centroiders(p_centroiders)
p_rtc.set_controllers(p_controllers)
