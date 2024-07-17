import shesha.config as conf

simul_name = ""
# loop
p_loop = conf.ParamLoop()

p_loop.set_niter(20000)
p_loop.set_ittime(1 / 100)  # =1/500

# geom
p_geom = conf.ParamGeom()

p_geom.set_zenithangle(0.0)

# tel
p_tel = conf.ParamTel()

p_tel.set_diam(8.0)
p_tel.set_cobs(0.2)

# atmos
p_atmos = conf.ParamAtmos()

p_atmos.set_r0(0.12)
p_atmos.set_nscreens(3)
p_atmos.set_frac([0.5, 0.3, 0.2])
p_atmos.set_alt([0, 5000.0, 9000.0])
p_atmos.set_windspeed([5.0, 10.0, 15.0])
p_atmos.set_winddir([0, 120, 240])
p_atmos.set_L0([100, 100, 100])
p_atmos.set_seeds([1234, 1235, 1236])

# target
p_target = conf.ParamTarget()
p_targets = [p_target]
p_target.set_xpos(0.0)
p_target.set_ypos(0.0)
p_target.set_Lambda(1.65)
p_target.set_mag(10.0)

# wfs
# p_wfs0= conf.ParamWfs()
# p_wfs1= conf.ParamWfs()
p_wfss = [conf.ParamWfs(roket=True)]

for i in range(len(p_wfss)):
    p_wfss[i].set_type("sh")
    p_wfss[i].set_nxsub(40)
    p_wfss[i].set_npix(6)
    p_wfss[i].set_pixsize(0.5)
    p_wfss[i].set_fracsub(0.8)
    p_wfss[i].set_xpos(5.0)
    p_wfss[i].set_ypos(0.0)
    p_wfss[i].set_Lambda(0.5)
    p_wfss[i].set_gsmag(9.0)
    p_wfss[i].set_optthroughput(0.5)
    p_wfss[i].set_zerop(3.0e10)
    p_wfss[i].set_noise(-1)
    p_wfss[i].set_atmos_seen(1)

# lgs parameters
# p_wfss[0].set_gsalt(90*1.e3)
# p_wfss[0].set_lltx(0)
# p_wfss[0].set_llty(0)
# p_wfss[0].set_laserpower(10)
# p_wfss[0].set_lgsreturnperwatt(1.e3)
# p_wfss[0].set_proftype("Exp")
# p_wfss[0].set_beamsize(0.8)

# dm
# p_dm0=conf.ParamDm()
# p_dm1=conf.ParamDm()
p_dms = [conf.ParamDm(), conf.ParamDm()]
p_dms[0].set_type("pzt")
nact = p_wfss[0].nxsub + 1
p_dms[0].set_nact(nact)
p_dms[0].set_alt(0.0)
p_dms[0].set_thresh(0.3)
p_dms[0].set_coupling(0.2)
p_dms[0].set_unitpervolt(1.0)
p_dms[0].set_push4imat(1.0)

p_dms[1].set_type("tt")
p_dms[1].set_alt(0.0)
p_dms[1].set_unitpervolt(1.0)
p_dms[1].set_push4imat(1.0)

# centroiders
# p_centroider0=conf.ParamCentroider()
p_centroiders = [conf.ParamCentroider()]

for i in range(len(p_centroiders)):
    p_centroiders[i].set_nwfs(i)
    p_centroiders[i].set_type("cog")
    # p_centroiders[i].set_nmax(8)
    p_centroiders[i].set_thresh(0)

# p_centroider0.set_type("corr")
# p_centroider0.set_type_fct("model")

# controllers
p_controller1 = conf.ParamController()
p_controllers = [p_controller1]

p_controller1.set_type("ls")
p_controller1.set_nwfs([0])
p_controller1.set_ndm([0, 1])
p_controller1.set_maxcond(20)
p_controller1.set_delay(0)
p_controller1.set_gain(0.3)
