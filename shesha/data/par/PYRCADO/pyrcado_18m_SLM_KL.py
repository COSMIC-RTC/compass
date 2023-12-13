import shesha.config as conf

simul_name = ""

#loop
p_loop = conf.ParamLoop()

p_loop.set_niter(1000)  #500Hz: 1mn = 30000, 1kH = 60000
p_loop.set_ittime(1 / 500.)  #=1/500

#geom
p_geom = conf.ParamGeom()
p_geom.set_zenithangle(0.)

#tel
p_tel = conf.ParamTel()
p_tel.set_diam(18.)

#atmos
p_atmos = conf.ParamAtmos()

p_atmos.set_r0(0.129)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.])
p_atmos.set_winddir([45.])
p_atmos.set_L0([25.])  # Not simulated in Yorick?

# No target

# Fake WFS for raytracing
p_wfs0 = conf.ParamWfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")

p_wfs0.set_nxsub(38)

p_wfs0.set_npix(4)
p_wfs0.set_pixsize(0.3)
p_wfs0.set_fracsub(0.8)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(1.00)
p_wfs0.set_gsmag(15.5)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(1.e11)
p_wfs0.set_noise(-1)  # in electrons units
p_wfs0.set_atmos_seen(1)
p_wfs0.set_fstop("round")
p_wfs0.set_fssize(2.4)

#dm
p_dm0 = conf.ParamDm()
p_dms = [p_dm0]
p_dm0.set_type("kl")
p_dm0.set_nkl(1000)
p_dm0.set_kl_type("kolmo")

p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(1.)

p_controller0 = conf.ParamController()
p_controllers = [p_controller0]
p_controller0.set_ndm([0])
