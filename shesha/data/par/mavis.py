import shesha.config as conf
import numpy as np
simul_name = "mavis"

# loop
p_loop = conf.Param_loop()

p_loop.set_niter(5000)
p_loop.set_ittime(0.0007)  # =1/500

# geom
p_geom = conf.Param_geom()

p_geom.set_zenithangle(0.)

# tel
p_tel = conf.Param_tel()

p_tel.set_diam(8.0)
p_tel.set_cobs(0.12)

# atmos
p_atmos = conf.Param_atmos()

p_atmos.set_r0(0.12)
p_atmos.set_nscreens(4)
p_atmos.set_frac([0.5, 0.2, 0.2, 0.1])
p_atmos.set_alt([0., 1000.0, 4500., 9000.])
p_atmos.set_windspeed([10.0, 10.0, 10.0, 10.0])
p_atmos.set_winddir([0., 10., 20., 25.])
p_atmos.set_L0([25., 25., 25., 25.])

# target
p_target = conf.Param_target()

p_target.set_ntargets(49)
xpos, ypos = np.meshgrid(np.linspace(-15, 15, 7), np.linspace(-15, 15, 7))
p_target.set_xpos(xpos.flatten())
p_target.set_ypos(ypos.flatten())
p_target.set_Lambda([0.65])
p_target.set_mag([10.])

# wfs
p_wfs_lgs = [conf.Param_wfs()] * 6
p_wfs_ngs = [conf.Param_wfs()] * 3
p_wfss = p_wfs_lgs + p_wfs_ngs

for p_wfs in p_wfs_lgs:
    p_wfs.set_type("sh")
    p_wfs.set_nxsub(40)
    p_wfs.set_npix(8)
    p_wfs.set_pixsize(0.5)
    p_wfs.set_fracsub(0.8)

    p_wfs.set_Lambda(0.589)
    p_wfs.set_gsmag(8.)
    p_wfs.set_optthroughput(0.5)
    p_wfs.set_zerop(1.e11)
    p_wfs.set_noise(0.3)
    p_wfs.set_atmos_seen(1)

    p_wfs.set_gsalt(90 * 1.e3)
    p_wfs.set_lltx(0.)
    p_wfs.set_llty(0.)
    p_wfs.set_laserpower(10)
    p_wfs.set_lgsreturnperwatt(1.e3)
    p_wfs.set_proftype("Gauss1")
    p_wfs.set_beamsize(0.8)

x = np.linspace(0, 2 * np.pi, 7)[:-1]
radius = 15
lgs_xpos = radius * np.cos(x)
lgs_ypos = radius * np.sin(x)
for k in range(len(p_wfs_lgs)):
    p_wfs_lgs[k].set_xpos(lgs_xpos[k])
    p_wfs_lgs[k].set_ypos(lgs_ypos[k])

xpos_ngs = [-20, 20, 20]
ypos_ngs = [0, 20, -20]
k = 0
for p_wfs in p_wfs_ngs:
    p_wfs.set_type("sh")
    p_wfs.set_nxsub(2)
    p_wfs.set_npix(20)
    p_wfs.set_pixsize(0.030)
    p_wfs.set_fracsub(0.8)
    p_wfs.set_xpos(xpos_ngs[k])
    p_wfs.set_ypos(ypos_ngs[k])

    p_wfs.set_Lambda(1.5)
    p_wfs.set_gsmag(8.)
    p_wfs.set_optthroughput(0.5)
    p_wfs.set_zerop(1.e11)
    p_wfs.set_noise(0.3)
    p_wfs.set_atmos_seen(1)
    k += 1

# dm
p_dm0 = conf.Param_dm()
p_dm1 = conf.Param_dm()
p_dm2 = conf.Param_dm()
p_dms = [p_dm0, p_dm1, p_dm2]

p_dm0.set_type("pzt")
p_dm0.set_nact(41)
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.3)
p_dm0.set_coupling(0.3)
p_dm0.set_unitpervolt(0.01)
p_dm0.set_push4imat(100.)

p_dm1.set_type("pzt")
p_dm1.set_nact(41)
p_dm1.set_alt(4000.)
p_dm1.set_thresh(0.)
p_dm1.set_coupling(0.3)
p_dm1.set_unitpervolt(0.01)
p_dm1.set_push4imat(100.)

p_dm2.set_type("pzt")
p_dm2.set_nact(61)
p_dm2.set_alt(12000.)
p_dm2.set_thresh(0.)
p_dm2.set_coupling(0.3)
p_dm2.set_unitpervolt(0.01)
p_dm2.set_push4imat(100.)

# centroiders
p_centroiders = [conf.Param_centroider()] * 9

k = 0
for p_centroider in p_centroiders:
    p_centroider.set_nwfs(k)
    p_centroider.set_type("cog")
    k += 1

# controllers
p_controller0 = conf.Param_controller()
p_controllers = [p_controller0]

p_controller0.set_type("generic")
p_controller0.set_nwfs(np.arange(9))
p_controller0.set_ndm([0, 1, 2])
p_controller0.set_maxcond(150.)
p_controller0.set_delay(1.)
p_controller0.set_gain(0.3)
