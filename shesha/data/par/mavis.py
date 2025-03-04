import shesha.config as conf
import numpy as np

simul_name = "mavis"

# loop
p_loop = conf.ParamLoop()

p_loop.set_niter(5000)
p_loop.set_ittime(0.0007)  # =1/1500

# geom
p_geom = conf.ParamGeom()

p_geom.set_zenithangle(0.0)

# tel
p_tel = conf.ParamTel()

p_tel.set_diam(8.0)
p_tel.set_cobs(0.12)

# atmos
p_atmos = conf.ParamAtmos()

p_atmos.set_r0(0.12)
p_atmos.set_nscreens(4)
p_atmos.set_frac([0.5, 0.2, 0.2, 0.1])
p_atmos.set_alt([0.0, 1000.0, 4500.0, 9000.0])
p_atmos.set_windspeed([10.0, 10.0, 10.0, 10.0])
p_atmos.set_winddir([0.0, 10.0, 20.0, 25.0])
p_atmos.set_L0([25.0, 25.0, 25.0, 25.0])

# target
ntargets = 49
p_targets = [conf.ParamTarget() for _ in range(49)]

xpos, ypos = np.meshgrid(np.linspace(-15, 15, 7), np.linspace(-15, 15, 7))
xpos = xpos.flatten()
ypos = ypos.flatten()
k = 0
for p_target in p_targets:
    p_target.set_xpos(xpos[k])
    p_target.set_ypos(ypos[k])
    p_target.set_Lambda(0.65)
    p_target.set_mag(10.0)
    k += 1
# wfs
n_lgs = 6
n_ngs = 3

p_wfs_lgs = [conf.ParamWfs() for _ in range(n_lgs)]
p_wfs_ngs = [conf.ParamWfs() for _ in range(n_ngs)]
p_wfss = p_wfs_lgs + p_wfs_ngs

for p_wfs in p_wfs_lgs:
    p_wfs.set_type("sh")
    p_wfs.set_nxsub(40)
    p_wfs.set_npix(8)
    p_wfs.set_pixsize(0.5)
    p_wfs.set_fracsub(0.8)

    p_wfs.set_Lambda(0.589)
    p_wfs.set_gsmag(8.0)
    p_wfs.set_optthroughput(0.5)
    p_wfs.set_zerop(1.0e11)
    p_wfs.set_noise(0.3)
    p_wfs.set_atmos_seen(1)

    p_wfs.set_gsalt(90 * 1.0e3)
    p_wfs.set_lltx(0.0)
    p_wfs.set_llty(0.0)
    p_wfs.set_laserpower(10)
    p_wfs.set_lgsreturnperwatt(1.0e3)
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
    p_wfs.set_fracsub(0.0)
    p_wfs.set_xpos(xpos_ngs[k])
    p_wfs.set_ypos(ypos_ngs[k])

    p_wfs.set_Lambda(1.5)
    p_wfs.set_gsmag(8.0)
    p_wfs.set_optthroughput(0.5)
    p_wfs.set_zerop(1.0e11)
    p_wfs.set_noise(0.3)
    p_wfs.set_atmos_seen(1)
    p_wfs.set_is_low_order(True)
    k += 1

# dm
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dm2 = conf.ParamDm()
p_dms = [p_dm0, p_dm1, p_dm2]

p_dm0.set_type("pzt")
p_dm0.set_nact(41)
p_dm0.set_alt(0.0)
p_dm0.set_thresh(0.3)
p_dm0.set_coupling(0.3)
p_dm0.set_unitpervolt(0.01)
p_dm0.set_push4imat(100.0)

p_dm1.set_type("pzt")
p_dm1.set_nact(41)
p_dm1.set_alt(4000.0)
p_dm1.set_thresh(0.0)
p_dm1.set_coupling(0.3)
p_dm1.set_unitpervolt(0.01)
p_dm1.set_push4imat(100.0)

p_dm2.set_type("pzt")
p_dm2.set_nact(61)
p_dm2.set_alt(12000.0)
p_dm2.set_thresh(0.0)
p_dm2.set_coupling(0.3)
p_dm2.set_unitpervolt(0.01)
p_dm2.set_push4imat(100.0)

# centroiders
p_centroiders = [conf.ParamCentroider() for _ in range(len(p_wfss))]

k = 0
for p_centroider in p_centroiders:
    p_centroider.set_nwfs(k)
    p_centroider.set_type("cog")
    k += 1

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("mv")
p_controller0.set_nwfs(np.arange(len(p_wfss)))
p_controller0.set_ndm([0, 1, 2])
p_controller0.set_maxcond(150.0)
p_controller0.set_delay(1.0)
p_controller0.set_gain(0.3)
