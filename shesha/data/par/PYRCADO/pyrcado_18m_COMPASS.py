import shesha.config as conf
import shesha.constants as scons
import numpy as np

simul_name = ""

# loop
p_loop = conf.ParamLoop()

p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.0)  # =1/500
p_loop.set_devices([4, 5, 6, 7])

# geom
p_geom = conf.ParamGeom()
p_geom.set_zenithangle(0.0)
# p_geom.set_pupdiam(304)

# tel
p_tel = conf.ParamTel()
p_tel.set_diam(18.0)

# atmos
p_atmos = conf.ParamAtmos()

p_atmos.set_r0(0.129)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.0])
p_atmos.set_winddir([45.0])
p_atmos.set_L0([25.0])  # Not simulated in Yorick?

# Lambda target(s)
p_targets = [conf.ParamTarget() for _ in range(2)]
Lambda = [0.658, 1.65]
k = 0
for p_target in p_targets:
    p_target.set_xpos(0.0)
    p_target.set_ypos(0.0)
    p_target.set_Lambda(Lambda[k])
    p_target.set_mag(4.0)
    k += 1

# wfs
p_wfs0 = conf.ParamWfs()
p_wfss = [p_wfs0]

p_wfs0.set_type(scons.WFSType.PYRHR)
p_wfs0.set_nxsub(38)
p_wfs0.set_fracsub(0.0001)
p_wfs0.set_xpos(0.0)
p_wfs0.set_ypos(0.0)
p_wfs0.set_Lambda(0.658)
p_wfs0.set_gsmag(16)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(2.6e10)  # 2.6e10 ph/s/m**2 computed by Rico in R band for MOSAIC
p_wfs0.set_noise(0.3)  # in electrons units
p_wfs0.set_atmos_seen(1)
p_wfs0.set_fstop(scons.FieldStopType.SQUARE)
p_wfs0.set_fssize(1.6)

rMod = 8

p_wfs0.set_pyr_npts(int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.0) * 4))
p_wfs0.set_pyr_ampl(rMod)
p_wfs0.set_pyr_pup_sep((32))

# dm
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type(scons.DmType.PZT)
nact = p_wfs0.nxsub + 1
p_dm0.set_margin_out(0.3)
p_dm0.set_nact(nact)
p_dm0.set_alt(0.0)
p_dm0.set_thresh(0.3)  # fraction units
# !!!!!!!!!!!!!!!!!!!!!!!!! attention pas autre chose que 0.2 !!!!!!!!!
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(1.0)

p_dm1.set_type(scons.DmType.TT)
p_dm1.set_alt(0.0)
p_dm1.set_unitpervolt(658e-9 / 18.0 * scons.CONST.RAD2ARCSEC)
p_dm1.set_push4imat(0.1)

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
p_centroider0.set_type(scons.CentroiderType.PYR)

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

# p_controller0.set_type("ls")
p_controller0.set_type(scons.ControllerType.GENERIC)

p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(150.0)
p_controller0.set_delay(1)
p_controller0.set_gain(1)


# PyrCentro adhoc
class p_pyrCentro:
    rebin = 5  # Rebin factor of HR pixels
    radiusForce = -1  # Force fitted HR pupil radius to given value; -1 bypasses
    radiusGuard = 0  # Add this number of HR pixels to value fitted
    controlMethod = "ABCD"
    latency = 2
    rMod = 2.0
