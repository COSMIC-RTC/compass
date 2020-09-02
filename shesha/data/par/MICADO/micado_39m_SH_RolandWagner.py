import shesha.config as ao
import shesha.constants as scons
import numpy as np

simul_name = "micado_39m_SH"

# loop
p_loop = ao.Param_loop()
p_loop.set_niter(30000)
p_loop.set_ittime(1 / 500.)  #=1/500

# geom
p_geom = ao.Param_geom()
p_geom.set_zenithangle(0.)

# tel
p_tel = ao.Param_tel()
p_tel.set_diam(39.0)
p_tel.set_cobs(0.30)

# atmos
p_atmos = ao.Param_atmos()
p_atmos.set_r0(0.129)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.])
p_atmos.set_winddir([45.])
p_atmos.set_L0([1e6])

# 3 targets
p_targets = [ao.Param_target() for _ in range(3)]
Lambda = [1.2, 1.65, 2.2]
k = 0
for p_target in p_targets:
    p_target.set_xpos(0.)
    p_target.set_ypos(0.)
    p_target.set_Lambda(Lambda[k])
    p_target.set_mag(4.)

#wfs
p_wfs0 = ao.Param_wfs(roket=True)
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")
p_wfs0.set_nxsub(78)
p_wfs0.set_npix(6)
p_wfs0.set_pixsize(0.3)
p_wfs0.set_fracsub(0.8)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.5)
p_wfs0.set_gsmag(11)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(1.e11)
p_wfs0.set_noise(0.1)  # in electrons units
p_wfs0.set_atmos_seen(1)
p_wfs0.set_fstop("round")
p_wfs0.set_fssize(2.4)

#dm
p_dm0 = ao.Param_dm()
p_dm1 = ao.Param_dm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type("pzt")
nact = p_wfs0.nxsub + 1

p_dm0.set_nact(nact)
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.3)  # fraction units
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(1.)

p_dm1.set_type("tt")
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(1)
p_dm1.set_push4imat(0.05)

#centroiders
p_centroider0 = ao.Param_centroider()
p_centroiders = [p_centroider0]
p_centroider0.set_nwfs(0)
p_centroider0.set_type("bpcog")
p_centroider0.set_nmax(10)

#controllers
p_controller0 = ao.Param_controller()
p_controllers = [p_controller0]

p_controller0.set_type("ls")
p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(150.)
p_controller0.set_delay(1)
p_controller0.set_gain(0.3)
