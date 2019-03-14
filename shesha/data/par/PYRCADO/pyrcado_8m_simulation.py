import shesha.config as conf
import shesha.constants as scons
import numpy as np

# FOR BOTH SIMU AND RTC STANDALONE
simul_name = "PYRCADO_8m_20190115"

p_loop = conf.Param_loop()
p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.)  # =1/500
p_loop.set_devices([4, 5, 6, 7])

# FOR THE SIMULATION ONLY
# geom
p_geom = conf.Param_geom()
p_geom.set_zenithangle(0.)

# tel
p_tel = conf.Param_tel()
p_tel.set_diam(8.)

# atmos
p_atmos = conf.Param_atmos()

p_atmos.set_r0(0.129)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.])
p_atmos.set_winddir([45.])
p_atmos.set_L0([25.])

# Lambda target(s)
p_targets = [conf.Param_target() for _ in range(2)]
Lambda = [0.658, 1.65]
for k, p_target in enumerate(p_targets):
    p_target.set_xpos(0.)
    p_target.set_ypos(0.)
    p_target.set_Lambda(Lambda[k])
    p_target.set_mag(4.)

# wfs
p_wfs0 = conf.Param_wfs()
p_wfss = [p_wfs0]

p_wfs0.set_type(scons.WFSType.PYRHR)
p_wfs0.set_nxsub(24)  # No more soft binning
p_wfs0.set_fracsub(0.0001)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.658)
p_wfs0.set_gsmag(13)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(2.6e10)  # 2.6e10 ph/s/m**2 computed by Rico in R band for MOSAIC
p_wfs0.set_noise(-1)  # in electrons units
p_wfs0.set_atmos_seen(1)
p_wfs0.set_fstop(scons.FieldStopType.SQUARE)
p_wfs0.set_fssize(1.6)

rMod = 4

p_wfs0.set_pyr_npts(int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.) * 4))
p_wfs0.set_pyr_ampl(rMod)
p_wfs0.set_pyr_pup_sep((20))

# dm
p_dm0 = conf.Param_dm()
p_dm1 = conf.Param_dm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type(scons.DmType.PZT)
nact = 17
p_dm0.set_margin_out(0.3)
p_dm0.set_nact(nact)
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.3)  # fraction units
# !!!!!!!!!!!!!!!!!!!!!!!!! attention pas autre chose que 0.2 !!!!!!!!!
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(1.)

p_dm1.set_type(scons.DmType.TT)
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(p_wfs0.Lambda * 1e-6 / p_tel.diam * scons.CONST.RAD2ARCSEC)
# -> Such that we talk to TT mirror in l/D units
p_dm1.set_push4imat(0.1)
