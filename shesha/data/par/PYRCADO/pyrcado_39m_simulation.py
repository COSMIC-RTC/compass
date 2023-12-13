import shesha.config as conf
import shesha.constants as scons
import numpy as np

# FOR BOTH SIMU AND RTC STANDALONE
simul_name = "PYRCADO_39m_20190122"

p_loop = conf.ParamLoop()
p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.)  # =1/500
p_loop.set_devices([4, 5, 6, 7])

# FOR THE SIMULATION ONLY
# geom
p_geom = conf.ParamGeom()
p_geom.set_zenithangle(0.)

# tel
p_tel = conf.ParamTel()

#E_ELT PUPIL Rico like
p_tel.set_type_ap(scons.ApertureType.EELT)
p_tel.set_diam(38.542)
p_tel.set_pupangle(0)  #ELT pup rotation in degrees
#p_tel.set_t_spiders(0.51)  #Spider size in meters
p_tel.set_t_spiders(0.)

# atmos
p_atmos = conf.ParamAtmos()

p_atmos.set_r0(0.129)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.])
p_atmos.set_winddir([45.])
p_atmos.set_L0([25.])

# Lambda target(s)
p_targets = [conf.ParamTarget() for _ in range(2)]
Lambda = [0.658, 1.65]
for k, p_target in enumerate(p_targets):
    p_target.set_xpos(0.)
    p_target.set_ypos(0.)
    p_target.set_Lambda(Lambda[k])
    p_target.set_mag(4.)

# wfs
p_wfs0 = conf.ParamWfs()
p_wfss = [p_wfs0]

p_wfs0.set_type(scons.WFSType.PYRHR)
p_wfs0.set_nxsub(92)
# -> 92 sub aps for hexagonal grid of actuators eq. 78 subaps square grid. (pitch = 0.5m)
p_wfs0.set_fracsub(0.0001)  # Minimal illumination fraction
p_wfs0.set_xpos(0.)  # direction of guide star in X (arcsec)
p_wfs0.set_ypos(0.)  # direction of guide star in Y (arcsec)
p_wfs0.set_Lambda(0.658)  # wavelength (microns)
p_wfs0.set_gsmag(13)  # magnitude of guide star
p_wfs0.set_optthroughput(0.28)  # optical transmission
p_wfs0.set_zerop(2.6e10)  # 2.6e10 ph/s/m**2 computed by Rico in R band for MOSAIC
p_wfs0.set_noise(-1.)  # units: electrons/pixel
p_wfs0.set_atmos_seen(1)  # tell if atmos is seen or not
p_wfs0.set_fstop(scons.FieldStopType.SQUARE)
p_wfs0.set_fssize(1.6)  # size of field stop (arcsec)

rMod = 4  # Modulation radius, in lam/D units
nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.) * 4)
p_wfs0.set_pyr_npts(nbPtMod)  # nb pts modu around circle
p_wfs0.set_pyr_ampl(rMod)  # define modulation amplitude
"""
With 52 pixels of margin between 2 edges of pupils on a 240x240 detector and 92 pixels of pupil:
in Compass pupsep is separation between 1 pupil center and Half of detector
pupsep = 52/2+92/2 = 72
"""
p_wfs0.set_pyr_pup_sep(72)  # half pupil separation (center-to-center)

# dm
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dms = [p_dm0, p_dm1]

p_dm0.set_type(scons.DmType.PZT)
p_dm0.set_nact(73)  #73 actuators for a projected M4 pitch of 53cm
p_dm0.set_alt(0.)
p_dm0.set_margin_out(0.6)
p_dm0.set_thresh(0.6)  # fraction units
p_dm0.set_coupling(0.2)  # 0.2 only
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(0.01)
p_dm0.set_type_pattern(scons.PatternType.HEXAM4)
p_dm0.set_influ_type(scons.InfluType.RADIALSCHWARTZ)

p_dm1.set_type(scons.DmType.TT)
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(p_wfs0.Lambda * 1e-6 / p_tel.diam * scons.CONST.RAD2ARCSEC)
# -> Such that we talk to TT mirror in l/D units
p_dm1.set_push4imat(0.005)
