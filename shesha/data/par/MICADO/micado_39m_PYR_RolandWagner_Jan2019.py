import shesha.config as ao
import shesha.constants as scons
import numpy as np

simul_name = ""
layout = "layoutDeFab"

# loop
p_loop = ao.Param_loop()
p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.)  # =1/500Hz = 2ms
p_loop.set_devices([0, 1, 2, 3, 4, 5, 6, 7])

# geom
p_geom = ao.Param_geom()
p_geom.set_zenithangle(0.)

# tel
p_tel = ao.Param_tel()
p_tel.set_type_ap("EELT")
p_tel.set_diam(39)
p_tel.set_pupangle(0)  #ELT pup rotation in degrees
p_tel.set_t_spiders(0.51)  #Spider size in meters

# atmos params
# Q1=0.214
# Q2=0.163
# median=0.144
# Q3=0.127
# Q4=0.089
p_atmos = ao.Param_atmos()
p_atmos.set_r0(0.144)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([9.13])
p_atmos.set_winddir([45.])
p_atmos.set_L0([25.])

# 1 Lambda targets
p_target = ao.Param_target()
p_targets = [p_target]
p_target.set_xpos(0)
p_target.set_ypos(0)
p_target.set_Lambda(2.2)
p_target.set_mag(4)

# wfs params
p_wfs0 = ao.Param_wfs(roket=True)
#p_wfs0= ao.Param_wfs()
p_wfss = [p_wfs0]
p_wfs0.set_type("pyrhr")  # type of WFS: "sh", "pyrhr", "pyr"
p_wfs0.set_nxsub(92)  # 92 sub aps
p_wfs0.set_fracsub(0.1)  # Minimal illumination fraction
p_wfs0.set_xpos(0.)  # direction of guide star in X (arcsec)
p_wfs0.set_ypos(0.)  # direction of guide star in Y (arcsec)
p_wfs0.set_Lambda(0.7)  # wfs wavelength (microns)
p_wfs0.set_gsmag(15)  # magnitude of guide star
p_wfs0.set_optthroughput(0.28)  # optical transmission as computed during PDR
p_wfs0.set_zerop(2.6e10)  # 2.6e10 ph/s/m**2 computed in R band for MOSAIC
p_wfs0.set_noise(0.3)  # units: electrons/pixel OCAM =0.3e-
p_wfs0.set_atmos_seen(1)  # tell if atmos is seen or not
p_wfs0.set_fstop("square")  # shape of field stop, "round", "square"
p_wfs0.set_fssize(1.6)  # size of field stop (arcsec)
rMod = 3  # Modulation radius, in lam/D units
nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.) * 4)
p_wfs0.set_pyr_npts(nbPtMod)  # nb pts modu around circle
p_wfs0.set_pyr_ampl(rMod)  # define modulation amplitude
"""
With 52 pixels of margin between 2 edges of pupils on a 240x240 detector and 92 pixels of pupil:
in Compass pupsep is separation between 1 pupil center and Half of detector
pupsep = 52/2+92/2 = 72
"""
p_wfs0.set_pyr_pup_sep(72)  # half pupil separation (center-to-center)

# dms parameters
p_dm0 = ao.Param_dm()
p_dm1 = ao.Param_dm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type(scons.DmType.PZT)
nact = p_wfs0.nxsub + 1

p_dm0.set_nact(
        73
)  #73 actuators across the diameter for a projected M4 pitch of 53cm (computed from ESO data package)
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.8)  # fraction units (> = less)
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(0.01)
p_dm0.set_type_pattern("hexaM4")
p_dm0.set_influ_type("radialSchwartz")

p_dm1.set_type("tt")
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(1)
p_dm1.set_push4imat(0.005 * 2)

# centroiders params
p_centroider0 = ao.Param_centroider()
p_centroiders = [p_centroider0]
p_centroider0.set_nwfs(0)
p_centroider0.set_type(scons.CentroiderType.PYR)

# controllers params
p_controller0 = ao.Param_controller()
p_controllers = [p_controller0]

#p_controller0.set_type("ls")     # V(k) = V(k-1) + g.R.m(k)
p_controller0.set_type("generic")  # V(k) = a.E.V(k-1) + g.R.m(k)

p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(150.)
p_controller0.set_delay(
        1
)  # loop delay. without integration time "0 = 1 frame delay in total". "1 = 2 frames delay in total".
p_controller0.set_gain(0.7)
