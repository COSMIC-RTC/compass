import shesha.constants as scons
import shesha.config as conf
import numpy as np

simul_name = ""
layout = "layoutDeFab"

# LOOP
p_loop = conf.ParamLoop()
p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.0)
p_loop.set_devices([0, 1, 2, 3])

# GEOM
p_geom = conf.ParamGeom()
p_geom.set_zenithangle(0.0)

# TEL
p_tel = conf.ParamTel()
p_tel.set_type_ap("EELT")  # E_ELT PUPIL Rico like
p_tel.set_diam(40.0)  # causes rescaling
p_tel.set_pupangle(0.0)  # ELT pup rotation in degrees
p_tel.set_t_spiders(0.51)  # Spider size in meters

# ATMOS
p_atmos = conf.ParamAtmos()
p_atmos.set_r0(0.089)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.0])
p_atmos.set_winddir([45.0])
p_atmos.set_L0([25.0])  # Not simulated in Yorick?

# 1 LAMBDA TARGET
p_target = conf.ParamTarget()
p_targets = [p_target]
p_target.set_xpos(0)
p_target.set_ypos(0)
p_target.set_Lambda(2.2)
p_target.set_mag(4)

rMod = 3  # Modulation radius, in lam/D units
nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.0) * 4)

# WFS
p_wfs0 = conf.ParamWfs(roket=True)
p_wfss = [p_wfs0]
p_wfs0.set_type("pyrhr")  # type de WFS: "sh", "pyrhr", "pyr"
p_wfs0.set_nxsub(100)  # to get 95 after rescaling
p_wfs0.set_fracsub(0.1)  # Minimal illumination fraction
p_wfs0.set_xpos(0.0)  # direction of guide star in X (arcsec)
p_wfs0.set_ypos(0.0)  # direction of guide star in Y (arcsec)
p_wfs0.set_Lambda(0.7)  # wavelength (microns)
p_wfs0.set_gsmag(11)  # magnitude of guide star
p_wfs0.set_optthroughput(0.28)  # optical transmission
p_wfs0.set_zerop(2.6e10)  # ph/s/m**2 computed by Rico in R band for MOSAIC
p_wfs0.set_noise(0.3)  # electrons/pixel
p_wfs0.set_atmos_seen(1)  # tell if atmos is seen or not
p_wfs0.set_fstop("square")  # shape of field stop, "round", "square"
p_wfs0.set_fssize(1.6)  # size of field stop (arcsec)
p_wfs0.set_pyr_npts(nbPtMod)  # nb pts modu around circle
p_wfs0.set_pyr_ampl(rMod)  # define modulation amplitude
p_wfs0.set_pyr_pup_sep(72)  # half pupil separation (center-to-center)

# DMS
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type(scons.DmType.PZT)
p_dm0.set_thresh(0.0)  # fraction units
p_dm0.set_margin_out(0.6)
p_dm0.set_file_influ_fits("cropped_M4IF.fits")
p_dm0.set_diam_dm(40.0)
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(0.01)
p_dm1.set_type("tt")
p_dm1.set_alt(0.0)
p_dm1.set_unitpervolt(1)
p_dm1.set_push4imat(0.005)

# CENTROIDERS
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]
p_centroider0.set_nwfs(0)
p_centroider0.set_type("maskedpix")

# CONTROLLERS
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]
p_controller0.set_type("generic")  # V(k) = a.E.V(k-1) + g.R.m(k)
p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(150.0)
p_controller0.set_delay(1)  # loop delay. "0 = 1 frame delay".
p_controller0.set_gain(1)
