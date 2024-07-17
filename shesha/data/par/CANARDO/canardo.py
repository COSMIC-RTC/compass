import shesha.config as conf
import numpy as np

simul_name = "CANARDO"
layout = "layoutDeFab"

# loop
p_loop = conf.ParamLoop()

p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.0)  # =1/500
p_loop.set_devices([4, 5, 6, 7])

# geom
p_geom = conf.ParamGeom()
p_geom.set_zenithangle(0.0)

# tel
p_tel = conf.ParamTel()
p_tel.set_diam(4.2)

# atmos
p_atmos = conf.ParamAtmos()

p_atmos.set_r0(0.129)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.0])
p_atmos.set_winddir([45.0])
p_atmos.set_L0([1.0e5])

# target
p_targets = [conf.ParamTarget()]

p_targets[0].set_xpos(0.0)
p_targets[0].set_ypos(0.0)
p_targets[0].set_Lambda(1.65)
p_targets[0].set_mag(4.0)

rMod = 3  # Modulation radius, in lam/D units
nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.0) * 4)
# wfs
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
p_wfs0.set_noise(-1)  # electrons/pixel
p_wfs0.set_atmos_seen(1)  # tell if atmos is seen or not
p_wfs0.set_fstop("square")  # shape of field stop, "round", "square"
p_wfs0.set_fssize(16)  # size of field stop (arcsec)
p_wfs0.set_pyr_npts(nbPtMod)  # nb pts modu around circle
p_wfs0.set_pyr_ampl(rMod)  # define modulation amplitude
p_wfs0.set_pyr_pup_sep(72)  # half pupil separation (center-to-center)

# dm
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type("pzt")
p_dm0.set_file_influ_fits("simuALPAO4KDM.fits")
p_dm0.set_alt(0.0)
p_dm0.set_thresh(0.3)
p_dm0.set_unitpervolt(0.01)
p_dm0.set_push4imat(100.0)
p_dm0.set_diam_dm(0.093)

p_dm1.set_type("tt")
p_dm1.set_alt(0.0)
p_dm1.set_unitpervolt(0.0005)
p_dm1.set_push4imat(10.0)

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
p_centroider0.set_type("pyr")

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

# p_controller0.set_type("ls")     # V(k) = V(k-1) + g.R.m(k)
p_controller0.set_type("generic")  # V(k) = a.E.V(k-1) + g.R.m(k)
# p_controller0.set_type("geo")    # bypass the WFS (direct DM proj)

p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(150.0)
p_controller0.set_delay(1)  # loop delay. "0 = 1 frame delay".
p_controller0.set_gain(0.4)
