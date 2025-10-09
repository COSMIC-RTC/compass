""" SPHERE first stage
to be used with a twoStagesManager
"""
import shesha.config as conf
import numpy as np
# from data.par import atmosphere

simul_name = "sphere"
layout = "sphere_corono"

# loop
p_loop = conf.ParamLoop()
p_loop.set_devices([0, 1, 2, 3])
first_stage_frequency = 500  # [Hz]
p_loop.set_ittime(1 / first_stage_frequency)

# geom
p_geom = conf.ParamGeom()
p_geom.set_pupdiam(400)
p_geom.set_zenithangle(0.)

# tel
p_tel = conf.ParamTel()
p_tel.set_diam(8.0)            # /!\  VLT diameter
p_tel.set_type_ap("VLT")       # /!\  VLT pupil
p_tel.set_cobs(0.14)           # /!\  central obstruction
p_tel.set_t_spiders(0.00625)   # /!\  spider width = 5 cm

# atmos
p_atmos = conf.ParamAtmos()
p_atmos.set_r0(0.12)
p_atmos.set_nscreens(3)
p_atmos.set_frac([0.50, 0.35, 0.15])
p_atmos.set_alt([0, 4000, 10000])
p_atmos.set_windspeed([15, 15, 35])
p_atmos.set_winddir([0, 20, 180])
p_atmos.set_L0([30, 30, 30])
# atm = atmosphere.TurbulenceProfile()                     # /?\ instantiate 1L or 35L atmosphere
# conditions = 'median'                                    # median profile
# seeing = 1.2                                             # [arcsec] seeing
# tau0 = 3                                                 # [ms] coherence time
# r0 = 500e-9 / seeing * 180 * 3600 / np.pi                # [m] Fried parameter
# veff = 0.314 * r0 / tau0 * 1e3                           # [m/s] effective windspeed
# veff_median = np.sum(atm.wind_velocity**(5 / 3) * atm.fractional['median'])**(3 / 5)
# scaling_wind = veff / veff_median.value                  
# wind_profile = atm.wind_velocity * scaling_wind
# p_atmos.set_r0(r0)                                       # /!\ set r0
# p_atmos.set_nscreens(len(atm.height.value))              # /!\ set Nb of layers
# p_atmos.set_windspeed(wind_profile)                      # /!\ set wind speeds
# p_atmos.set_frac(atm.profile_data[conditions].tolist())  # /!\ set frac Cn2
# p_atmos.set_alt(atm.heights())                           # /!\ set heights
# p_atmos.set_winddir(atm.wind_dir.value)                  # /!\ set wind directions
# p_atmos.set_L0([25.] * len(atm.height.value))            # /!\ outer scale 25m

# target
p_target = conf.ParamTarget()
p_targets = [p_target]

p_target.set_xpos(0.)         # /!\ On axis
p_target.set_ypos(0.)         # /!\ On axis
p_target.set_Lambda(1.65)     # /!\ H Band
p_target.set_mag(6.)          # /!\

# wfs
p_wfs0 = conf.ParamWfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")         # /!\ Shack-Hartmann
p_wfs0.set_nxsub(40)          # /!\ nb of sub-apertures.
p_wfs0.set_npix(6)            # /!\ nb of pixels / sub-aperture.
p_wfs0.set_pixsize(0.36)      # /!\ Shannon at 700nm. No exact reference found
p_wfs0.set_fracsub(0.5)       # /!\ Select 1240 subapertures.
p_wfs0.set_xpos(0.)           # /!\ On axis
p_wfs0.set_ypos(0.)           # /!\ On axis
p_wfs0.set_Lambda(0.7)        # /!\ SAXO SH bandwidth : [475, 900] nm
p_wfs0.set_gsmag(1.5583496)   # /?\ guide star magnitude
p_wfs0.set_optthroughput(0.5) # still unknown
p_wfs0.set_zerop(1e11)        # zero point for guide star magnitude
p_wfs0.set_noise(0.1)         # EMCCD with < 0.1e- RON
p_wfs0.set_atmos_seen(1)      # /!\
# p_wfs0.set_fstop("square")    # /!\
# p_wfs0.set_fssize(0.89)       # [arcsec] small size : 0.82 arcsec
                                # medium : 0.89 arcsec ?
                                # large : 1.07 arcsec ?

# dm
p_dm0 = conf.ParamDm()       # /!\
p_dm1 = conf.ParamDm()       # /!\
p_dms = [p_dm0, p_dm1]        # /!\

p_dm0.set_type("pzt")         # /!\
p_dm0.set_thresh(-0.5)        # /!\ to get the SAXO 1377 active actuators
p_dm0.set_alt(0.)             # /!\
p_dm0.set_unitpervolt(1.)     # /!\
p_dm0.set_push4imat(0.180)    #     to displace ~ half a pixel
dm_data = '/55_v550_twostage/shesha/data/dm-data/'
p_dm0.set_file_influ_fits("SAXO_HODM_gauss_fitSPARTA.fits")

# tip-tilt
p_dm1.set_type("tt")         # /!\
p_dm1.set_alt(0.)            # /!\
p_dm1.set_unitpervolt(1.)    # /!\
p_dm1.set_push4imat(0.18)    #     to displace about half a pixel

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)        # /!\
p_centroider0.set_type("wcog")   # weighted center of gravity
p_centroider0.set_width(2)
p_centroider0.set_thresh(0)

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("generic")   # /?\ ls (classic easy simple) or generic
p_controller0.set_calpix_name("compass1_calPix")
p_controller0.set_loopdata_name("compass1_loopData")
p_controller0.set_nwfs([0])         # /!\
p_controller0.set_ndm([0, 1])       # /!\
p_controller0.set_gain(0.3)
p_controller0.set_delay(1.15 * first_stage_frequency / 1380)  # /!\ delay = 1.15 frame at 1.38 kHz

# coronagraphs
p_corono0 = conf.ParamCoronagraph()
p_coronos = [p_corono0]

p_corono0.set_type("perfect")      # coronagraph type : "perfect", "SPHERE_APLC", "custom"
p_corono0.set_wavelength_0(1.667)  # coronagraph wavelength in micron
p_corono0.set_image_sampling(1667e-9 / 8 * 180 / np.pi * 3600 * 1000 / 12.25) # to match SPHERE IRDIS pixel scale
p_corono0.set_dim_image(160)       # size of the science image in pixel
