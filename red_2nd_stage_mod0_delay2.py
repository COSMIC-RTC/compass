""" SPHERE second stage
to be used with a twoStagesManager
"""
import shesha.config as conf
import numpy as np

simul_name = "sphere+"
layout = "sphere+_corono"

# loop
p_loop = conf.ParamLoop()
p_loop.set_devices([0, 1, 2, 3])
second_stage_frequency = 2000  # [Hz]
p_loop.set_ittime(1 / second_stage_frequency)

# geom
p_geom = conf.ParamGeom()
p_geom.set_zenithangle(0.)

# tel
p_tel = conf.ParamTel()
p_tel.set_diam(8.0)            # /!\  VLT diameter
p_tel.set_cobs(0.14)           # /!\  central obstruction
p_tel.set_type_ap("VLT")       # /!\  VLT pupil
p_tel.set_t_spiders(0.00625)   # /!\  spider width = 5 cm

# atmos
p_atmos = conf.ParamAtmos()
p_atmos.set_r0(0.14)          #     Fried parameters @ 500 nm
p_atmos.set_nscreens(1)       # /!\ Number of layers
p_atmos.set_frac([1.0])       # /!\ Fraction of atmosphere (100% = 1)
p_atmos.set_alt([0.0])        # /!\ Altitude(s) in meters
p_atmos.set_windspeed([8])    #     Wind speed of layer(s) in m/s
p_atmos.set_winddir([45])     # /!\ Wind direction in degrees
p_atmos.set_L0([25])          #     Outer scale in meters

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

p_wfs0.set_type("pyrhr")        # /!\ pyramid
p_wfs0.set_nxsub(50)            #     number of pixels
p_wfs0.set_fracsub(0.0001)      #     threshold on illumination fraction for valid pixel
p_wfs0.set_Lambda(1.2)          #     wavelength
p_wfs0.set_gsmag(11.484375)     #     guide star magnitude
p_wfs0.set_zerop(1.e11)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_noise(0.8)           #     readout noise
p_wfs0.set_xpos(0.)             # /!\ On axis
p_wfs0.set_ypos(0.)             # /!\ On axis
rMod = 0.                       # Modulation radius, in lam/D units
p_wfs0.set_pyr_ampl(rMod)
if rMod == 0:
    nbPtMod = 1
else:
    nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.) * 4)
p_wfs0.set_pyr_npts(nbPtMod) 
p_wfs0.set_pyr_pup_sep(p_wfs0.nxsub) # separation between the 4 images of the pyramid 
p_wfs0.set_fstop("round")
p_wfs0.set_fssize(1.5)          # Size of the field stop
p_wfs0.set_atmos_seen(1)        # /!\

# dm
p_dm0 = conf.ParamDm()
p_dms = [p_dm0]

p_dm0.set_type("pzt")           # /!\
p_dm0.set_thresh(-200)          # /!\ to get all Boston actuators
p_dm0.set_alt(0.)               # /!\
p_dm0.set_unitpervolt(1.)       # /!\
p_dm0.set_push4imat(1.0e-3)
p_dm0.set_file_influ_fits("Boston28x28_SquareSchwartz.fits")
p_dm0.set_diam_dm((28 - 3) * 450e-6)

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)           # /!\
p_centroider0.set_type("maskedpix")

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("generic")   # /?\ generic => must do manual imat.
p_controller0.set_calpix_name("compass2_calPix")
p_controller0.set_loopdata_name("compass2_loopData")
p_controller0.set_nwfs([0])         # /!\
p_controller0.set_ndm([0])          # /!\
p_controller0.set_gain(0.3)
p_controller0.set_delay(2)          # 1 frame delay

p_hrtc = conf.ParamHrtc()

p_hrtc.set_frame_shm("UnscrambledWFSframeInterface")
p_hrtc.set_com_shm("IntegratedCommands")
p_hrtc.set_framesize(240)

p_corono0 = conf.ParamCoronagraph()
p_coronos = [p_corono0]

p_corono0.set_type('perfect')      # coronagraph type : 'perfect', 'SPHERE_APLC', 'custom'
p_corono0.set_wavelength_0(1.667)  # coronagraph central wavelength in micron
p_corono0.set_image_sampling(1.667e-6 / 8 * 180 / np.pi * 3600 * 1000 / 12.25)
p_corono0.set_dim_image(160)
