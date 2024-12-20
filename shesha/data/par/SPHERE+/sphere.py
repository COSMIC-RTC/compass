#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday 16th of September 2022

@author: Clementine Bechet for the SPHERE+ simulation group
"""

# /!\ : This symbol marks the parameters that must not be changed.
# /?\ : This symbol marks the parameters that must be chosen among a list.

import shesha.config as conf
import numpy as np

simul_name = "sphere"
layout = "layoutDeFab_SH"  # Reloads custom display layout from layoutDeFab_SH.area

# loop
p_loop = conf.ParamLoop()
p_loop.set_niter(5000)  #     number of loops
# /?\ second stage frequency
# p_loop.set_ittime(1./3000.)     # second loop at 3 kHz
p_loop.set_ittime(1.0 / 3000.0)  # second loop at 2 kHz
# p_loop.set_ittime(1./1000.)   # second loop at 1 kHz
p_loop.set_devices([0, 1, 2, 3])

# geom
p_geom = conf.ParamGeom()
p_geom.set_pupdiam(400)
p_geom.set_zenithangle(0.0)

# tel
p_tel = conf.ParamTel()
p_tel.set_diam(8.0)  # /!\  VLT diameter
p_tel.set_cobs(0.14)  # /!\  central obstruction
p_tel.set_type_ap("VLT")  # /!\  VLT pupil
p_tel.set_t_spiders(0.00625)  # /!\  spider width = 5 cm

# atmos
p_atmos = conf.ParamAtmos()
p_atmos.set_r0(0.14)  #     Fried parameters @ 500 nm
p_atmos.set_nscreens(1)  # /!\ Number of layers
p_atmos.set_frac([1.0])  # /!\ Fraction of atmosphere (100% = 1)
p_atmos.set_alt([0.0])  # /!\ Altitude(s) in meters
p_atmos.set_windspeed([8])  #     Wind speed of layer(s) in m/s
p_atmos.set_winddir([45])  # /!\ Wind direction in degrees
p_atmos.set_L0([25])  #     Outer scale in meters

# target
p_target = conf.ParamTarget()
p_targets = [p_target]

p_target.set_xpos(0.0)  # /!\ On axis
p_target.set_ypos(0.0)  # /!\ On axis
p_target.set_Lambda(1.65)  # /!\ H Band
p_target.set_mag(6.0)  # /!\

# wfs
p_wfs0 = conf.ParamWfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")  # /!\ Shack-Hartmann
p_wfs0.set_nxsub(40)  # /!\ nb of sub-apertures.
p_wfs0.set_npix(6)  # /!\ nb of pixels / sub-aperture.
p_wfs0.set_pixsize(0.36)  # /!\ Shannon at 700nm. No exact reference found
p_wfs0.set_fracsub(0.5)  # /!\ Select 1240 subapertures.
p_wfs0.set_xpos(0.0)  # /!\ On axis
p_wfs0.set_ypos(0.0)  # /!\ On axis
p_wfs0.set_Lambda(0.7)  # /!\ SAXO SH bandwidth : [475, 900] nm
# p_wfs0.set_gsmag(6.)
p_wfs0.set_gsmag(6.0 + 2.5 * np.log10(3 / 2))  # at 2 kHz
# p_wfs0.set_gsmag(6. + 2.5 * np.log10(3 / 1))  # at 1 kHz
p_wfs0.set_optthroughput(0.5)  # still unknown
p_wfs0.set_zerop(1e11)  # zero point for guide star magnitude
p_wfs0.set_noise(0.1)  # EMCCD with < 0.1e- RON
p_wfs0.set_atmos_seen(1)  # /!\
p_wfs0.set_fstop("square")  # /!\
# Choose one spatial filter or none.
# p_wfs0.set_fssize(0.79412)   # 1.1*lambda/dSubap
# p_wfs0.set_fssize(0.8663)    # 1.2*lambda/dSubap
# p_wfs0.set_fssize(0.9385)    # 1.3*lambda/dSubap
# p_wfs0.set_fssize(1.0107)    # 1.4*lambda/dSubap
p_wfs0.set_fssize(1.0829)  # 1.5*lambda/dSubap
# p_wfs0.set_fssize(1.227275)  # 1.7*lambda/dSubap
# p_wfs0.set_fssize(1.44385)   # 2*lambda/dSubap

# dm
p_dm0 = conf.ParamDm()  # /!\
p_dm1 = conf.ParamDm()  # /!\
p_dms = [p_dm0, p_dm1]  # /!\

p_dm0.set_type("pzt")  # /!\
p_dm0.set_thresh(-0.5)  # /!\ to get the SAXO 1377 active actuators
p_dm0.set_alt(0.0)  # /!\
p_dm0.set_unitpervolt(1.0)  # /!\
p_dm0.set_push4imat(0.180)  #     to displace ~ half a pixel
p_dm0.set_file_influ_fits("SAXO_HODM.fits")

# tip-tilt
p_dm1.set_type("tt")  # /!\
p_dm1.set_alt(0.0)  # /!\
p_dm1.set_unitpervolt(1.0)  # /!\
p_dm1.set_push4imat(0.18)  #     to displace about half a pixel

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)  # /!\
p_centroider0.set_type("bpcog")  # need to be replaced by WCOG at some point
p_centroider0.set_nmax(9)  # à regarder vite fait (à priori entre 8 et 12)
# p_centroider0.set_type_fct("model")

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("generic")  # /?\ ls (classic easy simple) or generic
# p_controller0.set_type("ls")
p_controller0.set_nwfs([0])  # /!\
p_controller0.set_ndm([0, 1])  # /!\
p_controller0.set_maxcond(1500)  #     determines the nb of modes to be filtered
p_controller0.set_delay(1.15 * 1000 / 1380)  # /!\ same delay in ms as in saxo.py
p_controller0.set_gain(0.3)
p_controller0.set_calpix_name("compass1_calPix")
p_controller0.set_loopdata_name("compass1_loopData")


"""
p_coronos = [conf.Param_corono()]

p_coronos[0].set_type("SPHERE_APLC")
p_coronos[0].set_wavelength_0(1.667)
p_coronos[0].set_delta_wav(0.054)
p_coronos[0].set_nb_wav(3)

"""
"""
CORONO PARFAIT

p_corono.set_type("perfect")  # coronagraph type : "perfect", "SPHERE_APLC", "custom"
p_corono.set_wavelength_0(1.667)  # coronagraph wavelength in micron
p_corono.set_dim_image(250)     # size of the science image in pixel
p_corono.set_image_sampling(2)  # number of pixel in lambda/D (2 : Shannon sampling)
# p_corono.set_delta_wav(0.054)  # spectral bandwidth in micron (for coronagraph images only)
# p_corono.set_nb_wav(3)        # number of simulated wavelengths in the coronagraph module

# p_corono.set_image_sampling(1667e-9 / 8 * 180 / np.pi * 3600 * 1000 / 12.25)  # to match SPHERE IRDIS pixel scale
"""

"""
CORONO APLC SPHERE

p_corono.set_type("SPHERE_APLC")
p_corono.set_wavelength_0(1.667)
"""

"""
CORONO MANUEL
"""

"""
p_corono = conf.Param_corono()

p_corono.set_type("custom")
p_corono.set_wavelength_0(1.667)
p_corono.set_apodizer_name('SPHERE_APLC_apodizer_APO1')  # 'SPHERE_APLC_apodizer_APO1' or path to fits file of size (pupdiam, pupdiam)
p_corono.set_focal_plane_mask_name('classical_Lyot')  # 'classical_Lyot', "SPHERE_APLC_fpm_ALC1" or 2 or 3
                                                      # or path to fits file of user defined size [dimx, dimy, nb_wav]
                                                      # or [dimx, dimy] if same for all wavelength
p_corono.set_lyot_fpm_radius(2.15)  # radius of the mask in lambda/D units (required)
# p_corono.set_fpm_sampling()  # lambda/D size in pixel (default = 20.)
                             # optional if classical lyot
p_corono.set_lyot_stop_name('SPHERE_APLC_Lyot_stop')  # 'SPHERE_APLC_Lyot_stop', or path to fits file of size (pupdiam, pupdiam)
                                                      # homogeneous to get_s_pupil() size in COMPASS
p_corono.set_dim_image(200)
p_corono.set_image_sampling(3)

# p_corono.set_babinet_trick() # True or False for enabling Babinet propagation method
"""
