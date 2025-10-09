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
import os
# import atmosphere
import sys
sys.path.append(os.path.abspath("/home/micado/saxoplus-reference-cases/scripts"))
atmosphere = __import__("atmosphere")
# from saxoPlusUtils import saxo_compute_delay
# import dictionary_handler as simu       # for dictionary of simulations

# ---- Directories around this script with the DM data are needed
# _my_cases_dir, _my_param_name = os.path.split(os.path.abspath(__file__))
# _my_dm_dir = os.path.abspath(os.path.join(_my_cases_dir, "../data"))

simul_name = "sphere"  # [?]

# FIXME: question "should us remove?" at mark: [?]

# /?\ Choose a first stage loop frequency. (you may create your own pair)
#saxo_frequency_delay_pair = [1380.,1.15] # from Cantalloube et al.
#saxo_frequency_delay_pair = [600., 0.5] # from Cantalloube et al.
#saxo_frequency_delay_pair = [300.,0.25] # from Cantalloube et al.

freq = 1000. # create your own pair here.
saxo_frequency_delay_pair = [freq, 0.84] # your custom pair

# loop
p_loop = conf.ParamLoop()
p_loop.set_niter(5000)                               # number of loops [?]
p_loop.set_ittime(1./saxo_frequency_delay_pair[0])   # first stage frequency [?]
p_loop.set_devices([0, 1, 2, 3])  # FIXME: only one device?

# geom
p_geom = conf.ParamGeom()
p_geom.set_pupdiam(400)        # /!\ value to get WFS pixsize = 0.338 (SAXO truth = 0.36)
p_geom.set_zenithangle(0.)

# tel
p_tel = conf.ParamTel()
p_tel.set_diam(8.0)            # /!\  VLT diameter
p_tel.set_cobs(0.14)           # /!\  central obstruction
p_tel.set_type_ap("VLT")       # /!\  VLT pupil
p_tel.set_t_spiders(0.00625)   # /!\  spider width = 5 cm

# 35 layer atmosphere by ESO
p_atmos = conf.ParamAtmos()
atm = atmosphere.TurbulenceProfile(single_layer=False)  # /?\ instantiate 1L or 35L atmosphere
profileCondition = 'median'    # /?\ 'median', 'Q1', 'Q2', 'Q3' or 'Q4'

# FIXME: remove r0 (set in run_script)? [?]
# r0 = (atm.r0[atm.profile_conditions.index(profileCondition)]).value # /?\ take ESO r0 or..
#r0 = 0.14                                                          # /?\ set your r0 value
r0=0.26
n_layers = 35
p_atmos.set_nscreens(n_layers)                         # /!\ set Nb of layers
p_atmos.set_r0(r0)                                                  # /!\ set r0 [?]
p_atmos.set_frac(atm.profile_data[profileCondition].tolist())       # /!\ set frac Cn2
p_atmos.set_alt(atm.height.value)                                   # /!\ set heights
p_atmos.set_windspeed(atm.wind_velocity.value)                      # /!\ set wind speeds
p_atmos.set_winddir(atm.wind_dir.value)                             # /!\ set wind directions
p_atmos.set_L0([25.] * len(atm.height.value))                       # /!\ outer scale 25m

coherenceTime = 9.
windSpeed = np.copy(p_atmos.windspeed) # [m/s]
equivSpeed = np.sum(p_atmos.frac*(np.abs(windSpeed))**(5./3))**(3./5)
tau0Seeing = 0.31 * r0 / equivSpeed 
speedFactor = coherenceTime * 1e-3 / tau0Seeing
p_atmos.set_windspeed(atm.wind_velocity.value/speedFactor)   

# target
p_target = conf.ParamTarget()
p_targets = [p_target]
p_target.set_xpos(0.)         # /!\ On axis
p_target.set_ypos(0.)         # /!\ On axis
p_target.set_Lambda(1.65) # /!\ H Band
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
p_wfs0.set_Lambda(0.690)      # /!\ SAXO SH bandwidth : [475, 900] nm

# gsmag, zerop, and optthroughput are not used. Photon flux is fixed, the 
# corresponding magnitude is computed taking zerop and optthroughput into 
# account and set, and compass will get the photons back from the magnitude 
# taking the same zerop and optthroughput into account.
p_wfs0.set_gsmag(15.58)          # Not used
p_wfs0.set_zerop(1e11)        # Not used (zero point for guide star magnitude)
p_wfs0.set_optthroughput(0.5) # Not used

# Standard deviation of electronic noise.
# If <0, no noise, if =0 photon noise enabled, if >0 add electronic noise
p_wfs0.set_noise(0.058)         # EMCCD with < 0.1e- RON

p_wfs0.set_atmos_seen(1)      # /!\
p_wfs0.set_fstop("square")    # /!\

# spatial filter
# The size is computed in arcsec for lambda = 690 nm.
# Real possible values are: "open", "small" (1.15*lambda/dSubap), "medium"
# (1.25*lambda/dSubap), or large (1.5*lambda/dSubap).
# Default value is medium.
                              # Choose one spatial filter or none.
#p_wfs0.set_fssize(0.8184)     # 1.15*lambda/dSubap, called SMALL in SAXO
p_wfs0.set_fssize(0.8895)     # 1.25*lambda/dSubap, called MEDIUM in SAXO
#p_wfs0.set_fssize(1.0674)     # 1.5*lambda/dSubap, called LARGE in SAXO

# dm
p_dm0 = conf.ParamDm()       # /!\
p_dm1 = conf.ParamDm()       # /!\
p_dms = [p_dm0, p_dm1]        # /!\

p_dm0.set_type("pzt")         # /!\
p_dm0.set_thresh(-0.5)        # /!\ to get the SAXO 1377 active actuators
p_dm0.set_alt(0.)             # /!\
p_dm0.set_unitpervolt(1.)     # /!\
p_dm0.set_push4imat(0.180)    # to displace about half a pixel
dm_data = 'data/dm-and-cmat-data/saxo/'
p_dm0.set_file_influ_fits('SAXO_HODM_gauss_fitSPARTA.fits') # /!\ to use a custom SAXO HO DM

# tip-tilt
p_dm1.set_type("tt")         # /!\
p_dm1.set_alt(0.)            # /!\
p_dm1.set_unitpervolt(1.)    # /!\
p_dm1.set_push4imat(0.18)    # to displace about half a pixel

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)     # /!\
p_centroider0.set_type("wcog")
p_centroider0.set_width(2) # size of the diffraction limit spot
p_centroider0.set_thresh(0)
# p_centroider0.set_type_fct("model")

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("generic")   # /?\ ls (classic easy simple) or generic
p_controller0.set_nwfs([0])         # /!\
p_controller0.set_ndm([0, 1])       # /!\
p_controller0.set_delay(saxo_frequency_delay_pair[1]) # /!\ see at the head of the file. [?]
p_controller0.set_gain(0.2)

# coronagraph
p_corono0 = conf.ParamCoronagraph()
p_coronos = [p_corono0]

p_corono0.set_type('perfect')
p_corono0.set_dim_image(200)    # size of the science image in pixel
p_corono0.set_wavelength_0(1.65)
                                # coronagraph wavelength in microns
                                # 1.667 μm = central wavelength of H3 IRDIS filter
                                # this value is just an example, not a fixed parameter

# this calculation is done to match IRDIS pixel scale, according to wavelength_0
IRDIS_pixel_scale = 12.25 / 1000 / 3600  # [rad] IRDIS pixel scale
VLT_diameter = 8                         # [m]
wavelength_0 = p_corono0.get_wavelength_0() * 1e-6
image_sampling = np.rad2deg(wavelength_0 / VLT_diameter) / IRDIS_pixel_scale
p_corono0.set_image_sampling(image_sampling)  # size of lambda/D in pixel

# for a polychromatic coronagraph (optional), set these two parameters:
# p_corono0.set_nb_wav(3)           # number of simulated wavelength
# p_corono0.set_delta_wav(0.054)    # spectral bandwidth in micron
                                    # 0.054 μm = bandwidth of H3 IRDIS filter
                                    # this value is just an example, not a fixed parameter
