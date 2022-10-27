#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday 7th of February 2020

@author: nour
"""

import shesha.config as conf

simul_name = "sphere"
layout = "layoutDeFab_SH" # Reloads custom display layout from layoutDeFab_SH.area

# loop
p_loop = conf.Param_loop()

p_loop.set_niter(5000)          # number of loop iterations
p_loop.set_ittime(1./3000.)     # assuming loop at 1.38 kHz
p_loop.set_devices([0, 1, 2, 3])
# geom
p_geom = conf.Param_geom()
p_geom.set_pupdiam(400)
p_geom.set_zenithangle(0.)

# tel
p_tel = conf.Param_tel()

p_tel.set_diam(8.0)            # VLT diameter
#p_tel.set_cobs(0.14)           # central obstruction
p_tel.set_type_ap("VLT")       # VLT pupil
#p_tel.set_spiders_type("four")
#p_tel.set_t_spiders(0.00625)

# atmos
p_atmos = conf.Param_atmos()

p_atmos.set_r0(0.14)       # Fried parameters @ 500 nm
p_atmos.set_nscreens(1)    # Number of layers
p_atmos.set_frac([1.0])    # Fraction of atmosphere (100% = 1)
p_atmos.set_alt([0.0])     # Altitude(s) in meters
p_atmos.set_windspeed([8]) # wind speed of layer(s) in m/s
p_atmos.set_winddir([45])  # wind direction in degrees
p_atmos.set_L0([25])       # in meters

# target
p_target = conf.Param_target()
p_targets = [p_target]
p_target.set_xpos(0.)
p_target.set_ypos(0.)
p_target.set_Lambda(1.65) # H Band
p_target.set_mag(6.)

# wfs
p_wfs0 = conf.Param_wfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")
p_wfs0.set_nxsub(40)
p_wfs0.set_npix(6)
p_wfs0.set_pixsize(0.3)
p_wfs0.set_fracsub(0.8)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.7)        # SAXO SH bandwidth : [475, 900] nm
p_wfs0.set_gsmag(6.)
p_wfs0.set_optthroughput(0.5) # still unknown
p_wfs0.set_zerop(1e11)        # zero point for guide star magnitude
p_wfs0.set_noise(0.2)         # -1 = No noise
p_wfs0.set_atmos_seen(1)
p_wfs0.set_fstop("round")
p_wfs0.set_fssize(1.1)

# dm (waiting for the custom HODM)
p_dm0 = conf.Param_dm()
p_dm1 = conf.Param_dm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type("pzt")
p_dm0.set_alt(0.)
p_dm0.set_thresh(-0.5)
p_dm0.set_coupling(0.3)
p_dm0.set_unitpervolt(1.)
p_dm0.set_push4imat(0.001)
p_dm0.set_file_influ_fits('SAXO_HODM.fits')

#p_dm0.set_margin_out(0.3)

p_dm1.set_type("tt")
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(1.)
p_dm1.set_push4imat(0.005)
# centroiders
p_centroider0 = conf.Param_centroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
p_centroider0.set_type("bpcog")
p_centroider0.set_nmax(6) # nb of brightest pixels
# p_centroider0.set_type("corr")

# controllers
p_controller0 = conf.Param_controller()
p_controllers = [p_controller0]

p_controller0.set_type("generic") # ls (classic easy simple) or generic
p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(1500)   # what determines the number of modes to be filtered
p_controller0.set_delay(1.5)
p_controller0.set_gain(0.3)
p_controller0.set_calpix_name("compass1_calPix")
p_controller0.set_loopdata_name("compass1_loopData")