#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday 7th of February 2020

@author: nour
"""

import shesha.config as conf
import numpy as np

simul_name = "sphere+2stages"
layout = "sphere+2Stages"
# loop
p_loop = conf.ParamLoop()

p_loop.set_niter(5000)  # number of loop iterations
p_loop.set_ittime(1.0 / 3000.0)  # =1/2000 - assuming loop at 2kHz
p_loop.set_devices([0, 1, 2, 3])  # ????
# geom
p_geom = conf.ParamGeom()

p_geom.set_zenithangle(0.0)

# tel
p_tel = conf.ParamTel()

p_tel.set_diam(8.0)  # Subaru diameter
p_tel.set_cobs(0.12)  # TBC (central obstruction)

# atmos
# here we simulate the first stage of correction of ao188
p_atmos = conf.ParamAtmos()

p_atmos.set_r0(0.12)  # Fried parameters @ 500 nm
p_atmos.set_nscreens(1)  # Number of layers
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([15.0])
p_atmos.set_winddir([45])
p_atmos.set_L0([25])  # in meters. here we simulate ao188's precorrection. Layers outer scale

# target
p_target = conf.ParamTarget()
p_targets = [p_target]
p_target.set_xpos(0.0)
p_target.set_ypos(0.0)
p_target.set_Lambda(1.65)
p_target.set_mag(10.0)

# wfs
p_wfs0 = conf.ParamWfs(roket=True)

p_wfs0.set_type("pyrhr")
p_wfs0.set_nxsub(
    50
)  # TBC Number of pixels along the pupil diameter, NB. need more subaperture than nactu.
p_wfs0.set_fssize(2)  # Size of the field stop
p_wfs0.set_fracsub(0.0001)  # was 0.8 before Vincent
p_wfs0.set_xpos(0.0)
p_wfs0.set_ypos(0.0)
p_wfs0.set_Lambda(1.2)  # pyramid wavelength
p_wfs0.set_gsmag(5.0)  # Guide star magnitude
p_wfs0.set_optthroughput(0.5)  # Optiical throughput coefficient
p_wfs0.set_zerop(1.0e11)
p_wfs0.set_noise(-1)
p_wfs0.set_fstop("round")
rMod = 20.0  # Modulation radius, in lam/D units
nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.0) * 4)
p_wfs0.set_pyr_npts(nbPtMod)  # Number of modulation point along the circle
p_wfs0.set_pyr_ampl(rMod)  # Pyramid modulation amplitude (pyramid only)
p_wfs0.set_pyr_pup_sep(p_wfs0.nxsub)  # separation between the 4 images of the pyramid
p_wfs0.set_atmos_seen(1)  # If False, the WFS don’t see the atmosphere layers
p_wfs0.set_dms_seen(np.array([0, 2]))  # If False, the WFS don’t see the atmosphere layers


# wfs
p_wfs1 = conf.ParamWfs(roket=True)

p_wfs1.set_type("pyrhr")
p_wfs1.set_nxsub(
    50
)  # TBC Number of pixels along the pupil diameter, NB. need more subaperture than nactu.
p_wfs1.set_fssize(1.5)  # Size of the field stop
p_wfs1.set_fracsub(0.0001)  # was 0.8 before Vincent
p_wfs1.set_xpos(0.0)
p_wfs1.set_ypos(0.0)
p_wfs1.set_Lambda(1.2)  # pyramid wavelength
p_wfs1.set_gsmag(5.0)  # Guide star magnitude
p_wfs1.set_optthroughput(0.5)  # Optiical throughput coefficient
p_wfs1.set_zerop(1.0e11)
p_wfs1.set_noise(-1)
p_wfs1.set_fstop("round")
rMod = 3.0  # Modulation radius, in lam/D units
nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.0) * 4)
p_wfs1.set_pyr_npts(nbPtMod)  # Number of modulation point along the circle
p_wfs1.set_pyr_ampl(rMod)  # Pyramid modulation amplitude (pyramid only)
p_wfs1.set_pyr_pup_sep(50)  # separation between the 4 images of the pyramid
p_wfs1.set_atmos_seen(1)  # If False, the WFS don’t see the atmosphere layers


p_wfss = [p_wfs0, p_wfs1]


# dm
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dm2 = conf.ParamDm()

p_dm0.set_type("pzt")
# nact = p_wfs0.nxsub + 1
nact = 40
p_dm0.set_nact(nact)
p_dm0.set_alt(0.0)  # Layers altitudes
p_dm0.set_thresh(
    0.3
)  # Threshold on response for selection of valid actuators. Expressed in fraction of the maximal response
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(1.0)
p_dm0.set_push4imat(1.0e-3)  # Nominal voltage for imat = integration matrix = response matrix
p_dm0.set_margin_out(0.3)  # pour adapter la taille de la pupille du DM a celle du WFS


p_dm1.set_type("pzt")
# nact = p_wfs0.nxsub + 1
nact = 32
p_dm1.set_nact(nact)
p_dm1.set_alt(0.0)  # Layers altitudes
p_dm1.set_thresh(
    0.3
)  # Threshold on response for selection of valid actuators. Expressed in fraction of the maximal response
p_dm1.set_coupling(0.2)
p_dm1.set_unitpervolt(1.0)
p_dm1.set_push4imat(1.0e-3)  # Nominal voltage for imat = integration matrix = response matrix
p_dm1.set_margin_out(0.3)  # pour adapter la taille de la pupille du DM a celle du WFS


p_dm2.set_type("tt")
p_dm2.set_alt(0.0)
p_dm2.set_unitpervolt(1.0)  # Influence function sensitivity
p_dm2.set_push4imat(0.01)

p_dms = [p_dm0, p_dm1, p_dm2]

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroider1 = conf.ParamCentroider()
p_centroiders = [p_centroider0, p_centroider1]

p_centroider0.set_nwfs(0)
p_centroider0.set_type("pyr")
p_centroider1.set_nwfs(1)
p_centroider1.set_type("pyr")
# p_centroider0.set_type("corr")
# p_centroider0.set_type_fct("model")

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("generic")  # ls (classic easy simple) or generic
p_controller0.set_nwfs([0, 1])
p_controller0.set_ndm([0, 1, 2])
p_controller0.set_maxcond(5.0)  # what determines the number of modes to be filtered
p_controller0.set_delay(1)
p_controller0.set_gain(0.4)
# p_controller0.set_nstates(6)
