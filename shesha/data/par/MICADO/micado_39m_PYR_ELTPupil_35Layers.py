import shesha.config as ao
import shesha.constants as scons
import numpy as np
simul_name = ""
layout = "layoutDeFab"

# loop
p_loop = ao.Param_loop()

p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.)  # =1/500
p_loop.set_devices([0, 1, 2, 3])
# geom
p_geom = ao.Param_geom()
p_geom.set_zenithangle(0.)

# tel
p_tel = ao.Param_tel()
p_tel.set_diam(38.542)
p_tel.set_cobs(0.28)

#E_ELT PUPIL

p_tel.set_type_ap("EELT-Nominal")
p_tel.set_spiders_type("six")
p_tel.set_pupangle(0)
p_tel.set_t_spiders(0)
"""
p_tel.set_nbrmissing(7)
p_tel.set_referr(0.01)
p_tel.set_std_tt(0.050) # microns
p_tel.set_std_piston(0.050) # microns
"""

# -----------------------------
# 1 layer profile
# -----------------------------

r01L = 0.129
alt1L = [0.0]
frac1L = [1.0]
wind1L = [10.0]

# -----------------------------
# 35 layers profile
# -----------------------------

r0med = 0.144

r0Q1 = 0.2147
r0Q2 = 0.1633
r0Q3 = 0.1275
r0Q4 = 0.089

zenithAngle = 30.
altESO = np.array([
        30, 90, 150, 200, 245, 300, 390, 600, 1130, 1880, 2630, 3500, 4500, 5500, 6500,
        7500, 8500, 9500, 10500, 11500, 12500, 13500, 14500, 15500, 16500, 17500, 18500,
        19500, 20500, 21500, 22500, 23500, 24500, 25500, 26500
]) / np.cos(zenithAngle * 2 * np.pi / 360)
altESO = altESO.astype(int)

fracmed = [
        24.2, 12, 9.68, 5.9, 4.73, 4.73, 4.73, 4.73, 3.99, 3.24, 1.62, 2.6, 1.56, 1.04,
        1, 1.2, 0.4, 1.4, 1.3, 0.7, 1.6, 2.59, 1.9, 0.99, 0.62, 0.4, 0.25, 0.22, 0.19,
        0.14, 0.11, 0.06, 0.09, 0.05, 0.04
]
fracQ1 = [
        22.6, 11.2, 10.1, 6.4, 4.15, 4.15, 4.15, 4.15, 3.1, 2.26, 1.13, 2.21, 1.33, 0.88,
        1.47, 1.77, 0.59, 2.06, 1.92, 1.03, 2.3, 3.75, 2.76, 1.43, 0.89, 0.58, 0.36,
        0.31, 0.27, 0.2, 0.16, 0.09, 0.12, 0.07, 0.06
]
fracQ2 = [
        25.1, 11.6, 9.57, 5.84, 3.7, 3.7, 3.7, 3.7, 3.25, 3.47, 1.74, 3, 1.8, 1.2, 1.3,
        1.56, 0.52, 1.82, 1.7, 0.91, 1.87, 3.03, 2.23, 1.15, 0.72, 0.47, 0.3, 0.25, 0.22,
        0.16, 0.13, 0.07, 0.11, 0.06, 0.05
]
fracQ3 = [
        25.5, 11.9, 9.32, 5.57, 4.5, 4.5, 4.5, 4.5, 4.19, 4.04, 2.02, 3.04, 1.82, 1.21,
        0.86, 1.03, 0.34, 1.2, 1.11, 0.6, 1.43, 2.31, 1.7, 0.88, 0.55, 0.36, 0.22, 0.19,
        0.17, 0.12, 0.1, 0.06, 0.08, 0.04, 0.04
]
fracQ4 = [
        23.6, 13.1, 9.81, 5.77, 6.58, 6.58, 6.58, 6.58, 5.4, 3.2, 1.6, 2.18, 1.31, 0.87,
        0.37, 0.45, 0.15, 0.52, 0.49, 0.26, 0.8, 1.29, 0.95, 0.49, 0.31, 0.2, 0.12, 0.1,
        0.09, 0.07, 0.06, 0.03, 0.05, 0.02, 0.02
]

windESO = [
        5.5, 5.5, 5.1, 5.5, 5.6, 5.7, 5.8, 6, 6.5, 7, 7.5, 8.5, 9.5, 11.5, 17.5, 23, 26,
        29, 32, 27, 22, 14.5, 9.5, 6.3, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10
]

# atmos
p_atmos = ao.Param_atmos()
"""
r0 = r01L
alt = alt1L
frac = frac1L
wind = wind1L
"""

r0 = r0Q3
alt = altESO
frac = fracQ3
wind = windESO

nbLayers = len(alt)
# p_atmos.set_r0(0.129)
p_atmos.set_r0(r0)
p_atmos.set_nscreens(nbLayers)
p_atmos.set_frac(frac)
p_atmos.set_alt(alt)
p_atmos.set_windspeed(wind)
p_atmos.set_winddir([45.] * nbLayers)
p_atmos.set_L0([25.] * nbLayers)  # Not simulated in Yorick?

# target
#p_target = ao.Param_target()
#p_target.set_nTargets(1)
#p_target.set_xpos([0])
#p_target.set_ypos([0.])
#p_target.set_Lambda([1.65])
#p_target.set_mag([4.])

# 3 Lambda targets
#p_target=ao.Param_target()
#p_target.set_nTargets(3)
#p_target.set_xpos([0, 0, 0])
#p_target.set_ypos([0, 0, 0])
#p_target.set_Lambda([1.2, 1.65, 2.2])
#p_target.set_mag([4, 4., 4])

# 1 Lambda targets
p_target = ao.Param_target()
p_targets = [p_target]
p_target.set_xpos(0)
p_target.set_ypos(0)
p_target.set_Lambda(2.2)
p_target.set_mag(4)
"""
p_target=ao.Param_target()
p_target.set_ntargets(11)
p_target.set_xpos([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
p_target.set_ypos([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
p_target.set_Lambda([2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2])
p_target.set_mag([4, 4., 4, 4., 4, 4., 4, 4., 4, 4., 4.])
# wfs
"""
p_wfs0 = ao.Param_wfs(roket=True)
#p_wfs0= ao.Param_wfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("pyrhr")

p_wfs0.set_nxsub(
        92
)  # 92 sub aps for hexagonal grid of actuators eq. 78 subaps square grid. (pitch = 0.5m)
p_wfs0.set_fracsub(0.8)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.7)
p_wfs0.set_gsmag(11)
p_wfs0.set_optthroughput(0.28)
p_wfs0.set_zerop(2.6e10)  # 2.6e10 ph/s/m**2 computed by Rico in R band for MOSAIC
p_wfs0.set_noise(0.2)  # in electrons units
p_wfs0.set_atmos_seen(1)
p_wfs0.set_fstop("square")
p_wfs0.set_fssize(1.6)
rMod = 3
p_wfs0.set_pyr_npts(int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.) * 4))
#p_wfs0.set_pyr_npts(31)
p_wfs0.set_pyr_ampl(rMod)
#p_wfs0.set_pyr_pup_sep(int(2 / 3. * p_wfs0.nxsub)) # diffraction effect
#p_wfs0.set_pyr_pup_sep((p_wfs0.nxsub))
"""
With 52 pixels of margin between 2 edges of pupils on a 240x240 detector and 92 pixels of pupil:
in Compass pupsep is separation between 1 pupil center and Half of detector
pupsep = 52/2+92/2 = 72
"""
p_wfs0.set_pyr_pup_sep((72))

# dm
p_dm0 = ao.Param_dm()
p_dm1 = ao.Param_dm()
p_dms = [p_dm0, p_dm1]
p_dm0.set_type(scons.DmType.PZT)
nact = p_wfs0.nxsub + 1
#nact = 9

#p_dm0.set_nact(nact)
p_dm0.set_nact(73)  #73 actuators for a projected M4 pitch of 53cm
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.1)  # fraction units
# !!!!!!!!!!!!!!!!!!!!!!!!! attention pas autre chose que 0.2 !!!!!!!!!
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(0.01)
p_dm0.set_type_pattern("hexa")
#p_dm0.set_influType("gaussian")
p_dm0.set_influType("radialSchwartz")
"""
p_dm0.set_file_influ_hdf5("/home/fvidal/compass/shesha/data/M4data/elt_influ_spider.h5")
p_dm0.set_center_name("center")
p_dm0.set_cube_name("m_influ")
p_dm0.set_x_name("xpos")
p_dm0.set_y_name("ypos")
p_dm0.set_influ_res("res")
p_dm0.set_diam_dm("diam")
p_dm0.set_diam_dm_proj("diam_projet")
"""

#p_dm0.set_gain(0.2)

p_dm1.set_type("tt")
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(1)
p_dm1.set_push4imat(0.005)
#p_dm1.set_gain(0.2)

# centroiders
p_centroider0 = ao.Param_centroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
p_centroider0.set_type("pyr")
# sp_centroider0.set_thresh(0.)

# p_centroider0.set_type("bpcog")
# p_centroider0.set_nmax(10)

# p_centroider0.set_type("corr")
# p_centroider0.set_type_fct("model")

# controllers
p_controller0 = ao.Param_controller()
p_controllers = [p_controller0]

#p_controller0.set_type("ls")
p_controller0.set_type("generic")

p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(150.)
p_controller0.set_delay(1)
p_controller0.set_gain(1)

# p_controller0.set_modopti(0)
# p_controller0.set_nrec(2048)
# p_controller0.set_nmodes(5064)
# p_controller0.set_gmin(0.001)
# p_controller0.set_gmax(0.5)
# p_controller0.set_ngain(500)

# rtc
#p_rtc = ao.Param_rtc()
#
#p_rtc.set_nwfs(1)
#p_rtc.set_centroiders(p_centroiders)
#p_rtc.set_controllers(p_controllers)
#
