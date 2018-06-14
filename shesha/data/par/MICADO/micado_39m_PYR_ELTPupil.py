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

#E_ELT PUPIL Rico like
p_tel.set_type_ap("EELT")
p_tel.set_diam(39)
p_tel.set_pupangle(0)  #ELT pup rotation in degrees
p_tel.set_t_spiders(0.51)  #Spider size in meters
"""
#E_ELT PUPIL Alexis like
p_tel.set_diam(38.542)
p_tel.set_cobs(0.28)
p_tel.set_type_ap("EELT-Nominal")
p_tel.set_spiders_type("six")
p_tel.set_pupangle(0)
p_tel.set_t_spiders(0.)
"""
"""
p_tel.set_nbrmissing(7)
p_tel.set_referr(0.01)
p_tel.set_std_tt(0.050) # microns
p_tel.set_std_piston(0.050) # microns
"""

# atmos
p_atmos = ao.Param_atmos()

# p_atmos.set_r0(0.129)
p_atmos.set_r0(0.129)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([10.])
p_atmos.set_winddir([45.])
p_atmos.set_L0([25.])  # Not simulated in Yorick?

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
p_target.set_ntargets(1)
p_target.set_xpos([0])
p_target.set_ypos([0])
p_target.set_Lambda([2.2])
p_target.set_mag([4])
# wfs
p_wfs0 = ao.Param_wfs(error_budget=True)
#p_wfs0= ao.Param_wfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("pyrhr")

p_wfs0.set_nxsub(
        92
)  # 92 sub aps for hexagonal grid of actuators eq. 78 subaps square grid. (pitch = 0.5m)
p_wfs0.set_fracsub(0.1)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.7)
p_wfs0.set_gsmag(11)
p_wfs0.set_optthroughput(0.28)
p_wfs0.set_zerop(2.6e10)  # 2.6e10 ph/s/m**2 computed by Rico in R band for MOSAIC
p_wfs0.set_noise(0.3)  # in electrons units
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
p_wfs0.set_pyr_pup_sep(72)

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
p_dm0.set_thresh(0.2)  # fraction units
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
#p_controller0.set_type("geo")

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
