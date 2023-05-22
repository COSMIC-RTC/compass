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
# p_geom.set_pupdiam(1104)

# tel
p_tel = ao.Param_tel()

#E_ELT PUPIL Rico like
p_tel.set_type_ap("EELT")
p_tel.set_diam(40)
p_tel.set_pupangle(0.)  #ELT pup rotation in degrees
# p_tel.set_t_spiders(0.51)  #Spider size in meters
p_tel.set_t_spiders(0.51)  #Spider size in meters
#p_tel.set_std_tt(0.03)
#p_tel.set_std_piston(0.02)

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

# -----------------------------
# 35 layers profile Q3
# -----------------------------

r0Q3 = 0.1275

zenithAngle = 30.
altESO = np.array([
        30, 90, 150, 200, 245, 300, 390, 600, 1130, 1880, 2630, 3500, 4500, 5500, 6500,
        7500, 8500, 9500, 10500, 11500, 12500, 13500, 14500, 15500, 16500, 17500, 18500,
        19500, 20500, 21500, 22500, 23500, 24500, 25500, 26500
]) / np.cos(zenithAngle * 2 * np.pi / 360)
altESO = altESO.astype(int)

fracQ3 = [
        25.5, 11.9, 9.32, 5.57, 4.5, 4.5, 4.5, 4.5, 4.19, 4.04, 2.02, 3.04, 1.82, 1.21,
        0.86, 1.03, 0.34, 1.2, 1.11, 0.6, 1.43, 2.31, 1.7, 0.88, 0.55, 0.36, 0.22, 0.19,
        0.17, 0.12, 0.1, 0.06, 0.08, 0.04, 0.04
]

windESO = [
        5.5, 5.5, 5.1, 5.5, 5.6, 5.7, 5.8, 6, 6.5, 7, 7.5, 8.5, 9.5, 11.5, 17.5, 23, 26,
        29, 32, 27, 22, 14.5, 9.5, 6.3, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10
]

r0 = r0Q3
alt = altESO
frac = fracQ3
wind = windESO

nbLayers = len(alt)

# atmos
p_atmos = ao.Param_atmos()

## 1 Layer

p_atmos.set_r0(0.129)
# p_atmos.set_r0(0.078)
# p_atmos.set_r0(0.215)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([9.1])
p_atmos.set_winddir([45.])
p_atmos.set_L0([25.])  # Not simulated in Yorick?

## 35 Layers

# p_atmos.set_r0(r0)
# p_atmos.set_nscreens(nbLayers)
# p_atmos.set_frac(frac)
# p_atmos.set_alt(alt)
# p_atmos.set_windspeed(wind)
# p_atmos.set_winddir([45.] * nbLayers)
# p_atmos.set_L0([25.] * nbLayers)  # Not simulated in Yorick?

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

# wfs
p_wfs0 = ao.Param_wfs()
#p_wfs0= ao.Param_wfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("pyrhr") # type de WFS: "sh", "pyrhr", "pyr"
# FDR design has been done with 96 subaps in 38.542 m. This leads to subaps
# of 0.40148 metres. 
p_wfs0.set_nxsub(100)
p_wfs0.set_fracsub(0.01) # 0.01 for pixels below spider, nominal = 0.1 Minimal illumination fraction
p_wfs0.set_xpos(0.)     # direction of guide star in X (arcsec)
p_wfs0.set_ypos(0.)     # direction of guide star in Y (arcsec)
p_wfs0.set_Lambda(0.7)  # wavelength (microns)
p_wfs0.set_gsmag(11)    # magnitude of guide star
p_wfs0.set_optthroughput(0.28) # optical transmission
p_wfs0.set_zerop(2.6e10)  # 2.6e10 ph/s/m**2 computed by Rico in R band for MOSAIC
p_wfs0.set_noise(0.3)   # units: electrons/pixel
p_wfs0.set_atmos_seen(1)   # tell if atmos is seen or not
p_wfs0.set_fstop("square") # shape of field stop, "round", "square"
p_wfs0.set_fssize(1.6)     # size of field stop (arcsec)
# p_wfs0.set_fssize(0.0185)
rMod = 3.  # Modulation radius, in lam/D units
nbPtMod = int(np.ceil(int(rMod * 2 * np.pi) / 4.) * 4)
p_wfs0.set_pyr_npts(nbPtMod) # nb pts modu around circle
p_wfs0.set_pyr_ampl(rMod)  # define modulation amplitude
#p_wfs0.set_nPupils(6)    #8, 24

# p_wfs0.set_pyr_compute_focalplane(True)

#p_wfs0.set_pyr_npts(1)
#p_wfs0.set_pyr_pos(np.zeros((1,2)))

#p_wfs0.set_pyr_pup_sep(int(2 / 3. * p_wfs0.nxsub)) # diffraction effect
#p_wfs0.set_pyr_pup_sep((p_wfs0.nxsub))

npts = 36    # number of modulation points
t = np.linspace(0, 2 * np.pi * (1-1/npts), npts)
#p_wfs0.set_pyr_npts(npts)

nb = 4    # number of branches / edges of the clover
k = 2
a_mod = 3    # modulation amplitude [lambda / D]

path = np.zeros((npts, 2))
n_samp = 4096    # * 2
r_fac = a_mod / (nb-1 + k) * 2 * 0.7e-6 * 3600 * 180 / n_samp / 40
temp1 = ((nb - 1) * np.cos(t) + k * np.cos((nb - 1) * t))
temp2 = ((nb - 1) * np.sin(t) - k * np.sin((nb - 1) * t))

# path[:, 0] = r_fac * temp1
# path[:, 1] = r_fac * temp2
theta = 0
path[:, 0] = r_fac * (temp1 * np.cos(theta) - temp2 * np.sin(theta))
path[:, 1] = r_fac * (temp1 * np.sin(theta) + temp2 * np.cos(theta))
#p_wfs0.set_pyr_pos(path)

"""
With 52 pixels of margin between 2 edges of pupils on a 240x240 detector and 92 pixels of pupil:
in Compass pupsep is separation between 1 pupil center and Half of detector
pupsep = 52/2+92/2 = 72
"""
p_wfs0.set_pyr_pup_sep(72)    # 72, 72//8 # half pupil separation (center-to-center)
# p_wfs0.set_pyr_pup_sep(2)
# dm
p_dm0 = ao.Param_dm()
p_dm1 = ao.Param_dm()
p_dm2 = ao.Param_dm()
p_dms = [p_dm0, p_dm2, p_dm1]
#p_dms = [p_dm0, p_dm1]
p_dm0.set_type(scons.DmType.PZT)
nact = p_wfs0.nxsub + 1
#nact = 9

"""
#p_dm0.set_nact(nact)
p_dm0.set_nact(75)  # 75 actuators on 40m for a projected M4 pitch of 54.05 cm
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.)  # fraction units
# !!!!!!!!!!!!!!!!!!!!!!!!! attention pas autre chose que 0.2 !!!!!!!!!
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(0.01)
p_dm0.set_type_pattern("hexaM4")
#p_dm0.set_influType("gaussian")
p_dm0.set_influType("radialSchwartz")
p_dm0.set_margin_out(0.6)
p_dm0.segmented_mirror = True
"""

p_dm0.set_unitpervolt(1)
p_dm0.set_thresh(0.)  # fraction units
p_dm0.set_margin_out(0.6)
#p_dm0.set_file_influ_fits('/home/abertrou/m4_eelt_compass/testJojo.fits')
p_dm0.set_file_influ_fits('cropped_M4IF.fits')
p_dm0.set_push4imat(0.01)
p_dm0.set_diam_dm(40.0)

#p_dm0.set_gain(0.2)

p_dm1.set_type("tt")
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(1)
p_dm1.set_push4imat(0.005)
#p_dm1.set_gain(0.2)

p_dm2.set_type('pzt')
p_dm2.set_alt(0.)
p_dm2.set_unitpervolt(1)
p_dm2.set_push4imat(0.01)
p_dm2.set_influ_type("petal")

# centroiders
p_centroider0 = ao.Param_centroider()

p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
# p_centroider0.set_method(2)
# p_centroider0.set_type("pyr")
p_centroider0.set_type("maskedpix")
# sp_centroider0.set_thresh(0.)

# p_centroider0.set_type("bpcog")
# p_centroider0.set_nmax(10)

# p_centroider0.set_type("corr")
# p_centroider0.set_type_fct("model")

# controllers
p_controller0 = ao.Param_controller()
p_controllers = [p_controller0]

#p_controller0.set_type("ls")     # V(k) = V(k-1) + g.R.m(k)
p_controller0.set_type("generic") # V(k) = a.E.V(k-1) + g.R.m(k)
#p_controller0.set_type("geo")    # bypass the WFS (direct DM proj)

p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 2])
p_controller0.set_maxcond(150.)
p_controller0.set_delay(1)   # loop delay. "0 = 1 frame delay".
p_controller0.set_gain(1)

"""
 PERFECT CORONO FOR MICADO

p_coronos = [ao.Param_corono()]


p_coronos[0].set_type("perfect")  # coronagraph type : "perfect", "SPHERE_APLC", "custom"
p_coronos[0].set_wavelength_0(1.667)  # coronagraph wavelength in micron
p_coronos[0].set_image_sampling(2)  # number of pixel in lambda/D (2 : Shannon sampling)
p_coronos[0].set_dim_image(250)     # size of the science image in pixel
p_coronos[0].set_delta_wav(0.054)
p_coronos[0].set_nb_wav(3)
"""




"""
CUSTOM CORONO FOR MICADO

p_coronos = [ao.Param_corono()]


p_coronos[0].set_type("custom")
p_coronos[0].set_wavelength_0(1.667)
p_coronos[0].set_apodizer_name('SPHERE_APLC_apodizer_APO1')  # 'SPHERE_APLC_apodizer_APO1' or path to fits file of size (pupdiam, pupdiam)
p_coronos[0].set_focal_plane_mask_name('classical_Lyot')  # 'classical_Lyot', "SPHERE_APLC_fpm_ALC1" or 2 or 3
                                                      # or path to fits file of user defined size [dimx, dimy, nb_wav]
                                                      # or [dimx, dimy] if same for all wavelength
p_coronos[0].set_lyot_fpm_radius(2.15)  # radius of the mask in lambda/D units (required)
# p_corono.set_fpm_sampling()  # lambda/D size in pixel (default = 20.)
                             # optional if classical lyot
p_coronos[0].set_lyot_stop_name("/home/fvidal/spupELT.fits")  # 'SPHERE_APLC_Lyot_stop', or path to fits file of size (pupdiam, pupdiam)
                                                      # homogeneous to get_s_pupil() size in COMPASS
p_coronos[0].set_dim_image(400)
p_coronos[0].set_image_sampling(3)

# p_corono.set_babinet_trick() # True or False for enabling Babinet propagation method



"""
# p_controller0.set_nstates(6)

# p_controller0.set_modopti(0)
# p_controller0.set_nrec(2048)
# p_controller0.set_nmodes(5064)
# p_controller0.set_gmin(0.001)
# p_controller0.set_gmax(0.5)
# p_controller0.set_ngain(500)
