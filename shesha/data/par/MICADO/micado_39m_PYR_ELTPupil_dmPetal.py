import shesha.config as ao
import shesha.constants as scons
import numpy as np

simul_name = ""
layout = "layoutDeArielle"

# loop
p_loop = ao.Param_loop()

p_loop.set_niter(1000)
p_loop.set_ittime(1 / 500.)  # =1/500
p_loop.set_devices([0,1,2,3])
# geom
p_geom = ao.Param_geom()
p_geom.set_zenithangle(0.)

# tel
p_tel = ao.Param_tel()

#E_ELT PUPIL Rico like
p_tel.set_type_ap("EELT")
p_tel.set_diam(40)
p_tel.set_pupangle(0.)  #ELT pup rotation in degrees
# p_tel.set_t_spiders(0.51)  #Spider size in meters
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
p_wfs0 = ao.Param_wfs(roket=True)
#p_wfs0= ao.Param_wfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("pyrhr") # type de WFS: "sh", "pyrhr", "pyr"

p_wfs0.set_nxsub(
        92
)  # 92 sub aps for hexagonal grid of actuators eq. 78 subaps square grid. (pitch = 0.5m)
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
rMod = 3.  # Modulation radius, in lam/D units
nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.) * 4)
p_wfs0.set_pyr_npts(nbPtMod) # nb pts modu around circle
p_wfs0.set_pyr_ampl(rMod)  # define modulation amplitude

clover_pos = np.array([[ 0.00000000e+00,  0.00000000e+00],       [ 4.01503344e-07,  2.53499235e-06],       [ 1.50860454e-06,  4.64300736e-06],       [ 3.05055438e-06,  5.98705007e-06],
       [ 4.64300736e-06,  6.39055139e-06],       [ 5.87299138e-06,  5.87299138e-06],       [ 6.39055139e-06,  4.64300736e-06],       [ 5.98705007e-06,  3.05055438e-06],
       [ 4.64300736e-06,  1.50860454e-06],       [ 2.53499235e-06,  4.01503344e-07],       [ 1.01715049e-21,  6.22825046e-38],       [-2.53499235e-06,  4.01503344e-07],
       [-4.64300736e-06,  1.50860454e-06],       [-5.98705007e-06,  3.05055438e-06],       [-6.39055139e-06,  4.64300736e-06],       [-5.87299138e-06,  5.87299138e-06],
       [-4.64300736e-06,  6.39055139e-06],       [-3.05055438e-06,  5.98705007e-06],       [-1.50860454e-06,  4.64300736e-06],       [-4.01503344e-07,  2.53499235e-06],
       [-2.49130019e-37,  2.03430098e-21],       [-4.01503344e-07, -2.53499235e-06],       [-1.50860454e-06, -4.64300736e-06],       [-3.05055438e-06, -5.98705007e-06],
       [-4.64300736e-06, -6.39055139e-06],       [-5.87299138e-06, -5.87299138e-06],       [-6.39055139e-06, -4.64300736e-06],       [-5.98705007e-06, -3.05055438e-06],
       [-4.64300736e-06, -1.50860454e-06],       [-2.53499235e-06, -4.01503344e-07],       [-3.05145147e-21, -5.60542542e-37],       [ 2.53499235e-06, -4.01503344e-07],
       [ 4.64300736e-06, -1.50860454e-06],       [ 5.98705007e-06, -3.05055438e-06],       [ 6.39055139e-06, -4.64300736e-06],       [ 5.87299138e-06, -5.87299138e-06],
       [ 4.64300736e-06, -6.39055139e-06],       [ 3.05055438e-06, -5.98705007e-06],       [ 1.50860454e-06, -4.64300736e-06],       [ 4.01503344e-07, -2.53499235e-06]])

hypo_pos2 = np.array([[ 8.30566406e-06,  0.00000000e+00],       [ 7.88220488e-06, -7.28701763e-07],       [ 6.69227230e-06, -1.14781454e-06],       [ 4.95995737e-06, -1.01894748e-06],
       [ 3.00501749e-06, -2.30494264e-07],       [ 1.17459828e-06,  1.17459828e-06],       [-2.30494264e-07,  3.00501749e-06],       [-1.01894748e-06,  4.95995737e-06],
       [-1.14781454e-06,  6.69227230e-06],       [-7.28701763e-07,  7.88220488e-06],       [-3.05145147e-22,  8.30566406e-06],       [ 7.28701763e-07,  7.88220488e-06],
       [ 1.14781454e-06,  6.69227230e-06],       [ 1.01894748e-06,  4.95995737e-06],       [ 2.30494264e-07,  3.00501749e-06],       [-1.17459828e-06,  1.17459828e-06],
       [-3.00501749e-06, -2.30494264e-07],       [-4.95995737e-06, -1.01894748e-06],       [-6.69227230e-06, -1.14781454e-06],       [-7.88220488e-06, -7.28701763e-07],
       [-8.30566406e-06, -6.10290295e-22],       [-7.88220488e-06,  7.28701763e-07],       [-6.69227230e-06,  1.14781454e-06],       [-4.95995737e-06,  1.01894748e-06],
       [-3.00501749e-06,  2.30494264e-07],       [-1.17459828e-06, -1.17459828e-06],       [ 2.30494264e-07, -3.00501749e-06],       [ 1.01894748e-06, -4.95995737e-06],
       [ 1.14781454e-06, -6.69227230e-06],       [ 7.28701763e-07, -7.88220488e-06],       [ 9.15435442e-22, -8.30566406e-06],       [-7.28701763e-07, -7.88220488e-06],
       [-1.14781454e-06, -6.69227230e-06],       [-1.01894748e-06, -4.95995737e-06],       [-2.30494264e-07, -3.00501749e-06],       [ 1.17459828e-06, -1.17459828e-06],
       [ 3.00501749e-06,  2.30494264e-07],       [ 4.95995737e-06,  1.01894748e-06],       [ 6.69227230e-06,  1.14781454e-06],       [ 7.88220488e-06,  7.28701763e-07]])

hypo_pos3 = np.array([[ 8.30566406e-06,  0.00000000e+00],       [ 7.80190421e-06, -1.23570023e-06],       [ 6.39055139e-06, -2.07641602e-06],       [ 4.34984649e-06, -2.21635749e-06],
       [ 2.07641602e-06, -1.50860454e-06],       [ 6.14742632e-22,  0.00000000e+00],       [-1.50860454e-06,  2.07641602e-06],       [-2.21635749e-06,  4.34984649e-06],
       [-2.07641602e-06,  6.39055139e-06],       [-1.23570023e-06,  7.80190421e-06],       [-5.08575245e-22,  8.30566406e-06],       [ 1.23570023e-06,  7.80190421e-06],
       [ 2.07641602e-06,  6.39055139e-06],       [ 2.21635749e-06,  4.34984649e-06],       [ 1.50860454e-06,  2.07641602e-06],       [ 6.14742632e-22,  1.22948526e-21],
       [-2.07641602e-06, -1.50860454e-06],       [-4.34984649e-06, -2.21635749e-06],       [-6.39055139e-06, -2.07641602e-06],       [-7.80190421e-06, -1.23570023e-06],
       [-8.30566406e-06, -1.01715049e-21],       [-7.80190421e-06,  1.23570023e-06],       [-6.39055139e-06,  2.07641602e-06],       [-4.34984649e-06,  2.21635749e-06],
       [-2.07641602e-06,  1.50860454e-06],       [-4.30319842e-21,  4.30319842e-21],       [ 1.50860454e-06, -2.07641602e-06],       [ 2.21635749e-06, -4.34984649e-06],
       [ 2.07641602e-06, -6.39055139e-06],       [ 1.23570023e-06, -7.80190421e-06],       [ 1.52572574e-21, -8.30566406e-06],       [-1.23570023e-06, -7.80190421e-06],
       [-2.07641602e-06, -6.39055139e-06],       [-2.21635749e-06, -4.34984649e-06],       [-1.50860454e-06, -2.07641602e-06],       [ 0.00000000e+00,  0.00000000e+00],
       [ 2.07641602e-06,  1.50860454e-06],       [ 4.34984649e-06,  2.21635749e-06],       [ 6.39055139e-06,  2.07641602e-06],       [ 7.80190421e-06,  1.23570023e-06]])

# a = np.zeros((20,2))
# scale = 1541.7678
# a[:,0] = np.concatenate((np.linspace(- 0.02561081/scale, 0.02561081/scale, 10), np.zeros(10)))
# a[:,1] = np.concatenate((np.zeros(10), np.linspace(- 0.02561081/scale, 0.02561081/scale, 10)))
# p_wfs0.set_pyr_pos(a)

# p_wfs0.set_pyr_npts(40)
# p_wfs0.set_pyr_pos(clover_pos*2)

# p_wfs0.set_pyr_npts(40)
# p_wfs0.set_pyr_pos(hypo_pos2*2)

# p_wfs0.set_pyr_npts(1)
# p_wfs0.set_pyr_pos(np.zeros((1,2)))

#p_wfs0.set_pyr_pup_sep(int(2 / 3. * p_wfs0.nxsub)) # diffraction effect
#p_wfs0.set_pyr_pup_sep((p_wfs0.nxsub))

"""
With 52 pixels of margin between 2 edges of pupils on a 240x240 detector and 92 pixels of pupil:
in Compass pupsep is separation between 1 pupil center and Half of detector
pupsep = 52/2+92/2 = 72
"""
p_wfs0.set_pyr_pup_sep(72) # half pupil separation (center-to-center)

# dm
p_dm0 = ao.Param_dm()
p_dm1 = ao.Param_dm()
p_dm2 = ao.Param_dm()
p_dms = [p_dm0, p_dm2, p_dm1]
#p_dms = [p_dm0, p_dm1]
p_dm0.set_type(scons.DmType.PZT)
nact = p_wfs0.nxsub + 1
#nact = 9

#p_dm0.set_nact(nact)
p_dm0.set_nact(75)  # 75 actuators on 40m for a projected M4 pitch of 54.05 cm
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.)  # fraction units
# !!!!!!!!!!!!!!!!!!!!!!!!! attention pas autre chose que 0.2 !!!!!!!!!
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(0.01)
p_dm0.set_type_pattern("hexaM4")
#p_dm0.set_influ_type("gaussian")
p_dm0.set_influ_type("radialSchwartz")
p_dm0.set_margin_out(0.6)
p_dm0.segmented_mirror = True

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
p_centroider0.set_type("pyr")
# p_centroider0.set_type("maskedpix")
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
p_controller0.set_delay(0)   # loop delay. "0 = 1 frame delay".
p_controller0.set_gain(1)
# p_controller0.set_nstates(6)

# p_controller0.set_modopti(0)
# p_controller0.set_nrec(2048)
# p_controller0.set_nmodes(5064)
# p_controller0.set_gmin(0.001)
# p_controller0.set_gmax(0.5)
# p_controller0.set_ngain(500)
