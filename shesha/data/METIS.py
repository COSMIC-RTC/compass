import shesha.config as conf

# simul_name="METIS_PyWFS_reg"

# loop
p_loop = conf.ParamLoop()
p_loop.set_niter(2000)  # Number of iterations
p_loop.set_ittime(0.001)  # Loop frequency [s] 1000Hz

# geom
p_geom = conf.ParamGeom()
p_geom.set_zenithangle(30.0)  # Zenith angle [deg] 0-30-60 deg
# p_geom.set_pupdiam(1440)

# tel
p_tel = conf.ParamTel()
p_tel.set_diam(37.0)  # Telescope diameter [m] -> METIS: cercle inscrit dans le M1
p_tel.set_cobs(0.3)  # 11.1/37.0 = 0.3 occ [/]
# p_tel.set_type_ap("EELT-Nominal")
p_tel.set_spiders_type("six")  # Nb of spiders [str]
p_tel.set_t_spiders(0.0162162)  # ratio spiders/diam =0.0162162(37m) or 0.62162/39=0.0159 [/]


# atmos
p_atmos = conf.ParamAtmos()
p_atmos.set_r0(0.160204)
# mseeing = 0.98*5.e-7/0.157*206265
# arc2rad = 1./(3600.*180./!pi)
# r0 = 0.5e-6 /(mseeing*arc2rad)
p_atmos.set_nscreens(35)
p_atmos.set_frac(
    [
        0.242,
        0.12,
        0.0968,
        0.059,
        0.0473,
        0.0473,
        0.0473,
        0.0473,
        0.0399,
        0.0324,
        0.0162,
        0.026,
        0.0156,
        0.0104,
        0.0100,
        0.012,
        0.004,
        0.0140,
        0.0130,
        0.007,
        0.016,
        0.0259,
        0.019,
        0.00990,
        0.0062,
        0.0040,
        0.00250,
        0.0022,
        0.00190,
        0.00140,
        0.0011,
        0.0006,
        0.00090,
        0.0005,
        0.0004,
    ]
)
p_atmos.set_alt(
    [
        30.0,
        90.0,
        150.0,
        200.0,
        245.0,
        300.0,
        390.0,
        600.0,
        1130.0,
        1880.0,
        2630.0,
        3500.0,
        4500.0,
        5500.0,
        6500.0,
        7500.0,
        8500.0,
        9500.0,
        10500.0,
        11500.0,
        12500.0,
        13500.0,
        14500.0,
        15500.0,
        16500.0,
        17500.0,
        18500.0,
        19500.0,
        20500.0,
        21500.0,
        22500.0,
        23500.0,
        24500.0,
        25500.0,
        26500.0,
    ]
)
p_atmos.set_windspeed(
    [
        5.5,
        5.5,
        5.1,
        5.5,
        5.6,
        5.7,
        5.8,
        6,
        6.5,
        7,
        7.5,
        8.5,
        9.5,
        11.5,
        17.5,
        23.0,
        26.0,
        29.0,
        32.0,
        27.0,
        22.0,
        14.5,
        9.5,
        6.3,
        5.5,
        6,
        6.5,
        7,
        7.5,
        8,
        8.5,
        9,
        9.5,
        10.0,
        10.0,
    ]
)
p_atmos.set_winddir(
    [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]
)
p_atmos.set_L0(
    [
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
        25.0,
    ]
)  # Outer scale Lo [m]

# target
p_target = conf.ParamTarget()
p_targets = [p_target]

p_target.set_xpos(0.0)
p_target.set_ypos(0.0)
p_target.set_Lambda(2.2)  # Sensing wavelength = K-band
p_target.set_mag(10.0)  # Star mag at the WFS wavelength

# wfs
p_wfs0 = conf.ParamWfs()
p_wfss = [p_wfs0]
p_wfs0.set_type("pyrhr")
p_wfs0.set_nxsub(74)  # Nb sub-apperture over Tel diameter
p_wfs0.set_fracsub(0.5)  # Valid subaperture [%] = fracIlum
p_wfs0.set_npix(4)  # nb of px per subaperture (for simu)
p_wfs0.set_pixsize(0.726053)  # taille d'1 des pixel en arcsec -> ? only for SH-WFS ? ***907566
# pixfov = 0.8 * (2.2 * 1.E-6 * 74. / 37. * 206265.)
#
p_wfs0.set_xpos(0.0)
p_wfs0.set_ypos(0.0)
p_wfs0.set_Lambda(2.2)  # analyse en K-band
p_wfs0.set_gsmag(10.0)  # mag de la GS
p_wfs0.set_gsalt(0.0)  # alt de la GS -> NGS ?
# Field stop:
p_wfs0.set_fssize(1.8)
p_wfs0.set_fstop("round")
# Modulation:
p_wfs0.set_pyr_npts(24)  # sampling [#]
p_wfs0.set_pyr_ampl(4.0)  # amplitude= radius [l/D]
# Noise:
p_wfs0.set_noise(5)  # Noise variance [electron/subaperture]
# Transmissio analyseur:
# =Tatm * Tcfo * Tdet = 0.9*0.25*0.7
p_wfs0.set_atmos_seen(1)
p_wfs0.set_optthroughput(0.1575)
# Zero point at lambda_wfs
p_wfs0.set_zerop(1.66e12)  # [ph/m2/s]

# dm (2 DM, un TT et un HO)
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dms = [p_dm0, p_dm1]

# DM1 - HO
p_dm0.set_type("pzt")
nact = p_wfs0.nxsub + 1
p_dm0.set_nact(nact)
p_dm0.set_alt(0.0)
p_dm0.set_unitpervolt(1.0)
p_dm0.set_push4imat(0.2)
p_dm0.set_coupling(0.2)  # pas toucher
p_dm0.set_thresh(0.3)

# DM2 - TT
p_dm1.set_type("tt")
p_dm1.set_alt(0.0)
p_dm1.set_gain(1.0)  # Gain*0.17
p_dm1.set_unitpervolt(0.0001)
p_dm1.set_push4imat(2000.0)

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]
p_centroider0.set_nwfs(0)
p_centroider0.set_type("pyr")

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]
p_controller0.set_type("ls")
p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(1500.0)  # default YAO = none
p_controller0.set_delay(1)
p_controller0.set_gain(0.17)

# rtc -> pas touche !
# p_rtc=conf.Param_rtc()

# p_rtc.set_nwfs(1)
# p_rtc.set_centroiders(p_centroiders)
# p_rtc.set_controllers(p_controllers)
