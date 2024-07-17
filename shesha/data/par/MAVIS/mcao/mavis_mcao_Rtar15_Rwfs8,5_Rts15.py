import shesha.config as conf
import numpy as np

simul_name = "mavis"
AOtype = "mcao"

NFILT = 137


# LOOP
p_loop = conf.ParamLoop()
p_loop.set_niter(1)
p_loop.set_ittime(0.0007)  # =1/1500


# GEOM
p_geom = conf.ParamGeom()
p_geom.set_zenithangle(0.0)


# TEL
p_tel = conf.ParamTel()
p_tel.set_diam(8.0)
p_tel.set_cobs(0.12)


# ATMOS
p_atmos = conf.ParamAtmos()
p_atmos.set_r0(0.12)
p_atmos.set_nscreens(4)
p_atmos.set_frac([0.5, 0.2, 0.2, 0.1])
p_atmos.set_alt([0.0, 4490.0, 4500.0, 9000.0])
p_atmos.set_windspeed([10.0, 10.0, 10.0, 10.0])
p_atmos.set_winddir([0.0, 10.0, 20.0, 25.0])
p_atmos.set_L0([25.0, 25.0, 25.0, 25.0])


# TARGET
p_targets = []
RADIUS_TAR = 15
NTAR_side = 10
NTAR = NTAR_side * NTAR_side
tar_xpos = np.array([0.0])
tar_ypos = np.array([0.0])
# NTAR=49
if NTAR > 1:
    tar_xpos, tar_ypos = np.meshgrid(
        np.linspace(-RADIUS_TAR, RADIUS_TAR, NTAR_side),
        np.linspace(-RADIUS_TAR, RADIUS_TAR, NTAR_side),
    )
for t in np.arange(NTAR):
    p_targets.append(conf.ParamTarget())
    p_targets[t].set_xpos(tar_xpos.flatten()[t])
    p_targets[t].set_ypos(tar_ypos.flatten()[t])
    # IR 1.6 / visible 0.65
    p_targets[t].set_Lambda(0.65)
    # p_targets[t].set_Lambda(1.6)
    p_targets[t].set_mag(10.0)
    # 1 DM per target: pzt
    p_targets[t].set_dms_seen([0, 1, 2, 3])
    # 2 DM per target: pzt + tt
    # p_targets[t].set_dms_seen(np.array([t,NTAR]))


# WFS
RADIUS = 8.5
FRACSUB = 0.8
NXSUB_LGS = 40
NXSUB_NGS = 40
NXSUB_TAR = max(NXSUB_LGS, NXSUB_NGS)
NLGS = 6
NNGS = 3
NTS_side = 5
NTS = NTS_side * NTS_side

p_wfs_lgs = []
p_wfs_ngs = []
p_wfs_ts = []

# lgs position
x = np.linspace(0, 2 * np.pi, NLGS + 1)[:-1]
lgs_xpos = RADIUS * np.cos(x)
lgs_ypos = RADIUS * np.sin(x)

# NGS asterism
# closest star from asterism F2 to axis
# asterism_x = [ -1.213]
# asterism_y = [ 24.719 ]
x = np.linspace(np.pi / 6.0, 2 * np.pi + np.pi / 6.0, NNGS + 1)[:-1]
asterism_x = RADIUS * np.cos(x)
asterism_y = RADIUS * np.sin(x)

# Truth Sensors position
radius = 15
lspace = np.linspace(-radius + 1, radius - 1, NTS_side)
mesh = np.meshgrid(lspace, lspace)
TS_xpos = mesh[0].flatten()  # radius * np.cos(x)
TS_ypos = mesh[1].flatten()  # radius * np.sin(x)

# add 1 position for target
asterism_x = np.append(asterism_x, [0.0])
asterism_y = np.append(asterism_y, [0.0])

# create wfs lists
# LGS
for i in range(NLGS):
    p_wfs_lgs.append(conf.ParamWfs())
# NGS
for i in range(NNGS + 1):
    p_wfs_ngs.append(conf.ParamWfs())
# TS
for i in range(NTS):
    p_wfs_ts.append(conf.ParamWfs())
# concatenate LGS and NGS
p_wfss = p_wfs_lgs + p_wfs_ngs

# setting LGS
for p_wfs in p_wfs_lgs:
    k = p_wfs_lgs.index(p_wfs)
    p_wfs.set_type("sh")
    p_wfs.set_nxsub(NXSUB_LGS)
    p_wfs.set_npix(8)
    p_wfs.set_pixsize(0.5)
    p_wfs.set_fracsub(FRACSUB)
    p_wfs.set_xpos(lgs_xpos[k])
    p_wfs.set_ypos(lgs_ypos[k])

    p_wfs.set_Lambda(0.589)
    p_wfs.set_gsmag(8.0)
    p_wfs.set_optthroughput(0.5)
    p_wfs.set_zerop(1.0e11)
    p_wfs.set_noise(0.3)
    p_wfs.set_atmos_seen(1)

    p_wfs.set_gsalt(90 * 1.0e3)
    p_wfs.set_lltx(0.0)
    p_wfs.set_llty(0.0)
    p_wfs.set_laserpower(10)
    p_wfs.set_lgsreturnperwatt(1.0e3)
    p_wfs.set_proftype("Gauss1")
    p_wfs.set_beamsize(0.8)

# setting NGS / target
for p_wfs in p_wfs_ngs:
    k = p_wfs_ngs.index(p_wfs)
    if k >= NNGS:
        nxsub = NXSUB_TAR
        p_wfs.set_fracsub(FRACSUB)
    else:
        nxsub = NXSUB_NGS
        p_wfs.set_fracsub(FRACSUB)
        # p_wfs.set_fracsub(0.)
        # p_wfs.set_is_low_order(True)

    p_wfs.set_type("sh")
    p_wfs.set_nxsub(nxsub)
    p_wfs.set_npix(8)
    p_wfs.set_pixsize(0.5)

    p_wfs.set_Lambda(0.589)
    p_wfs.set_gsmag(8.0)
    p_wfs.set_optthroughput(0.5)
    p_wfs.set_zerop(1.0e11)
    p_wfs.set_noise(0.3)
    p_wfs.set_atmos_seen(1)
    p_wfs.set_xpos(asterism_x[k])
    p_wfs.set_ypos(asterism_y[k])

# setting TS
for p_wfs in p_wfs_ts:
    k = p_wfs_ts.index(p_wfs)
    p_wfs.set_xpos(TS_xpos[k])
    p_wfs.set_ypos(TS_ypos[k])


# DM
p_dm0 = conf.ParamDm()
p_dm1 = conf.ParamDm()
p_dm2 = conf.ParamDm()
p_dm3 = conf.ParamDm()
p_dms = [p_dm0, p_dm1, p_dm2, p_dm3]
# adding target DMs
p_dm0.set_type("pzt")
p_dm0.set_nact(41)
p_dm0.set_alt(0.0)
p_dm0.set_thresh(0.3)  # change actuator threshold
p_dm0.set_coupling(0.3)
p_dm0.set_unitpervolt(1)
p_dm0.set_push4imat(1)

p_dm1.set_type("pzt")
p_dm1.set_nact(41)
p_dm1.set_alt(4500.0)
p_dm1.set_thresh(0.3)
p_dm1.set_coupling(0.3)
p_dm1.set_unitpervolt(0.1)
p_dm1.set_push4imat(1)

p_dm2.set_type("pzt")
p_dm2.set_nact(61)
p_dm2.set_alt(9000.0)
p_dm2.set_thresh(0.3)
p_dm2.set_coupling(0.3)
p_dm2.set_unitpervolt(1)
p_dm2.set_push4imat(1)

p_dm3.set_type("tt")
p_dm3.set_alt(0.0)
p_dm3.set_unitpervolt(0.0005)
p_dm3.set_push4imat(10.0)


# CENTROIDERS
p_centroiders = []
for i in range(NNGS + NLGS):
    p_centroiders.append(conf.ParamCentroider())

for p_centroider in p_centroiders:
    k = p_centroiders.index(p_centroider)
    p_centroider.set_nwfs(k)
    p_centroider.set_type("cog")
    if p_wfss[k].get_gsalt() > 0:
        p_centroider.set_filter_TT(True)


# CONTROLLERS
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("generic")
p_controller0.set_nwfs(np.arange(NLGS + NNGS))
p_controller0.set_ndm(list(range(len(p_dms))))
p_controller0.set_maxcond(150.0)
p_controller0.set_delay(1.0)
p_controller0.set_gain(0.3)
