import shesha.config as conf
import numpy as np
simul_name = "mavis"
AOtype="ltao" 

###################
# LOOP
p_loop = conf.Param_loop()
p_loop.set_niter(1)
p_loop.set_ittime(0.0007)  # =1/500

###################
# GEOM
p_geom = conf.Param_geom()
p_geom.set_zenithangle(0.)

###################
# TEL
p_tel = conf.Param_tel()
p_tel.set_diam(8.0)
p_tel.set_cobs(0.12)

###################
# ATMOS
p_atmos = conf.Param_atmos()
p_atmos.set_r0(0.12)
p_atmos.set_nscreens(4)
p_atmos.set_frac([0.5, 0.2, 0.2, 0.1])
p_atmos.set_alt([0., 1000.0, 4500., 9000.])
p_atmos.set_windspeed([10.0, 10.0, 10.0, 10.0])
p_atmos.set_winddir([0., 10., 20., 25.])
p_atmos.set_L0([25., 25., 25., 25.])

###################
# TARGET
p_targets=[]
NTAR= 1

tar_xpos=np.array([ 0.])
tar_ypos=np.array([ 0.])
#NTAR=49
#xpos, ypos = np.meshgrid(np.linspace(-15, 15, 7), np.linspace(-15, 15, 7))
for t in np.arange(NTAR):
    p_targets.append(conf.Param_target())
    p_targets[t].set_xpos(tar_xpos[t])
    p_targets[t].set_ypos(tar_ypos[t])
    #IR 1.6 / visible 0.65
    #p_targets[t].set_Lambda( 0.65)
    p_targets[t].set_Lambda( 1.6)
    p_targets[t].set_mag(10.)
    #1 DM per target: pzt
    p_targets[t].set_dms_seen(np.array([0]))
    #2 DM per target: pzt + tt
    #p_targets[t].set_dms_seen(np.array([t,NTAR]))


###################
# WFS
RADIUS      = 15
FRACSUB     = 0.55
NXSUB_LGS   = 40
NXSUB_NGS   = 40
NXSUB_TAR   = max(NXSUB_LGS,NXSUB_NGS)
NLGS        = 6
NNGS        = 1
NTS         = 0

p_wfs_lgs = []
p_wfs_ngs = []
p_wfs_ts  = []

#lgs position
x = np.linspace(0, 2 * np.pi, NLGS+1)[:-1]
lgs_xpos = RADIUS * np.cos(x)
lgs_ypos = RADIUS * np.sin(x)

#NGS asterism
# closest star from asterism F2 to axis
#asterism_x = [ -1.213]
#asterism_y = [ 24.719 ]
asterism_x = [ -1.213]
asterism_y = [ 8 ]

#Truth Sensors position
TS_xpos=[]
TS_ypos=[]

# add 1 position for target
asterism_x = asterism_x + [0.]
asterism_y = asterism_y + [0.]

#create wfs lists
#LGS
for i in range(NLGS):
    p_wfs_lgs.append(conf.Param_wfs())
#NGS
for i in range(NNGS+1):
    p_wfs_ngs.append(conf.Param_wfs())
#TS
for i in range(NTS):
    p_wfs_ts.append(conf.Param_wfs())
#concatenate LGS and NGS
p_wfss = p_wfs_lgs + p_wfs_ngs

#setting LGS
for p_wfs in p_wfs_lgs:
    k=p_wfs_lgs.index(p_wfs)
    p_wfs.set_type("sh")
    p_wfs.set_nxsub(NXSUB_LGS)
    p_wfs.set_npix(8)
    p_wfs.set_pixsize(0.5)
    p_wfs.set_fracsub(FRACSUB)
    p_wfs.set_xpos(lgs_xpos[k])
    p_wfs.set_ypos(lgs_ypos[k])

    p_wfs.set_Lambda(0.589)
    p_wfs.set_gsmag(8.)
    p_wfs.set_optthroughput(0.5)
    p_wfs.set_zerop(1.e11)
    p_wfs.set_noise(0.3)
    p_wfs.set_atmos_seen(1)

    p_wfs.set_gsalt(90 * 1.e3)
    p_wfs.set_lltx(0.)
    p_wfs.set_llty(0.)
    p_wfs.set_laserpower(10)
    p_wfs.set_lgsreturnperwatt(1.e3)
    p_wfs.set_proftype("Gauss1")
    p_wfs.set_beamsize(0.8)

#setting NGS / target
for p_wfs in p_wfs_ngs:
    k=p_wfs_ngs.index(p_wfs)
    if( k >= NNGS):
        nxsub=NXSUB_TAR
        p_wfs.set_fracsub(FRACSUB)
    else:
        nxsub=NXSUB_NGS
        p_wfs.set_fracsub(FRACSUB)
        #p_wfs.set_fracsub(0.)
        #p_wfs.set_is_low_order(True)

    p_wfs.set_type("sh")
    p_wfs.set_nxsub(nxsub)
    p_wfs.set_npix(8)
    p_wfs.set_pixsize(0.5)

    p_wfs.set_Lambda(0.589)
    p_wfs.set_gsmag(8.)
    p_wfs.set_optthroughput(0.5)
    p_wfs.set_zerop(1.e11)
    p_wfs.set_noise(0.3)
    p_wfs.set_atmos_seen(1)
    p_wfs.set_xpos(asterism_x[k])
    p_wfs.set_ypos(asterism_y[k])

#setting TS
for p_wfs in p_wfs_ts:
    k=p_wfs_ts.index(p_wfs)
    p_wfs.set_xpos(TS_xpos[k])
    p_wfs.set_ypos(TS_ypos[k])

###################
# DM
p_dms = [ ]
#adding target DMs
p_dms.append(conf.Param_dm())
p_dms[-1].set_type("pzt")
p_dms[-1].set_nact(41)
p_dms[-1].set_alt(0.)
p_dms[-1].set_thresh(0.3) # change actuator threshold
p_dms[-1].set_coupling(0.3)
p_dms[-1].set_unitpervolt(1)
p_dms[-1].set_push4imat(1)

##tip-tilt
#p_dms.append(conf.Param_dm())
#p_dms[-1].set_type("tt")
#p_dms[-1].set_alt(0.)
#p_dms[-1].set_unitpervolt(0.0005)
#p_dms[-1].set_push4imat(10.)

###################
# CENTROIDERS
p_centroiders = []
for i in range( NNGS+NLGS):
    p_centroiders.append(conf.Param_centroider())

for p_centroider in p_centroiders:
    k=p_centroiders.index(p_centroider)
    p_centroider.set_nwfs(k)
    p_centroider.set_type("cog")

###################
# CONTROLLERS
p_controller0 = conf.Param_controller()
p_controllers = [p_controller0]

p_controller0.set_type("generic")
p_controller0.set_nwfs(np.arange(NLGS+NNGS))
p_controller0.set_ndm(list(range(len(p_dms))))
p_controller0.set_maxcond(150.)
p_controller0.set_delay(1.)
p_controller0.set_gain(0.3)
