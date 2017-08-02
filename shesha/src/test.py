import shesha_config as conf
import naga
p_loop = conf.Param_loop()
p_loop.niter = 1000
p_loop.ittime = 1 / 500.

p_geom = conf.Param_geom()
p_geom.zenithangle = 0.
#p_geom.pupdiam = 500

p_tel = conf.Param_tel()
p_tel.set_diam(8.)

p_atmos = conf.Param_atmos()
p_atmos.r0 = 0.16
p_atmos.nscreens = 1
p_atmos.frac = 1.
p_atmos.alt = 0.
p_atmos.windspeed = 10.
p_atmos.winddir = 45.
p_atmos.L0 = 25.

#wfs
p_wfs0= conf.Param_wfs()
p_wfs1= conf.Param_wfs()
p_wfss=[p_wfs0, p_wfs1]

p_wfs0.set_type("sh")
p_wfs0.set_nxsub(16)
p_wfs0.set_npix(8)
p_wfs0.set_pixsize(0.3)
p_wfs0.set_fracsub(0.8)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.5)
p_wfs0.set_gsmag(3.)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(1.e11)
p_wfs0.set_noise(-1.)
p_wfs0.set_atmos_seen(1)

p_wfs1.set_type("sh")
p_wfs1.set_nxsub(16)
p_wfs1.set_npix(8)
p_wfs1.set_pixsize(0.3)
p_wfs1.set_fracsub(0.8)
p_wfs1.set_xpos(0.)
p_wfs1.set_ypos(0.)
p_wfs1.set_Lambda(0.5)
p_wfs1.set_gsmag(3.)
p_wfs1.set_optthroughput(0.5)
p_wfs1.set_zerop(1.e11)
p_wfs1.set_noise(-1.)
p_wfs1.set_atmos_seen(1)
p_wfs1.set_gsalt(90*1.e3)
p_wfs1.set_lltx(0.)
p_wfs1.set_llty(0.)
p_wfs1.set_laserpower(10.)
p_wfs1.set_lgsreturnperwatt(1.e3)
p_wfs1.set_proftype("Exp")
p_wfs1.set_beamsize(0.8)


#p_geom.geom_init(p_tel)

from shesha_init.atmos_init import atmos_init
from shesha_init.wfs_init import wfs_init
from shesha_init.geom_init import tel_init

c = naga.naga_context(0)

Tel = tel_init(c, p_geom, p_tel, p_atmos, p_loop, p_wfss)
Atmos = atmos_init(c, p_atmos, p_tel, p_geom, p_loop)
Wfs = wfs_init(c, p_wfss, p_atmos, p_tel, p_geom, p_loop, Tel)
