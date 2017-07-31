import shesha_config as conf

p_loop = conf.Param_loop()
p_loop.niter = 1000
p_loop.ittime = 1 / 500.

p_geom = conf.Param_geom()
p_geom.zenithangle = 0.
p_geom.pupdiam = 500

p_tel = conf.Param_tel()
p_tel.set_diam(8.)

p_atmos = conf.Param_atmos()
p_atmos.r0 = 0.129
p_atmos.nscreens = 1
p_atmos.frac = 1.
p_atmos.alt = 0.
p_atmos.windspeed = 10.
p_atmos.winddir = 45.
p_atmos.L0 = 25.

p_geom.geom_init(p_tel)

from shesha_init.atmos_init import atmos_init

Atmos = atmos_init(p_atmos, p_tel, p_geom, p_loop)
