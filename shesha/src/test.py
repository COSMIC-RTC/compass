import shesha_config as conf
import naga
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

p_loop = conf.Param_loop()
p_loop.niter = 1000
p_loop.ittime = 1 / 500.

p_geom = conf.Param_geom()
p_geom.zenithangle = 0.
#p_geom.pupdiam = 500

p_tel = conf.Param_tel()
p_tel.set_diam(8.)
p_tel.set_cobs(0.12)

p_atmos = conf.Param_atmos()
p_atmos.r0 = 0.16
p_atmos.nscreens = 1
p_atmos.frac = 1.
p_atmos.alt = 0.
p_atmos.windspeed = 10.
p_atmos.winddir = 45.
p_atmos.L0 = 25.

p_target = conf.Param_target()

p_target.set_ntargets(1)
p_target.set_xpos([0])
p_target.set_ypos([0.])
p_target.set_Lambda([1.65])
p_target.set_mag([10])

# wfs
p_wfs0 = conf.Param_wfs()
p_wfs1 = conf.Param_wfs()
p_wfss = [p_wfs0]

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

"""
p_wfs0.set_type("pyrhr")
p_wfs0.set_nxsub(16)
p_wfs0.set_fssize(1.5)
p_wfs0.set_fracsub(0.8)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.5)
p_wfs0.set_gsmag(5.)
p_wfs0.set_optthroughput(0.5)
p_wfs0.set_zerop(1.e11)
p_wfs0.set_noise(-1.)
p_wfs0.set_fstop("round")
p_wfs0.set_pyr_npts(16)
p_wfs0.set_pyr_ampl(3.)
p_wfs0.set_pyr_pup_sep(p_wfs0.nxsub)
p_wfs0.set_atmos_seen(1)
"""
"""
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
"""
p_dm0 = conf.Param_dm()
p_dm1 = conf.Param_dm()
p_dm2 = conf.Param_dm()
p_dms = [p_dm0, p_dm1, p_dm2]
p_dm0.set_type_dm("pzt")
nact = p_wfs0.nxsub + 1
p_dm0.set_nact(nact)
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.3)
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(0.01)
p_dm0.set_push4imat(100.)

p_dm1.set_type_dm("tt")
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(0.0005)
p_dm1.set_push4imat(10.)

p_dm2.set_type_dm("kl")
p_dm2.set_nkl(100)
p_dm2.set_unitpervolt(0.0005)
p_dm2.set_push4imat(10.)


# p_geom.geom_init(p_tel)

from shesha_init.atmos_init import atmos_init
from shesha_init.wfs_init import wfs_init
from shesha_init.geom_init import tel_init
from shesha_init.dm_init import dm_init
from shesha_init.target_init import target_init

c = naga.naga_context(0)

Tel = tel_init(c, p_geom, p_tel, p_atmos, p_loop, p_wfss)
Atmos = atmos_init(c, p_atmos, p_tel, p_geom, p_loop)
Dms = dm_init(c, p_dms, p_wfss, p_geom, p_tel)
Tar = target_init(c, Tel, p_target, p_atmos, p_geom, p_tel, p_dms)
Wfs = wfs_init(c, p_wfss, p_dms, p_atmos, p_tel, p_geom, p_loop, Tel)

Atmos.move_atmos()
Wfs.raytrace(0, b"atmos", Tel, Atmos)
Wfs.comp_img(0)
plt.matshow(Atmos.get_screen(0.))
plt.title("Atmos")
plt.matshow(Wfs.get_phase(0))
plt.title("WFS Phase")
plt.matshow(Wfs.get_binimg(0))
plt.title("WFS img")
