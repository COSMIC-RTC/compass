import sys
sys.path.insert(0,"../chakra/")
import chakra_ao as ao
import chakra as ch
import numpy as np


p_geom=ao.param_geom()
p_atmos=ao.param_atmos()
p_tel=ao.param_tel()
p_target=ao.param_target()

#geom
p_geom.set_zenithangle(10.)
#large test
p_geom.set_pupdiam(1016)
#small test
#p_geom.set_pupdiam(144)

#tel
p_tel.set_diam(4.)
p_tel.set_cobs(0.1)


#atmos
p_atmos.set_r0(0.16)
p_atmos.set_nscreens(1)
p_atmos.set_frac([1.0])
p_atmos.set_alt([0.0])
p_atmos.set_windspeed([20.0])
p_atmos.set_winddir([45])
#large test
p_atmos.set_dim_screens([1024])
p_atmos.set_L0([1.e5])
#small test
"""
p_atmos.set_dim_screens([152])
p_atmos.set_L0([1.e3])
"""

#target
p_target.set_nTargets(1)
p_target.set_xpos([0])
p_target.set_ypos([0.])
p_target.set_Lambda([1.65])
p_target.set_mag([10])



c=ch.chakra_context()

#initialisation:
#   geom
p_geom.geom_init(p_tel,p_geom.pupdiam)
#   atmos
atm=p_atmos.atmos_init(p_tel,p_geom,1./500.)
#   target
tar=p_target.target_init(c,p_atmos,p_geom)

#generate turbulence, raytrace through the atmos
# and display:
#   the turbulence
#   the target's image 
ao.see_atmos_target(50,atm,tar,alt=0,n_tar=0,f=1)
