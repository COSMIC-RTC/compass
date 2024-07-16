## @package   shesha.scripts.dm_standalone
## @brief     Python dm standalone script
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @date      2022/01/24
## @copyright 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>
#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>

# COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation, either version 3 of the 
# License, or any later version.

# COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
# If not, see <https://www.gnu.org/licenses/>

# Copyright (C) 2011-2024 COSMIC Team <https//://github.com/COSMIC-RTC/compass>


# import cProfile
# import pstats as ps

# import numpy as np
import carma as ch
import shesha.config as conf
import numpy as np
import matplotlib.pyplot as plt
from shesha.init.geom_init import geom_init_generic
from shesha.init.dm_init import dm_init_standalone

#geom
p_geom = conf.ParamGeom()
geom_init_generic(p_geom, 500)

#dm
p_dm0 = conf.ParamDm()
p_dms = [p_dm0]
p_dm0.set_type("pzt")
# p_dm0.set_pattern("hexa")
p_dm0.set_nact(80)
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.3)
p_dm0.set_coupling(0.2)
p_dm0.set_unitpervolt(0.01)
p_dm0.set_pzt_extent(0)

#   context
c = ch.context.get_instance_1gpu(0)

#   dm
print("->dm")
dms = dm_init_standalone(c, p_dms, p_geom)

print("====================")
print("init done")
print("====================")
print("objects initialzed on GPU:")
print("--------------------------------------------------------")
print(dms)

cmd = np.zeros(5268)
cmd[1111] = 1
dms.set_full_com(cmd)
plt.matshow(dms.d_dms[0].d_shape)
