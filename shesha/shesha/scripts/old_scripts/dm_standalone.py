# import cProfile
# import pstats as ps

import sys
import os
# import numpy as np
import carmaWrap as ch
import shesha.config as conf
import time
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

if (len(sys.argv) == 2):
    # get parameters from file
    param_file = sys.argv[1]
    if (param_file.split('.')[-1] == b"py"):
        filename = param_file.split('/')[-1]
        param_path = param_file.split(filename)[0]
        sys.path.insert(0, param_path)
        exec("import %s as config" % filename.split(".py")[0])
        sys.path.remove(param_path)
    else:
        raise ValueError("Parameter file extension must be .py")

    # init geom
    conf.Param_geom.geom_init(config.p_geom, config.p_tel, 750,
                              config.p_geom.apod)  #apod = apodizer
else:

    class config:
        #geom
        p_geom = conf.Param_geom()
        p_geom.geom_init_generic(500)

        #dm
        p_dm0 = conf.Param_dm()
        p_dms = [p_dm0]
        p_dm0.set_type("pzt")
        # p_dm0.set_pattern("hexa")
        p_dm0.set_nact(80)
        p_dm0.set_alt(0.)
        p_dm0.set_thresh(0.3)
        p_dm0.set_coupling(0.2)
        p_dm0.set_unitpervolt(0.01)


#   context
c = ch.carmaWrap_context(0)
# c = ch.carmaWrap_context(devices=np.array([0,1], dtype=np.int32))
# c.set_activeDevice(0) #useful only if you use ch.carmaWrap_context()
# c = ch.carmaWrap_context(devices=config.p_loop.devices)

config.p_dm0.set_pzt_extent(0)
#   dm
print("->dm")
if config.p_tel:
    dms = conf.dm_init_standalone(config.p_dms, config.p_geom, config.p_tel.diam,
                                  config.p_tel.cobs)
else:
    dms = conf.dm_init_standalone(config.p_dms, config.p_geom)

print("====================")
print("init done")
print("====================")
print("objects initialzed on GPU:")
print("--------------------------------------------------------")
print(dms)
