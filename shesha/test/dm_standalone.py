
# import cProfile
# import pstats as ps

import sys
import os
# import numpy as np
import naga as ch
import shesha as ao
import time
import matplotlib.pyplot as plt
import hdf5_utils as h5u
plt.ion()

if(len(sys.argv) != 2):
    error = 'command line should be:"python -i test.py parameters_filename"\n with "parameters_filename" the path to the parameters file'
    raise StandardError(error)

# get parameters from file
param_file = sys.argv[1]
if(param_file.split('.')[-1] == "py"):
    filename = param_file.split('/')[-1]
    param_path = param_file.split(filename)[0]
    sys.path.insert(0, param_path)
    exec("import %s as config" % filename.split(".py")[0])
    sys.path.remove(param_path)
else:
    raise ValueError("Parameter file extension must be .py")

# initialisation:

#   context
c = ch.naga_context(0)
# c = ch.naga_context(devices=np.array([0,1], dtype=np.int32))
# c.set_activeDevice(0) #useful only if you use ch.naga_context()
# c = ch.naga_context(devices=config.p_loop.devices)

config.p_dm0.set_pzt_extent(3)
config.p_geom.geom_init(config.p_tel, 500)

#   dm
print "->dm"
dms = ao.dm_init_2(config.p_dms, config.p_geom)

print "===================="
print "init done"
print "===================="
print "objects initialzed on GPU:"
print "--------------------------------------------------------"
print dms
