# import cProfile
# import pstats as ps

import sys
import os
from shesha_sim.simulator import Simulator

if (len(sys.argv) != 2):
    error = 'command line should be:"python -i test.py parameters_filename"\n with "parameters_filename" the path to the parameters file'
    raise Exception(error)

# get parameters from file
param_file = sys.argv[1]

sim = Simulator(param_file)
sim.init_sim()
sim.loop(sim.config.p_loop.niter)
