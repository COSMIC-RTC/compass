# import cProfile
# import pstats as ps

import sys
import os
from shesha_sim.simulator import Simulator, Bench

if not (len(sys.argv) == 2 or (len(sys.argv) == 3 and sys.argv[1] == 'bench')):
    error = 'Command line should be:"python -i test.py parameters_filename"\n'+\
            '    or python -i closed_loop.py bench parameters_filename for a timed call\n' +\
            ' with "parameters_filename" the path to the parameters file'
    raise Exception(error)

# Get parameters from file
if sys.argv[1] == 'bench':
    param_file = sys.argv[2]
    sim = Bench(param_file)
else:
    param_file = sys.argv[1]
    sim = Simulator(param_file)

sim.init_sim()
sim.loop(sim.config.p_loop.niter)
