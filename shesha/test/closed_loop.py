"""script test to simulate a closed loop

Usage:
  closed_loop.py <parameters_filename> [--brama] [--bench]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  --brama            Distribute data with BRAMA
  --bench            For a timed call
"""

from docopt import docopt

import sys
import os
from shesha_sim.simulator import Simulator, SimulatorBrama, Bench

arguments = docopt(__doc__)
param_file = arguments["<parameters_filename>"]

# Get parameters from file
if arguments["--bench"]:
    sim = Bench(param_file)
elif arguments["--brama"]:
    sim = SimulatorBrama(param_file)
else:
    sim = Simulator(param_file)

sim.init_sim()
sim.loop(sim.config.p_loop.niter)
