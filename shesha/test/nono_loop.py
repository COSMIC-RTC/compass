"""script test to simulate a closed loop

Usage:
  closed_loop.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  --brama            Distribute data with BRAMA
  --bench            For a timed call
  -d, --devices devices      Specify the devices
"""

from docopt import docopt

import shesha_sim
from shesha_constants import ControllerType

arguments = docopt(__doc__)
param_file = arguments["<parameters_filename>"]

# Get parameters from file
# if arguments["--bench"]:
#     sim = shesha_sim.Bench(param_file)
# elif arguments["--brama"]:
#     sim = shesha_sim.SimulatorBrama(param_file)
# else:
#     sim = shesha_sim.Simulator(param_file)

sim = shesha_sim.BenchBrama(param_file)
sim.config.p_controller0.set_type(ControllerType.GENERIC)

if arguments["--devices"]:
    devices = []
    for k in range(len(arguments["--devices"])):
        devices.append(int(arguments["--devices"][k]))
    sim.config.p_loop.set_devices(devices)

sim.init_sim()
sim.loop(monitoring_freq=100)
