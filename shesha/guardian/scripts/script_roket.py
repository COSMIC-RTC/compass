"""script for ROKET

Usage:
  script_roket.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help                   Show this help message and exit
  -s, --savefile savename     Set the name of the ouput h5 file that will be saved in $DATA_GUARDIAN
  -d, --diam diam             Set the telescope diameter [m]
  --niter niter               Set the number of iterations
  --nssp nxsub                Set the number of subapertures of the WFS. Number of actuators is actualized to nxsub+1
  --npix npix                 Set the number of pixels per subap.
  --pixsize pixsize           Set the WFS pixel size [arcsec]
  --nfilt nfilt               Set the number of filtered modes

Usage with Ipython: ipython [-i] script_roket.py -- [options]
"""

from docopt import docopt

import sys
import os
from guardian.roket import Roket

arguments = docopt(__doc__)
param_file = arguments["<parameters_filename>"]
print(arguments)
# Get parameters from file
if arguments["--savefile"]:
    savefile = arguments["--savefile"]
else:
    savefile = "roket_default.h5"

roket = Roket(param_file)

if arguments["--diam"]:
    roket.config.p_tel.set_diam(float(arguments["--diam"]))
if arguments["--niter"]:
    roket.n = int(arguments["--niter"]) + roket.N_preloop
    roket.config.p_loop.set_niter(roket.n)
if arguments["--nssp"]:
    roket.config.p_wfss[0].set_nxsub(int(arguments["--nssp"]))
    roket.config.p_dms[0].set_nact(int(arguments["--nssp"]) + 1)
if arguments["--npix"]:
    roket.config.p_wfss[0].set_npix(int(arguments["--npix"]))
if arguments["--pixsize"]:
    roket.config.p_wfss[0].set_pixsize(float(arguments["--pixsize"]))
if arguments["--nfilt"]:
    roket.config.p_controllers[0].set_maxcond(float(arguments["--nfilt"]))

roket.init_sim()
roket.loop()
roket.save_in_hdf5(savefile)