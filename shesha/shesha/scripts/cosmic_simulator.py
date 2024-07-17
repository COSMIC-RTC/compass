#!/usr/bin/env python

#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
#
# COMPASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COMPASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COMPASS. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2024 COSMIC Team

"""
script test to simulate a closed loop

Usage:
  cosmic_simulator.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -n --niter niter   Number of iterations
  -g --generic       Use generic controller
  -f --fast          Compute PSF only during monitoring
  -d --devices dev   Use devices (comma separated)
"""

from shesha.config import ParamConfig
from docopt import docopt
from shesha.supervisor.cosmicSupervisor import CosmicSupervisor as Supervisor

if __name__ == "__main__":
    arguments = docopt(__doc__)

    param_file = arguments["<parameters_filename>"]
    compute_tar_psf = not arguments["--fast"]

    config = ParamConfig(param_file)

    # Get parameters from file
    if arguments["--devices"]:
        config.p_loop.set_devices([int(device) for device in arguments["--devices"].split(",")])

    if arguments["--generic"]:
        config.p_controllers[0].set_type("generic")
        print("Using GENERIC controller...")

    if arguments["--niter"]:
        config.p_loop.set_niter(int(arguments["--niter"]))

    supervisor = Supervisor(config)

    supervisor.loop(supervisor.config.p_loop.niter, compute_tar_psf=compute_tar_psf)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename

        embed(basename(__file__), locals())
