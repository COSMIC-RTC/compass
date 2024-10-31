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
script test to simulate a closed loop for MICADO

Usage:
  micado_loop.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  --bench            For a timed call
  -i --interactive   keep the script interactive
  -d --devices devices Specify the devices
  -n --niter niter   Number of iterations
  -g --generic       Use generic controller
  -f --fast          Compute PSF only during monitoring
"""

from shesha.config import ParamConfig
from docopt import docopt
import numpy as np

if __name__ == "__main__":
    arguments = docopt(__doc__)

    param_file = arguments["<parameters_filename>"]
    compute_tar_psf = not arguments["--fast"]

    config = ParamConfig(param_file)

    # Get parameters from file
    if arguments["--bench"]:
        from shesha.supervisor.benchSupervisor import (
            BenchSupervisor as Supervisor,
        )
    else:
        from shesha.supervisor.compassSupervisor import (
            CompassSupervisor as Supervisor,
        )

    if arguments["--devices"]:
        config.p_loop.set_devices([int(device) for device in arguments["--devices"].split(",")])

    if arguments["--generic"]:
        config.p_controllers[0].set_type("generic")
        print("Using GENERIC controller...")

    if arguments["--niter"]:
        config.p_loop.set_niter(int(arguments["--niter"]))

    supervisor = Supervisor(config)

    # Specific MICADO stuff
    # Index of the actuators in and outside the pupil
    ipos_in, ipos_out = supervisor.basis.compute_ipos_in_pupil(d_obs=11.4, d_pup=37.)
    # Index of the actuators along the spiders
    ipos_spi1, ipos_spi2 = supervisor.basis.compute_ipos_spider()
    # Number of M4 actuators
    Nactu = supervisor.config.p_dms[0]._ntotact
    # ???
    ipos_out2 = np.unique(np.r_[ipos_out, ipos_spi1])
    ipos_in2 = np.arange(Nactu)[np.where(np.isin(np.arange(Nactu), ipos_out2) == False)]
    # Normalization matrix
    IFdelta = supervisor.basis.compute_influ_delta(0)
    # Gendron basis
    Bg, Bgext = supervisor.basis.compute_Bg(ipos_in=ipos_in2, IFdelta=IFdelta)
    # Continuous basis = Modal basis used
    Bc, Bcext = supervisor.basis.compute_Br(Bg, ipos_in2, ipos_out2, IFdelta=IFdelta)
    # Zeros padding to integrate TT
    Bext = np.zeros((Nactu+2,Nactu+2))
    nmodes = Bc.shape[1]-1
    Bext[-2:,0:2] = 0.02 * np.eye(2)
    Bext[:-2, 2:nmodes] = Bc[:, 2:-1].copy()    
    # Push4imat vector
    pushDMMic = 0.01
    pushTTArcsec = 0.005
    modesAmpli = np.ones(Bext.shape[0]) # shape[1] ?
    modesAmpli[0:nmodes - 2] = pushDMMic
    modesAmpli[nmodes - 2:] = pushTTArcsec
    # Compute imat
    imat = supervisor.calibration.do_imat_modal(0, modesAmpli, Bext[:, :nmodes], nmodes_max=0, noise=False, with_turbu=False, push_pull=True)
    # Compute cmat
    rmatc = np.linalg.pinv(imat)
    Nslopes = 26912
    Cext = np.zeros((Nactu + 2, Nslopes))
    Cext[:nmodes, :] = rmatc.copy()
    # CLOSE config
    gain = 0.5 ; qp = 0.05 ; qm = 2 * qp ; trgt = 0.0
    mask = 1.0 * np.r_[np.ones(nmodes), np.zeros(Nactu+2 - nmodes, dtype=np.float32)]
    # Loop and CLOSE setup
    supervisor.rtc.set_modal_integrator_law(0)
    supervisor.modalgains.set_modal_basis(Bext)
    supervisor.modalgains.set_cmat_modal(Cext)
    supervisor.modalgains.set_mask(mask)
    supervisor.modalgains.set_config(0.3, qm, qp, trgt, 1)
    supervisor.modalgains.adapt_modal_gains(True)
    supervisor.rtc.set_gain(0, gain)
    # AO loop
    supervisor.loop(supervisor.config.p_loop.niter, compute_tar_psf=compute_tar_psf)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename

        embed(basename(__file__), locals())
