#!/usr/bin/env python
## @package   shesha.tests
## @brief     Runs a set of tests
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
## @version   5.5.0
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


"""script test to simulate a closed loop

Usage:
  check.py <parameters_filename> [options]

where parameters_filename is the path to the parameters file

Options:
  -h --help                    Show this help message and exit
  -d, --devices devices        Specify the devices
  --displayResult              Just print the results of the check process
  --repportResult=<repport.md> Save the results of the check process into a md_file
"""

from docopt import docopt
import time

if __name__ == "__main__":
    import pandas
    from shesha.supervisor.compassSupervisor import CompassSupervisor
    from shesha.config import ParamConfig

    arguments = docopt(__doc__)

    if arguments["--displayResult"]:
        from os import remove
        from tabulate import tabulate
        from datetime import datetime
        df = pandas.read_hdf("check.h5")
        print(tabulate(df, tablefmt="github", headers="keys"))
        if arguments["--repportResult"]:
            with open(arguments["--repportResult"], 'w') as the_file:
                the_file.write('# E2E Test Report\n')
                the_file.write('\n')
                the_file.write(datetime.now().strftime(
                        '*Report generated on %d-%b-%Y %H:%M:%S by checkCompass.sh*\n'))
                the_file.write('\n')
                the_file.write('[Unit Tests report](report_unit_test.html)\n')
                the_file.write('\n')
                the_file.write('## Summary\n')
                the_file.write('\n')
                the_file.write(str(tabulate(df, tablefmt="github", headers="keys")))
        remove("check.h5")
    else:
        # Get parameters from file
        param_file = arguments["<parameters_filename>"]
        config = ParamConfig(param_file)

        if arguments["--devices"]:
            config.p_loop.set_devices([
                    int(device) for device in arguments["--devices"].split(",")
            ])

        try:
            t0 = time.perf_counter()
            supervisor = CompassSupervisor(config)
            t_init = time.perf_counter() - t0
            is_init = supervisor.is_init
        except BaseException:
            supervisor = None
            is_init = False
            t_init = 0
            SR = "N/A"
            t_loop = 0

        if is_init:    
            try:
                t0 = time.perf_counter()
                supervisor.loop(supervisor.config.p_loop.niter)
                t_loop = time.perf_counter() - t0
                SR = supervisor.target.get_strehl(0)[1]
            except BaseException:
                SR = "N/A"
                t_loop = 0
            
        try:
            df = pandas.read_hdf("check.h5")
        except FileNotFoundError:
            # columns = ["Test name", "Init", "T Init", "SR@100iter", "T Loop"]
            df = pandas.DataFrame()

        idx = len(df.index)
        df.loc[idx, "Test name"] = param_file.split('/')[-1]
        df.loc[idx, "Init"] = str(is_init)
        df.loc[idx, "SR@100iter"] = str(SR)
        df.loc[idx, "T Init"] = str(t_init)
        df.loc[idx, "T Loop"] = str(t_loop/config.p_loop.niter)

        df.to_hdf("check.h5", key="check")
