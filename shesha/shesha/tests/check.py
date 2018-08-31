#!/usr/bin/env python
"""script test to simulate a closed loop

Usage:
  closed_loop.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -d, --devices devices      Specify the devices
  --displayResult    Just print the results of the check process
"""

from docopt import docopt

if __name__ == "__main__":
    import pandas
    from shesha.supervisor.compassSupervisor import CompassSupervisor

    arguments = docopt(__doc__)
    param_file = arguments["<parameters_filename>"]

    if arguments["--displayResult"]:
        import os
        df = pandas.read_hdf("check.h5")
        print(df)
        os.remove("check.h5")
    else:
        # Get parameters from file
        supervisor = CompassSupervisor(param_file)

        if arguments["--devices"]:
            supervisor.config.p_loop.set_devices([
                    int(device) for device in arguments["--devices"].split(",")
            ])
        try:
            supervisor.initConfig()
            isInit = supervisor.isInit()
        except:
            isInit = False
            SR = "N/A"
        try:
            supervisor.loop(supervisor.config.p_loop.niter)
            SR = supervisor.getStrehl(0)[1]
        except:
            SR = "N/A"

        try:
            df = pandas.read_hdf("check.h5")
        except FileNotFoundError:
            columns = ["Test name", "Init", "SR@100"]
            df = pandas.DataFrame(columns=columns)

        idx = len(df.index)
        df.loc[idx, "Test name"] = param_file.split('/')[-1]
        df.loc[idx, "Init"] = isInit
        df.loc[idx, "SR@100"] = SR

        df.to_hdf("check.h5", "check")