## @package   guardians.misc
## @brief     Miscellaneous roket scripts
## @author    Florian Ferreira <florian.ferreira@obspm.fr>
## @version   5.5.0
## @date      2019/01/24
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
import h5py
import numpy as np
import glob


def validfile(filename):
    f = h5py.File(filename)
    if (list(f.attrs.keys()).count("target.Lambda")):
        Lambda = f.attrs["target.Lambda"][0]
    else:
        Lambda = 1.65
    # nactus = f["noise"][:].shape[0]
    niter = f["noise"][:].shape[1]
    P = f["P"][:]
    nmodes = P.shape[0]
    data = np.zeros((nmodes, niter))
    error_list = [
            "noise", "aliasing", "tomography", "filtered modes", "bandwidth",
            "non linearity"
    ]

    for i in error_list:
        data += np.dot(P, f[i][:])

    data = np.var(data, axis=1)
    data = np.sum(data)
    data = np.exp(-data * (2 * np.pi / Lambda)**2)
    data *= np.exp(-f["fitting"].value)

    SR2 = f["SR2"].value
    SR = f["SR"].value

    if (np.abs(data - SR) < 0.05 or np.abs(data - SR2) < 0.05):
        f.attrs["validity"] = True


datapath = "/home/fferreira/Data/correlation/"
filenames = glob.glob(datapath + "roket_8m*.h5")
atm_layer = len(filenames)
ind = 0
for f in filenames:
    validfile(f)
    ind += 1
    print(" reading : %d/%d\r" % (ind, atm_layer), end=' ')
print("read")
