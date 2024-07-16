## @package   guardians.misc
## @brief     Miscellaneous roket scripts
## @author    Florian Ferreira <florian.ferreira@obspm.fr>
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
import numpy as np
import matplotlib.pyplot as plt
import h5py
from guardians import groot, gamora
import os

filename = os.getenv("DATA_GUARDIAN") + "roket_8m_LE.h5"
Cab = groot.compute_Cerr(filename)
_, _, psfModel, _ = gamora.psf_rec_Vii(filename, fitting=False,
                                       cov=Cab.astype(np.float32))

f = h5py.File(filename, 'r')
tb = f["tomography"][:] + f["bandwidth"][:]

for k in range(10000, 201000, 10000):
    C = tb[:, :k].dot(tb[:, :k].T) / k
    _, _, psfC, _ = gamora.psf_rec_Vii(filename, fitting=False,
                                       covmodes=C.astype(np.float32))
    plt.matshow(
            np.log10(np.abs(psfC - psfModel)), vmin=np.log10(np.abs(psfModel)).min(),
            vmax=np.log10(np.abs(psfModel)).max())
    plt.title("niter = %d" % k)
    plt.colorbar()
