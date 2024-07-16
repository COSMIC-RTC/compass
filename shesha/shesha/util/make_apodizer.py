## @package   shesha.util.make_apodizer
## @brief     make_apodizer function
## @author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
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


import numpy as np
from scipy.ndimage import interpolation as interp

from . import utilities as util


def make_apodizer(dim, pupd, filename, angle):
    """TODO doc

    Args:

        (int) : im:

        (int) : pupd:

        (str) : filename:

        (float) : angle:
    """

    print("Opening apodizer")
    print("reading file:", filename)
    pup = np.load(filename)
    A = pup.shape[0]

    if (A > dim):
        raise ValueError("Apodizer dimensions must be smaller.")

    if (A != pupd):
        # use misc.imresize (with bilinear)
        print("TODO pup=bilinear(pup,pupd,pupd)")

    if (angle != 0):
        # use ndimage.interpolation.rotate
        print("TODO pup=rotate2(pup,angle)")
        pup = interp.rotate(pup, angle, reshape=False, order=2)

    reg = np.where(util.dist(pupd) > pupd / 2.)
    pup[reg] = 0.

    pupf = np.zeros((dim, dim), dtype=np.float32)

    if (dim != pupd):
        if ((dim - pupd) % 2 != 0):
            pupf[(dim - pupd + 1) / 2:(dim + pupd + 1) / 2, (dim - pupd + 1) /
                 2:(dim + pupd + 1) / 2] = pup

        else:
            pupf[(dim - pupd) / 2:(dim + pupd) / 2, (dim - pupd) / 2:(dim + pupd) /
                 2] = pup

    else:
        pupf = pup

    pupf = np.abs(pupf).astype(np.float32)

    return pupf
