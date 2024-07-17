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


import numpy as np
from shesha.sutra_wrap import Sensors, Telescope, Rtc_FFF as Rtc, Atmos


def open_loopSlp(
    tel: Telescope,
    atmos: Atmos,
    wfs: Sensors,
    rtc: Rtc,
    nrec: int,
    ncontrol: int,
    p_wfss: list,
):
    """Return a set of recorded open-loop slopes, usefull for initialize modal control optimization

    Args:

        tel: (Telescope) : Telescope object

        atmos: (Atmos) : Atmos object

        wfs: (Sensors) : Sensors object

        rtc: (Rtc) : Rtc object

        nrec: (int) : number of samples to record

        ncontrol: (int) : controller's index

        p_wfss: (list of ParamWfs) : wfs settings
    """
    # TEST IT
    ol_slopes = np.zeros(
        (sum([2 * p_wfss[i]._nvalid for i in range(len(p_wfss))]), nrec),
        dtype=np.float32,
    )

    print("Recording " + str(nrec) + " open-loop slopes...")
    for i in range(nrec):
        atmos.move_atmos()

        if p_wfss is not None and wfs is not None:
            for j in range(len(p_wfss)):
                wfs.d_wfs[j].d_gs.raytrace(atmos)
                wfs.d_wfs[j].d_gs.raytrace(tel)
                wfs.d_wfs[j].comp_image()
                rtc.d_centro[j].get_cog()
                ol_slopes[
                    j * p_wfss[j]._nvalid * 2 : (j + 1) * p_wfss[j]._nvalid * 2,
                    i,
                ] = np.array(wfs.d_wfs[j].d_slopes)
    print("Done")
    return ol_slopes
