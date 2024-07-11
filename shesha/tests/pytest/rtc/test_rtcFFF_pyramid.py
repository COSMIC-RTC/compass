## @package   shesha.tests
## @brief     Tests the RTC module
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


import numpy as np
import os
from shesha.sutra_wrap import Rtc_FFF as Rtc
from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor
from shesha.config import ParamConfig

precision = 1e-2

config = ParamConfig(os.path.dirname(__file__) + "/../par/test_pyrhr.py")
sup = Supervisor(config)
sup.next()
sup.rtc.open_loop(0)
sup.rtc.close_loop(0)
sup.rtc.do_control(0)
rtc = Rtc()
rtc.add_centroider(sup.context, sup.config.p_wfss[0]._nvalid * sup.config.p_wfss[0].nPupils, 0, sup.config.p_wfss[0].pixsize,
                   False, 0, "maskedpix")
rtc.add_controller(sup.context, "generic", 0, sup.config.p_controllers[0].delay,
                   sup.config.p_controllers[0].nslope, sup.config.p_controllers[0].nactu,
                    idx_centro=np.zeros(1), ncentro=1)
centro = rtc.d_centro[0]
control = rtc.d_control[0]
rtc.d_centro[0].set_npix(sup.config.p_wfss[0].npix)
xvalid = np.array(sup.rtc._rtc.d_centro[0].d_validx)
yvalid = np.array(sup.rtc._rtc.d_centro[0].d_validy)
rtc.d_centro[0].load_validpos(xvalid, yvalid, xvalid.size)
cmat = sup.rtc.get_command_matrix(0)
rtc.d_control[0].set_cmat(cmat)
rtc.d_control[0].set_gain(sup.config.p_controllers[0].gain)
frame = sup.wfs.get_wfs_image(0)
frame /= frame.max()
rtc.d_centro[0].load_img(frame, frame.shape[0])
rtc.d_centro[0].calibrate_img()

rtc.do_centroids(0)
rtc.do_control(0)

dark = np.random.random(frame.shape)
flat = np.random.random(frame.shape)
centro.set_dark(dark, frame.shape[0])
centro.set_flat(flat, frame.shape[0])


def relative_array_error(array1, array2):
    return np.abs((array1 - array2) / array2.max()).max()


def test_doCentroids_maskedPix():
    binimg = np.array(centro.d_img)
    slopes = np.zeros(xvalid.size)
    psum = binimg[xvalid, yvalid].sum() / slopes.size
    for k in range(slopes.size):
        slopes[k] = binimg[xvalid[k], yvalid[k]] / psum - 1 # -1 for ref slopes
    assert (relative_array_error(np.array(control.d_centroids), slopes) <
            precision)
