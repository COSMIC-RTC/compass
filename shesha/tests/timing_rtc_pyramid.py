## @package   shesha.tests
## @brief     Timing of te RTC pyramid
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


import os
import numpy as np
import matplotlib.pyplot as plt
from shesha.sutra_wrap import Rtc_FFF as Rtc, Rtc_FHF as RtcH
from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor
import carma
from rich.progress import track

plt.ion()
dec = 5
sup = Supervisor(
        os.getenv("COMPASS_ROOT") + "/shesha/data/par/MICADO/micado_39m_PYR_ELTPupil.py")
sup.config.p_controllers[0].set_type("generic")
sup.config.p_centroiders[0].set_type("maskedpix")
sup.config.p_wfss[0].roket = False
sup.init()
sup.next()
xvalid = np.array(sup.rtc._rtc.d_centro[0].d_validx)
yvalid = np.array(sup.rtc._rtc.d_centro[0].d_validy)
cmat = sup.rtc.get_command_matrix(0)
frame = sup.wfs.get_wfs_image(0)
frame /= frame.max()

rtc = Rtc()
rtc.add_centroider(sup.context, sup.config.p_wfss[0]._nvalid,
                   sup.config.p_wfss[0].npix / 2 - 0.5, sup.config.p_wfss[0].pixsize, False, 0,
                   "maskedpix")
rtc.add_controller(sup.context, "generic", 0, sup.config.p_controllers[0].delay,
                   sup.config.p_controllers[0].nslope, sup.config.p_controllers[0].nactu,
                   idx_centro=np.zeros(1), ncentro=1)
rtc.d_centro[0].set_npix(sup.config.p_wfss[0].npix)
rtc.d_centro[0].load_validpos(xvalid, yvalid, xvalid.size)
rtc.d_control[0].set_cmat(cmat)
rtc.d_control[0].set_gain(sup.config.p_controllers[0].gain)
rtc.d_centro[0].load_img(frame, frame.shape[0])

rtcH = RtcH()
rtcH.add_centroider(sup.context, sup.config.p_wfss[0]._nvalid,
                    sup.config.p_wfss[0].npix / 2 - 0.5, sup.config.p_wfss[0].pixsize, False, 0,
                    "maskedpix")
rtcH.add_controller(sup.context, "generic", 0, sup.config.p_controllers[0].delay,
                    sup.config.p_controllers[0].nslope, sup.config.p_controllers[0].nactu,
                    idx_centro=np.zeros(1), ncentro=1)
rtcH.d_centro[0].set_npix(sup.config.p_wfss[0].npix)
rtcH.d_centro[0].load_validpos(xvalid, yvalid, xvalid.size)
rtcH.d_control[0].set_cmat(cmat)
rtcH.d_control[0].set_gain(sup.config.p_controllers[0].gain)
rtcH.d_centro[0].load_img(frame, frame.shape[0])

timer = carma.timer()
niter = 100000
FP32 = np.zeros(niter)
FP16 = np.zeros(niter)
FP16TC = np.zeros(niter)

timer.start()
rtc.do_centroids(0)
rtc.do_control(0)
timer.stop()
timer.reset()

for k in track(range(niter)):
    timer.start()
    rtc.do_centroids(0)
    rtc.do_control(0)
    timer.stop()
    FP32[k] = timer.total_time
    timer.reset()

timer.start()
rtcH.do_centroids(0)
rtcH.do_control(0)
timer.stop()
timer.reset()

for k in track(range(niter)):
    timer.start()
    rtcH.do_centroids(0)
    rtcH.do_control(0)
    timer.stop()
    FP16[k] = timer.total_time
    timer.reset()

sup._sim.c.activate_tensor_cores(True)
timer.start()
rtcH.do_centroids(0)
rtcH.do_control(0)
timer.stop()
timer.reset()

for k in track(range(niter)):
    timer.start()
    rtcH.do_centroids(0)
    rtcH.do_control(0)
    timer.stop()
    FP16TC[k] = timer.total_time
    timer.reset()
