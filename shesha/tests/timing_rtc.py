## @package   shesha.tests
## @brief     Timing of the RTC module
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   5.0.0
## @date      2020/05/18
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
#  All rights reserved.
#  Distributed under GNU - LGPL
#
#  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
#  General Public License as published by the Free Software Foundation, either version 3 of the License,
#  or any later version.
#
#  COMPASS: End-to-end AO simulation tool using GPU acceleration
#  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
#
#  The final product includes a software package for simulating all the critical subcomponents of AO,
#  particularly in the context of the ELT and a real-time core based on several control approaches,
#  with performances consistent with its integration into an instrument. Taking advantage of the specific
#  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
#  conduct large simulation campaigns called to the ELT.
#
#  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
#  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
#  various systems configurations such as multi-conjugate AO.
#
#  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
#  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
#  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.

import numpy as np
import naga as ng
import matplotlib.pyplot as plt
plt.ion()
import os
from shesha.sutra_wrap import Rtc_FFF as Rtc, Rtc_FHF as RtcH
from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor
import carmaWrap
from tqdm import tqdm

dec = 5
sup = Supervisor(
        os.getenv("COMPASS_ROOT") + "/shesha/data/par/par4bench/scao_sh_80x80_8pix.py")
sup.config.p_controllers[0].set_type("generic")
sup.init()
sup.next()
xvalid = np.array(sup.rtc._rtc.d_centro[0].d_validx)
yvalid = np.array(sup.rtc._rtc.d_centro[0].d_validy)
cmat = sup.rtc.get_command_matrix(0)
frame = sup.wfs.get_wfs_image(0)
frame /= frame.max()

rtc = Rtc()
rtc.add_centroider(sup.context, sup.config.p_wfss[0]._nvalid,
                   sup.config.p_wfss[0].npix / 2 - 0.5, sup.config.p_wfss[0].pixsize, False, 0, "cog")
rtc.add_controller(sup.context, sup.config.p_wfss[0]._nvalid, sup.config.p_wfss[0]._nvalid * 2,
                   sup.config.p_controllers[0].nactu, sup.config.p_controllers[0].delay, 0,
                   "generic", idx_centro=np.zeros(1), ncentro=1)
rtc.d_centro[0].set_npix(sup.config.p_wfss[0].npix)
rtc.d_centro[0].load_validpos(xvalid, yvalid, xvalid.size)
rtc.d_control[0].set_cmat(cmat)
rtc.d_control[0].set_gain(sup.config.p_controllers[0].gain)
rtc.d_centro[0].load_img(frame, frame.shape[0])

rtcH = RtcH()
rtcH.add_centroider(sup.context, sup.config.p_wfss[0]._nvalid,
                    sup.config.p_wfss[0].npix / 2 - 0.5, sup.config.p_wfss[0].pixsize, False, 0,
                    "cog")
rtcH.add_controller(sup.context, sup.config.p_wfss[0]._nvalid, sup.config.p_wfss[0]._nvalid * 2,
                    sup.config.p_controllers[0].nactu, sup.config.p_controllers[0].delay, 0,
                    "generic", idx_centro=np.zeros(1), ncentro=1)
rtcH.d_centro[0].set_npix(sup.config.p_wfss[0].npix)
rtcH.d_centro[0].load_validpos(xvalid, yvalid, xvalid.size)
rtcH.d_control[0].set_cmat(cmat)
rtcH.d_control[0].set_gain(sup.config.p_controllers[0].gain)
rtcH.d_centro[0].load_img(frame, frame.shape[0])

timer = carmaWrap.timer()
niter = 100000
FP32 = np.zeros(niter)
FP16 = np.zeros(niter)
FP16TC = np.zeros(niter)

timer.start()
rtc.do_centroids(0)
rtc.do_control(0)
timer.stop()
timer.reset()

for k in tqdm(range(niter)):
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

for k in tqdm(range(niter)):
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

for k in tqdm(range(niter)):
    timer.start()
    rtcH.do_centroids(0)
    rtcH.do_control(0)
    timer.stop()
    FP16TC[k] = timer.total_time
    timer.reset()
