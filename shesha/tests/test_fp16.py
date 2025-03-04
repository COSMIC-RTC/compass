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
import os
from shesha.sutra_wrap import Rtc_FHF as Rtc
from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor
from rich.progress import track

precision = 1e-2
sup = Supervisor(os.getenv("COMPASS_ROOT") + "/shesha/data/par/par4bench/scao_sh_16x16_8pix.py")
sup.config.p_controllers[0].delay = 0
sup.init()
# sup.next()
sup.rtc.open_loop(0)
sup.rtc.close_loop(0)
sup2 = Supervisor(os.getenv("COMPASS_ROOT") + "/shesha/data/par/par4bench/scao_sh_16x16_8pix.py")
sup2.config.p_controllers[0].delay = 0
sup2.init_config()
# sup2.single_next()
sup2.open_loop()
sup2.close_loop()
rtc = Rtc()
rtc.add_centroider(
    sup.context,
    sup.config.p_wfss[0]._nvalid,
    sup.config.p_wfss[0].npix / 2 - 0.5,
    sup.config.p_wfss[0].pixsize,
    False,
    0,
    "cog",
)
rtc.add_controller(
    sup.context,
    "generic",
    0,
    sup.config.p_controllers[0].delay,
    sup.config.p_wfss[0]._nvalid * 2,
    sup.config.p_controllers[0].nactu,
    idx_centro=np.zeros(1),
    ncentro=1,
)
centro = rtc.d_centro[0]
control = rtc.d_control[0]
rtc.d_centro[0].set_npix(sup.config.p_wfss[0].npix)
xvalid = np.array(sup.rtc._rtc.d_centro[0].d_validx)
yvalid = np.array(sup.rtc._rtc.d_centro[0].d_validy)
rtc.d_centro[0].load_validpos(xvalid, yvalid, xvalid.size)
cmat = sup.rtc.get_command_matrix(0)
rtc.d_control[0].set_cmat(cmat)
rtc.d_control[0].set_gain(sup.config.p_controllers[0].gain)
# frame = sup.wfs.get_wfs_image(0)
# frame /= frame.max()
# rtc.d_centro[0].load_img(frame, frame.shape[0])
# rtc.d_centro[0].calibrate_img()
# rtc.do_centroids(0)
# rtc.do_control(0)

SR32 = []
SR16 = []


def goAhead():
    sup.next()
    sup2.single_next()
    sup.dms.set_command(np.array(sup._sim.rtc.d_control[0].d_voltage))
    sup._sim.raytraceTar(0, "all")
    sup._sim.comp_tar_image(0, compLE=False)
    sup._sim.compStrehl(0)
    SR32.append(sup.get_strehl(0)[0])
    frame = sup2.get_wfs_image(0)
    frame /= frame.max()
    rtc.d_centro[0].load_img(frame, frame.shape[0])
    rtc.d_centro[0].calibrate_img()
    rtc.do_centroids(0)
    rtc.do_control(0)
    rtc.comp_voltage(0)
    sup2.dms.set_command(np.array(rtc.d_control[0].d_voltage))
    sup2._sim.raytraceTar(0, "all")
    sup2._sim.comp_tar_image(0, compLE=False)
    sup2._sim.compStrehl(0)
    SR16.append(sup2.get_strehl(0)[0])


def loop(niter):
    for _ in track(range(niter)):
        goAhead()
