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
