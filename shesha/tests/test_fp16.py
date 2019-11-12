import numpy as np
import naga as ng
import os
from shesha.sutra_wrap import Rtc_FHF as Rtc
from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor
from scipy.ndimage.measurements import center_of_mass
from tqdm import tqdm

precision = 1e-2
sup = Supervisor(
        os.getenv("COMPASS_ROOT") + "/shesha/data/par/par4bench/scao_sh_16x16_8pix.py")
sup.config.p_controllers[0].delay = 0
sup.initConfig()
# sup.singleNext()
sup.openLoop()
sup.closeLoop()
sup2 = Supervisor(
        os.getenv("COMPASS_ROOT") + "/shesha/data/par/par4bench/scao_sh_16x16_8pix.py")
sup2.config.p_controllers[0].delay = 0
sup2.initConfig()
# sup2.singleNext()
sup2.openLoop()
sup2.closeLoop()
rtc = Rtc()
rtc.add_centroider(sup._sim.c, sup.config.p_wfs0._nvalid,
                   sup.config.p_wfs0.npix / 2 - 0.5, sup.config.p_wfs0.pixsize, False, 0, "cog")
rtc.add_controller(sup._sim.c, sup.config.p_wfs0._nvalid, sup.config.p_wfs0._nvalid * 2,
                   sup.config.p_controller0.nactu, sup.config.p_controller0.delay, 0,
                   "generic", idx_centro=np.zeros(1), ncentro=1)
centro = rtc.d_centro[0]
control = rtc.d_control[0]
rtc.d_centro[0].set_npix(sup.config.p_wfs0.npix)
xvalid = np.array(sup._sim.rtc.d_centro[0].d_validx)
yvalid = np.array(sup._sim.rtc.d_centro[0].d_validy)
rtc.d_centro[0].load_validpos(xvalid, yvalid, xvalid.size)
cmat = sup.getCmat(0)
rtc.d_control[0].set_cmat(cmat)
rtc.d_control[0].set_gain(sup.config.p_controller0.gain)
# frame = sup.getWfsImage()
# frame /= frame.max()
# rtc.d_centro[0].load_img(frame, frame.shape[0])
# rtc.d_centro[0].calibrate_img()
# rtc.do_centroids(0)
# rtc.do_control(0)

SR32 = []
SR16 = []

def goAhead():
    sup.singleNext()
    sup2.singleNext()
    sup.setDmShapeFrom(np.array(sup._sim.rtc.d_control[0].d_voltage))
    sup._sim.raytraceTar(0, "all")
    sup._sim.compTarImage(0,compLE=False)
    sup._sim.compStrehl(0)
    SR32.append(sup.getStrehl(0)[0])
    frame = sup2.getWfsImage()
    frame /= frame.max()
    rtc.d_centro[0].load_img(frame, frame.shape[0])
    rtc.d_centro[0].calibrate_img()
    rtc.do_centroids(0)
    rtc.do_control(0)  
    rtc.comp_voltage(0)  
    sup2.setDmShapeFrom(np.array(rtc.d_control[0].d_voltage))
    sup2._sim.raytraceTar(0, "all")
    sup2._sim.compTarImage(0,compLE=False)
    sup2._sim.compStrehl(0)
    SR16.append(sup2.getStrehl(0)[0])

def loop(niter):
    for _ in tqdm(range(niter)):
        goAhead()
