import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import naga as ng
from shesha.sutra_wrap import Rtc_UHF as Rtc
supervisor.openLoop()
supervisor.closeLoop()
supervisor._sim.doControl(0)
wfs = supervisor._sim.wfs.d_wfs[0]
wfs.set_fakecam(True)
wfs.set_maxFluxPerPix(int(supervisor.config.p_wfs0._nphotons // 2))
wfs.set_maxPixValue(2**12 - 1)
wfs.comp_image()

rtc = Rtc()
rtc.add_centroider(supervisor._sim.c, supervisor.config.p_wfs0._nvalid,
                   supervisor.config.p_wfs0.npix / 2 - 0.5,
                   supervisor.config.p_wfs0.pixsize, 0, "cog")
rtc.add_controller(supervisor._sim.c, supervisor.config.p_wfs0._nvalid,
                   supervisor.config.p_wfs0._nvalid * 2,
                   supervisor.config.p_controller0.nactu, 1., 0, "generic")
rtc.d_centro[0].set_npix(supervisor.config.p_wfs0.npix)
xvalid = np.array(supervisor._sim.rtc.d_centro[0].d_validx)
yvalid = np.array(supervisor._sim.rtc.d_centro[0].d_validy)
rtc.d_centro[0].load_validpos(xvalid, yvalid, xvalid.size)
cmat = supervisor.getCmat(0)
rtc.d_control[0].set_cmat(cmat)
rtc.d_control[0].set_gain(0.4)
frame = np.array(wfs.d_camimg)
rtc.d_centro[0].load_img(frame, frame.shape[0])
rtc.d_centro[0].calibrate_img()

rtc.do_centroids(0)
slp = ng.array(rtc.d_control[0].d_centroids)
rtc.do_control(0)
com = ng.array(rtc.d_control[0].d_com)

plt.figure()
plt.plot(supervisor.getSlope())
plt.plot(slp.toarray())
plt.legend(["RTC FP32->FP32", "RTC UINT16->FP16"])
plt.title("Slopes")

plt.figure()
plt.plot(supervisor.getCom(0))
plt.plot(com.toarray())
plt.legend(["RTC FP32->FP32", "RTC UINT16->FP16"])
plt.title("Commands")
