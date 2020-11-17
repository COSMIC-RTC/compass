import numpy as np
import naga as ng
import os
from shesha.supervisor.components import RtcStandalone
from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor
from scipy.ndimage.measurements import center_of_mass
from shesha.config import ParamConfig

precision = 1e-2

config = ParamConfig(os.getenv("COMPASS_ROOT") +
        "/shesha/tests/pytest/par/test_pyrhr.py")
sup = Supervisor(config)
sup.next()
sup.rtc.open_loop(0)
sup.rtc.close_loop(0)
sup.rtc._rtc.do_control(0)
rtc = RtcStandalone(sup.context, sup.config, 1, [sup.config.p_wfss[0]._nvalid], sup.config.p_controllers[0].nactu,
                    ["maskedpix"], [sup.config.p_controllers[0].delay], [0], [1])
centro = rtc._rtc.d_centro[0]
control = rtc._rtc.d_control[0]
rtc._rtc.d_centro[0].set_npix(sup.config.p_wfss[0].npix)
xvalid = np.array(sup.rtc._rtc.d_centro[0].d_validx)
yvalid = np.array(sup.rtc._rtc.d_centro[0].d_validy)
rtc._rtc.d_centro[0].load_validpos(xvalid, yvalid, xvalid.size)
cmat = sup.rtc.get_command_matrix(0)
rtc._rtc.d_control[0].set_cmat(cmat)
rtc._rtc.d_control[0].set_gain(sup.config.p_controllers[0].gain)
frame = sup.wfs.get_wfs_image(0)
frame /= frame.max()
rtc._rtc.d_centro[0].load_img(frame, frame.shape[0])
rtc._rtc.d_centro[0].calibrate_img()

rtc._rtc.do_centroids(0)
slp = ng.array(rtc._rtc.d_control[0].d_centroids)
rtc._rtc.do_control(0)
com = ng.array(rtc._rtc.d_control[0].d_com)

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
        slopes[k] = binimg[xvalid[k], yvalid[k]] / psum
    assert (relative_array_error(ng.array(control.d_centroids).toarray(), slopes) <
            precision)
            
def test_calibrate_img_validPix():
    centro.calibrate_img_validPix()
    validx = np.array(centro.d_validx)
    validy = np.array(centro.d_validy)
    valid_mask = frame*0
    valid_mask[validx, validy] = 1
    imgCal = (frame - dark) * flat * valid_mask
    assert (relative_array_error(np.array(centro.d_img), imgCal) < precision)
