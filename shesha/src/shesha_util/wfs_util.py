#!/usr/local/bin/python3.6
# encoding: utf-8
'''
Created on 1 aout 2017

@author: fferreira
'''
import numpy as np

import shesha_constants as scons

from shesha_config.PWFS import Param_wfs
from Sensors import Sensors

from . import utilities as util


def comp_new_fstop(
        wfs: Sensors, n: int, p_wfs: Param_wfs, fssize: float, fstop: bytes):
    """
        Compute a new field stop for pyrhr WFS

    :parameters:
        n : (int) : WFS index
        wfs : (Param_wfs) : WFS parameters
        fssize : (float) : field stop size [arcsec]
        fstop : (string) : "square" or "round" (field stop shape)
    """
    fsradius_pixels = int(fssize / p_wfs._qpixsize / 2.)
    if (fstop == scons.FieldStopType.ROUND):
        p_wfs.fstop = fstop
        focmask = util.dist(
                p_wfs._Nfft,
                xc=p_wfs._Nfft / 2. + 0.5,
                yc=p_wfs._Nfft / 2. + 0.5) < (fsradius_pixels)
        # fstop_area = np.pi * (p_wfs.fssize/2.)**2. #UNUSED
    elif (p_wfs.fstop == scons.FieldStopType.SQUARE):
        p_wfs.fstop = fstop
        x, y = util.indices(p_wfs._Nfft)
        x -= (p_wfs._Nfft + 1.) / 2.
        y -= (p_wfs._Nfft + 1.) / 2.
        focmask = (np.abs(x) <= (fsradius_pixels)) * \
            (np.abs(y) <= (fsradius_pixels))
        # fstop_area = p_wfs.fssize**2. #UNUSED
    else:
        msg = "p_wfs " + str(n) + ". fstop must be round or square"
        raise ValueError(msg)

    # pyr_focmask = np.roll(focmask,focmask.shape[0]/2,axis=0)
    # pyr_focmask = np.roll(pyr_focmask,focmask.shape[1]/2,axis=1)
    pyr_focmask = focmask * 1.0  # np.fft.fftshift(focmask*1.0)
    p_wfs._submask = np.fft.fftshift(pyr_focmask).astype(np.float32)
    p_wfs_fssize = fssize
    wfs.set_submask(n, p_wfs._submask)
