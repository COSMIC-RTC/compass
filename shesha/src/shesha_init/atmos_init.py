#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Created on 13 juil. 2017

@author: vdeo
'''

from naga import naga_context

import shesha_config as conf
from shesha_constants import CONST

from Atmos import Atmos

import numpy as np


def atmos_init(context: naga_context, p_atmos: conf.Param_atmos, p_tel: conf.Param_tel,
               p_geom: conf.Param_geom, ittime=None, p_wfss=None, sensors=None,
               p_target=None):
    """
        TODO: docstring
    """
    if not p_geom.is_init:
        raise RuntimeError("Cannot init atmosphere with uninitialized p_geom.")

    # Deleted naga_context : get the singleton

    if p_atmos.r0 is None:
        p_atmos.r0 = 0.

    if ittime is None:
        ittime = 1.
    # Adjust layers alt using zenith angle
    p_atmos.alt = p_atmos.alt / np.cos(p_geom.zenithangle * CONST.DEG2RAD)
    # Pixel size in meters
    p_atmos.pupixsize = p_tel.diam / p_geom.pupdiam

    # Off-axis wavefront sensors and targets
    # Note : p_wfss is a list of single-WFS configs
    #        but p_target groups several targets
    #        hence different xpos, ypos syntaxes
    norms = [0.]
    if p_wfss is not None:
        norms += [(w.xpos**2 + w.ypos**2)**0.5 for w in p_wfss]
    if p_target is not None:
        norms += list((p_target.xpos**2 + p_target.ypos**2)**0.5)

    max_size = max(norms)

    # Meta-pupil diameter for all layers depending on altitude
    patch_diam = (p_geom._n + 2 *
                  (max_size * CONST.ARCSEC2RAD * p_atmos.alt) / p_atmos.pupixsize +
                  4).astype(np.int64)
    p_atmos.dim_screens = patch_diam + patch_diam % 2

    # Phase screen speeds
    lin_delta = p_geom.pupdiam / p_tel.diam * p_atmos.windspeed * \
        np.cos(CONST.DEG2RAD * p_geom.zenithangle) * ittime
    p_atmos._deltax = -lin_delta * np.sin(CONST.DEG2RAD * p_atmos.winddir)
    p_atmos._deltay = -lin_delta * np.cos(CONST.DEG2RAD * p_atmos.winddir)

    # Fraction normalization
    p_atmos.frac /= np.sum(p_atmos.frac)

    if p_atmos.L0 is None:
        # Set almost infinite L0
        p_atmos.L0 = np.ones(p_atmos.nscreens, dtype=np.float32) * 1e5
    L0_pix = p_atmos.L0 * p_geom.pupdiam / p_tel.diam

    if p_atmos.seeds is None:
        p_atmos.seeds = (np.arange(p_atmos.nscreens, dtype=np.int64) + 1) * 1234

    return Atmos(context, p_atmos.nscreens, p_atmos.r0, L0_pix, p_atmos.pupixsize,
                 p_atmos.dim_screens, p_atmos.frac, p_atmos.alt, p_atmos.windspeed,
                 p_atmos.winddir, p_atmos._deltax, p_atmos._deltay, p_atmos.seeds, 0)
