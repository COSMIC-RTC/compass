#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Created on 13 juil. 2017

@author: vdeo
'''

import os
try:
    shesha_dir = os.environ['SHESHA_ROOT']
    os.environ["PATH"] += shesha_dir + '/src'
except KeyError as err:
    raise EnvironmentError(
            "Environment variable 'SHESHA_ROOT' must be defined")

import shesha_config as conf

from shesha import Atmos

import numpy as np


def atmos_init(
        p_atmos: conf.Param_atmos,
        p_tel: conf.Param_tel,
        p_geom: conf.Param_geom,
        p_loop: conf.Param_loop,
        p_wfss=None,
        sensors=None,
        p_target=None,
        rank=0,
        clean=1,
        load={}):

    if not p_geom.isInit:
        raise RuntimeError("Cannot init atmosphere with uninitialized p_geom.")

    # Deleted naga_context : get the singleton

    if p_atmos.r0 == None:  # ?
        p_atmos.r0 = 0.

    # Adjust layers alt using zenith angle
    p_atmos.alt = p_atmos.alt / np.cos(p_geom.zenithangle * conf.DEG2RAD)
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
    patch_diam = (p_geom._n + \
        2 * (max_size * conf.ARCSEC2RAD * p_atmos.alt) / p_atmos.pupixsize + 4
        ).astype(np.int64)
    p_atmos.dim_screens = patch_diam + patch_diam % 2

    # Phase screen speeds
    lin_delta = p_geom.pupdiam / p_tel.diam * p_atmos.windspeed * \
        np.cos(conf.DEG2RAD * p_geom.zenithangle) * p_loop.ittime
    p_atmos.deltax = -lin_delta * np.sin(conf.DEG2RAD * p_atmos.winddir)
    p_atmos.deltay = -lin_delta * np.cos(conf.DEG2RAD * p_atmos.winddir)

    # Fraction normalization
    p_atmos.frac /= np.sum(p_atmos.frac)

    if p_atmos.L0 is None:
        # Set almost infinite L0
        p_atmos.L0 = np.ones(p_atmos.nscreens, dtype=np.float32) * 1e5
    L0_pix = p_atmos.L0 * p_geom.pupdiam / p_tel.diam

    if p_atmos.seeds is None:
        p_atmos.seeds = (
                np.arange(p_atmos.nscreens, dtype=np.int64) + 1) * 1234

    type_target = b"atmos"  # FIXME

    if p_wfss is not None:
        for i in range(len(p_wfss)):
            p_wfs = p_wfss[i]
            if p_wfs.gsalt > 0:
                gsalt = 1. / p_wfs.gsalt
            else:
                gsalt = 0

            if p_wfs.atmos_seen:
                for j in range(p_atmos.nscreens):
                    xoff = (gsalt * p_atmos.alt[j] * p_tel.diam / 2. + \
                            p_wfs.xpos * conf.ARCSEC2RAD * p_atmos.alt[j]) / \
                            p_atmos.pupixsize
                    yoff = (gsalt * p_atmos.alt[j] * p_tel.diam / 2. + \
                            p_wfs.ypos * conf.ARCSEC2RAD * p_atmos.alt[j]) / \
                            p_atmos.pupixsize
                    xoff = xoff + (
                            p_atmos.dim_screens[j] - p_geom.get_n()) / 2.
                    yoff = yoff + (
                            p_atmos.dim_screens[j] - p_geom.get_n()) / 2.
                    sensors.sensors.d_wfs[i].d_gs.add_layer(
                            type_target, p_atmos.alt[j], xoff, yoff)

    return Atmos(
            None, p_atmos.nscreens, p_atmos.r0, L0_pix, p_atmos.pupixsize,
            p_atmos.dim_screens, p_atmos.frac, p_atmos.alt, p_atmos.windspeed,
            p_atmos.winddir, p_atmos.deltax, p_atmos.deltay, p_atmos.seeds,
            rank, clean, load)
