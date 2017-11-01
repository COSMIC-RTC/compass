# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 08:53:38 2015

@author: fferreira
"""

import cProfile
import pstats as ps

import sys, os
import numpy as np
import naga as ch
import shesha as ao
import time
import matplotlib.pyplot as pl

#get parameters from file
param_file = os.environ["SHESHA"] + "/data/par/p_aliasing_covariance.py"
filename = param_file.split('/')[-1]
param_path = param_file.split(filename)[0]
sys.path.insert(0, param_path)
exec("import %s as config" % filename.split(".py")[0])
sys.path.remove(param_path)

print("param_file is", param_file)

if (hasattr(config, "simul_name")):
    if (config.simul_name is None):
        simul_name = ""
    else:
        simul_name = config.simul_name
else:
    simul_name = ""

if (simul_name == b""):
    clean = 1
else:
    clean = 0

#initialisation:
#   context
c = ch.naga_context()
c.set_activeDevice(0)

#    wfs
print("->wfs")
wfs = ao.wfs_init(config.p_wfss, config.p_atmos, config.p_tel, config.p_geom,
                  config.p_target, config.p_loop, 1, 0, config.p_dms)

#   atmos
print("->atmos")
atm = ao.atmos_init(c, config.p_atmos, config.p_tel, config.p_geom, config.p_loop,
                    config.p_wfss, wfs, config.p_target, rank=0)

#   dm
print("->dm")
dms = ao.dm_init(config.p_dms, config.p_wfss[0], config.p_geom, config.p_tel)

#   target
print("->target")
tar = ao.target_init(c, config.p_target, config.p_atmos, config.p_geom, config.p_tel,
                     config.p_dms)

print("->rtc")
#   rtc
rtc = ao.rtc_init(wfs, config.p_wfss, dms, config.p_dms, config.p_geom, config.p_rtc,
                  config.p_atmos, atm, config.p_tel, config.p_loop, tar, config.p_target,
                  clean=clean, simul_name=simul_name)

print("====================")
print("init done")
print("====================")
print("objects initialzed on GPU:")
print("--------------------------------------------------------")
print(atm)
print(wfs)
print(dms)
print(tar)
print(rtc)

print("----------------------------------------------------")
print("iter# | S.E. SR | L.E. SR | Est. Rem. | framerate")
print("----------------------------------------------------")

###############################################################################################################################################################
# Here, we want to get centroids obtain from a SH wfs seeing an atmos screen
# and compare it to the centroids obtain from a wfs seeing the parallel component
# and a wfs seeing the orthogonal component
# For this, 1 wfs is seeing the atmos in open loop, an other sees atmos+dm after
# it shape has been computed with a geo controller, and the last one just sees
# the dm
###############################################################################################################################################################


def loop(n, r, e):
    """ Here, we want to get centroids obtain from a SH wfs seeing an atmos screen
    and compare it to the centroids obtain from a wfs seeing the parallel component
    and a wfs seeing the orthogonal component
    For this, 1 wfs is seeing the atmos in open loop, an other sees atmos+dm after
    it shape has been computed with a geo controller, and the last one just sees
    the dm
    """
    pl.ion()
    reset_com = np.zeros(rtc.getCom(0).size, dtype=np.float32)
    for i in range(n):
        for j in range(100):  # to get very different realisation of input phase screen
            atm.move_atmos()

        rtc.setCom(0, reset_com)  # Reset the command vector
        # Remove low frequency component with geo controller
        tar.atmos_trace(0, atm)
        rtc.docontrol_geo(1, dms, tar, 0)
        rtc.applycontrol(1, dms)
        wfs.sensors_trace(0, "all", atm,
                          dms)  # Aliasing component since no noise in system

        wfs.sensors_compimg(0)
        rtc.docentroids(0)
        rtc.docontrol(0)
        r[:, i] = rtc.getCom(0)

        rtc.setCom(0, reset_com)
        wfs.sensors_trace(0, "atmos", atm)
        wfs.sensors_compimg(0)
        rtc.docentroids(0)
        rtc.docontrol(0)
        e[:, i] = rtc.getCom(0)

        print(" Recording... %d%%\r" % (i * 100 / n), end=' ')
    print("Recorded")


def covariance_analysis(r):
    pass


#############################################################################################################
################################################# MAIN ANALYSIS #############################################
#############################################################################################################
#Allocate buffers
niter = 10000
r = np.zeros((rtc.getCom(0).size, niter))
e = np.copy(r)
loop(niter, r, e)
