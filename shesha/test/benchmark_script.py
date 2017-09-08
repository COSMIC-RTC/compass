import os
import sys
import numpy as np
import naga as ch
import shesha as ao
import time
import datetime
from subprocess import check_output
import pandas
import hdf5_utils as h5u
import platform
import re


SHESHA = os.environ.get('SHESHA_ROOT')
if(SHESHA is None):
    raise EnvironmentError("Environment variable 'SHESHA_ROOT' must be define")

SHESHA_SAVEPATH = SHESHA + "/data/"
PARPATH = SHESHA_SAVEPATH + "par/par4bench/"
BENCH_SAVEPATH = SHESHA_SAVEPATH + "bench-results/"


def get_processor_name():
    command = "cat /proc/cpuinfo"
    all_info = check_output(command, shell=True).strip()
    nb_cpu = 0
    cpu = []
    for line in all_info.split("\n"):
        if "model name" in line:
            cpu.append(re.sub(".*model name.*:", "", line, 1))
            nb_cpu += 1
    return nb_cpu, cpu


def script4bench(param_file, centroider, controller, devices, fwrite=True):
    """

    :parameters:
        param_file: (str) : parameters filename

        centroider: (str) : centroider type

        controller: (str) : controller type
    """

    c = ch.naga_context(devices=np.array(devices, dtype=np.int32))
#     c.set_activeDevice(device)

    timer = ch.naga_timer()

    # times measured
    synctime = 0.
    move_atmos_time = 0.
    t_raytrace_atmos_time = 0.
    t_raytrace_dm_time = 0.
    s_raytrace_atmos_time = 0.
    s_raytrace_dm_time = 0.
    comp_img_time = 0.
    docentroids_time = 0.
    docontrol_time = 0.
    applycontrol_time = 0.

    # reading parfile
    filename = param_file.split('/')[-1]
    param_path = param_file.split(filename)[0]
    sys.path.insert(0, param_path)
    exec("import %s as config" % filename.split(".py")[0])
    sys.path.remove(param_path)

    # set simulation name
    # if(hasattr(config,"simul_name")):
    #     if(config.simul_name is None):
    #         simul_name=""
    #     else:
    #         simul_name=config.simul_name
    # else:
    #     simul_name=""
    simul_name = ""
    matricesToLoad = {}
    config.p_centroiders[0].set_type(centroider)

    if(centroider == "tcog"):
        config.p_centroiders[0].set_thresh(9.)
    elif(centroider == "bpcog"):
        config.p_centroiders[0].set_nmax(16)
    elif(centroider == "geom"):
        config.p_centroiders[0].set_type("cog")
    elif(centroider == "wcog"):
        config.p_centroiders[0].set_type_fct("gauss")
        config.p_centroiders[0].set_width(2.0)
    elif(centroider == "corr"):
        config.p_centroiders[0].set_type_fct("gauss")
        config.p_centroiders[0].set_width(2.0)

    config.p_controllers[0].set_type(controller)
    if(controller == "modopti"):
        config.p_controllers[0].set_type("ls")
        config.p_controllers[0].set_modopti(1)

    config.p_loop.set_niter(2000)
    if(simul_name == ""):
        clean = 1
    else:
        clean = 0
        param_dict = h5u.params_dictionary(config)
        matricesToLoad = h5u.checkMatricesDataBase(
            os.environ["SHESHA_ROOT"] + "/data/", config, param_dict)

    ch.threadSync()
    timer.start()
    ch.threadSync()
    synctime = timer.stop()
    timer.reset()

    # init system
    timer.start()
    wfs, tel = ao.wfs_init(config.p_wfss, config.p_atmos, config.p_tel,
                           config.p_geom, config.p_target, config.p_loop, config.p_dms)
    ch.threadSync()
    wfs_init_time = timer.stop() - synctime
    timer.reset()

    timer.start()
    atm = ao.atmos_init(c, config.p_atmos, config.p_tel, config.p_geom, config.p_loop,
                        config.p_wfss, wfs, config.p_target, clean=clean, load=matricesToLoad)
    ch.threadSync()
    atmos_init_time = timer.stop() - synctime
    timer.reset()

    timer.start()
    dms = ao.dm_init(config.p_dms, config.p_wfss,
                     wfs, config.p_geom, config.p_tel)
    ch.threadSync()
    dm_init_time = timer.stop() - synctime
    timer.reset()

    timer.start()
    target = ao.target_init(c, tel, config.p_target, config.p_atmos,
                            config.p_geom, config.p_tel, config.p_dms)
    ch.threadSync()
    target_init_time = timer.stop() - synctime
    timer.reset()

    timer.start()
    rtc = ao.rtc_init(tel, wfs, config.p_wfss, dms, config.p_dms, config.p_geom, config.p_rtc, config.p_atmos,
                      atm, config.p_tel, config.p_loop, clean=clean, simul_name=simul_name, load=matricesToLoad)
    ch.threadSync()
    rtc_init_time = timer.stop() - synctime
    timer.reset()

    print "... Done with inits !"
    # h5u.validDataBase(os.environ["SHESHA_ROOT"]+"/data/",matricesToLoad)

    strehllp = []
    strehlsp = []
############################################################
#                  _         _
#                 (_)       | |
#  _ __ ___   __ _ _ _ __   | | ___   ___  _ __
# | '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \
# | | | | | | (_| | | | | | | | (_) | (_) | |_) |
# |_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/
#                                         | |
#                                         |_|
###########################################################
    if(controller == "modopti"):
        for zz in xrange(2048):
            atm.move_atmos()

    for cc in xrange(config.p_loop.niter):
        ch.threadSync()
        timer.start()
        atm.move_atmos()
        ch.threadSync()
        move_atmos_time += timer.stop() - synctime
        timer.reset()

        if(config.p_controllers[0].type_control != "geo"):
            if((config.p_target is not None) and (rtc is not None)):
                for i in xrange(config.p_target.ntargets):
                    timer.start()
                    target.atmos_trace(i, atm, tel)
                    ch.threadSync()
                    t_raytrace_atmos_time += timer.stop() - synctime
                    timer.reset()

                    if(dms is not None):
                        timer.start()
                        target.dmtrace(i, dms)
                        ch.threadSync()
                        t_raytrace_dm_time += timer.stop() - synctime
                        timer.reset()

            if(config.p_wfss is not None and wfs is not None):
                for i in xrange(len(config.p_wfss)):
                    timer.start()
                    wfs.sensors_trace(i, "atmos", tel, atm)
                    ch.threadSync()
                    s_raytrace_atmos_time += timer.stop() - synctime
                    timer.reset()

                    if(not config.p_wfss[i].openloop and dms is not None):
                        timer.start()
                        wfs.sensors_trace(i, "dm", dms=dms)
                        ch.threadSync()
                        s_raytrace_dm_time += timer.stop() - synctime
                        timer.reset()

                    timer.start()
                    wfs.sensors_compimg(i)
                    ch.threadSync()
                    comp_img_time += timer.stop() - synctime
                    timer.reset()

            if(config.p_rtc is not None and
               rtc is not None and
               config.p_wfss is not None and
               wfs is not None):
                if(centroider == "geom"):
                    timer.start()
                    rtc.docentroids_geom(0)
                    ch.threadSync()
                    docentroids_time += timer.stop() - synctime
                    timer.reset()
                else:
                    timer.start()
                    rtc.docentroids(0)
                    ch.threadSync()
                    docentroids_time += timer.stop() - synctime
                    timer.reset()

            if(dms is not None):
                timer.start()
                rtc.docontrol(0)
                ch.threadSync()
                docontrol_time += timer.stop() - synctime
                timer.reset()

                timer.start()
                rtc.applycontrol(0, dms)
                ch.threadSync()
                applycontrol_time += timer.stop() - synctime
                timer.reset()

        else:
            if(config.p_target is not None and target is not None):
                for i in xrange(config.p_target.ntargets):
                    timer.start()
                    target.atmos_trace(i, atm, tel)
                    ch.threadSync()
                    t_raytrace_atmos_time += timer.stop() - synctime
                    timer.reset()

                    if(dms is not None):
                        timer.start()
                        rtc.docontrol_geo(0, dms, target, i)
                        ch.threadSync()
                        docontrol_time += timer.stop() - synctime
                        timer.reset()

                        timer.start()
                        rtc.applycontrol(0, dms)
                        ch.threadSync()
                        applycontrol_time += timer.stop() - synctime
                        timer.reset()

                        timer.start()
                        target.dmtrace(i, dms)
                        ch.threadSync()
                        t_raytrace_dm_time += timer.stop() - synctime
                        timer.reset()

        strehltmp = target.get_strehl(0)
        strehlsp.append(strehltmp[0])
        if(cc > 50):
            strehllp.append(strehltmp[1])

    print "\n done with simulation \n"
    print "\n Final strehl : \n", strehllp[len(strehllp) - 1]
###################################################################
#  _   _
# | | (_)
# | |_ _ _ __ ___   ___ _ __ ___
# | __| | '_ ` _ \ / _ \ '__/ __|
# | |_| | | | | | |  __/ |  \__ \
#  \__|_|_| |_| |_|\___|_|  |___/
###################################################################

    move_atmos_time /= config.p_loop.niter / 1000.
    t_raytrace_atmos_time /= config.p_loop.niter / 1000.
    t_raytrace_dm_time /= config.p_loop.niter / 1000.
    s_raytrace_atmos_time /= config.p_loop.niter / 1000.
    s_raytrace_dm_time /= config.p_loop.niter / 1000.
    comp_img_time /= config.p_loop.niter / 1000.
    docentroids_time /= config.p_loop.niter / 1000.
    docontrol_time /= config.p_loop.niter / 1000.
    applycontrol_time /= config.p_loop.niter / 1000.

    time_per_iter = move_atmos_time + t_raytrace_atmos_time +\
        t_raytrace_dm_time + s_raytrace_atmos_time +\
        s_raytrace_dm_time + comp_img_time +\
        docentroids_time + docontrol_time +\
        applycontrol_time

###########################################################################
#  _         _  __ _____
# | |       | |/ _| ____|
# | |__   __| | |_| |__    ___  __ ___   _____
# | '_ \ / _` |  _|___ \  / __|/ _` \ \ / / _ \
# | | | | (_| | |  ___) | \__ \ (_| |\ V /  __/
# |_| |_|\__,_|_| |____/  |___/\__,_| \_/ \___|
###############################################################################

    if(config.p_wfss[0].gsalt > 0):
        stype = "lgs "
    else:
        stype = "ngs "

    if(config.p_wfss[0].gsmag > 3):
        stype += "noisy "

    stype += config.p_wfss[0].type_wfs

    if(controller == "modopti"):
        G = np.mean(rtc.get_mgain(0))
    else:
        G = 0.

    date = datetime.datetime.now()
    date = [date.year, date.month, date.day]

    version = ao.__version__

    # version=str(check_output(["svnversion",os.getenv("COMPASS_ROOT")]).replace("\n",""))
    hostname = check_output("hostname").replace("\n", "")
    nb_cpu, cpu = get_processor_name()
    imat = rtc.get_imat(0)
    keys_dict = {"date": date,
                 "simulname": config.simul_name,
                 "hostname": hostname,
                 "ndevices": c.get_ndevice(),
                 "device": c.get_device_names()[0],
                 "cuda_version": c.get_cudaRuntimeGetVersion(),
                 "magma_version": c.get_magma_info(),
                 "platform": platform.platform(),
                 "ncpu": nb_cpu,
                 "processor": cpu[0],
                 "tel.diam": config.p_tel.diam,
                 "sensor_type": config.p_wfss[0].type_wfs,
                 "nslopes": imat.shape[0],
                 "nactus": imat.shape[1],
                 "LGS": config.p_wfss[0].gsalt > 0,
                 "noisy": config.p_wfss[0].gsmag > 3,
                 "nxsub": config.p_wfss[0].nxsub,
                 "npix": config.p_wfss[0].npix,
                 "nphotons": config.p_wfss[0]._nphotons,
                 "controller": controller,
                 "centroider": centroider,
                 "finalSRLE": strehllp[len(strehllp) - 1],
                 "rmsSRLE": np.std(strehllp),
                 "wfs_init": wfs_init_time,
                 "atmos_init": atmos_init_time,
                 "dm_init": dm_init_time,
                 "target_init": target_init_time,
                 "rtc_init": rtc_init_time,
                 "move_atmos": move_atmos_time,
                 "target_trace_atmos": t_raytrace_atmos_time,
                 "target_trace_dm": t_raytrace_dm_time,
                 "sensor_trace_atmos": s_raytrace_atmos_time,
                 "sensor_trace_dm": s_raytrace_dm_time,
                 "comp_img": comp_img_time,
                 "docentroids": docentroids_time,
                 "docontrol": docontrol_time,
                 "applycontrol": applycontrol_time,
                 "iter_time": time_per_iter,
                 "Avg.gain": G,
                 "ResidualPhase:": target.get_phase(0)}

    store = pandas.HDFStore(BENCH_SAVEPATH + "benchmarks.h5")
    try:
        df = store.get(version)
    except KeyError:
        df = pandas.DataFrame(columns=keys_dict.keys(), dtype=object)

    ix = len(df.index)

    if(fwrite):
        print "writing files"
        for i in keys_dict.keys():
            df.loc[ix, i] = keys_dict[i]
        store.put(version, df)
        store.close()
#############################################################
#                 _
#                | |
#   ___ _ __   __| |
#  / _ \ '_ \ / _` |
# |  __/ | | | (_| |
#  \___|_| |_|\__,_|
#############################################################


if(len(sys.argv) < 4 or len(sys.argv) > 6):
    error = "wrong number of argument. Got %d (expect 4)\ncommande line should be: 'python benchmark_script.py <filename> <centroider> <controller>" % len(
        sys.argv)
    raise StandardError(error)

filename = PARPATH + sys.argv[1]
centroider = sys.argv[2]
controller = sys.argv[3]
device = 5
fwrite = True
if(len(sys.argv) > 4):
    devices = []
    if(len(sys.argv[4]) > 1):
        for k in range(len(sys.argv[4])):
            devices.append(int(sys.argv[4][k]))
    else:
        devices.append(int(sys.argv[4]))
if (len(sys.argv) == 6):
    fwrite = int(sys.argv[5])

script4bench(filename, centroider, controller, devices, fwrite)
