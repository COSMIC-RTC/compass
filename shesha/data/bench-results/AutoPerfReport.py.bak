# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 15:41:57 2016

@author: fferreira

Module pylatex needed:
    pip install pylatex

LaTeX distribution needed + lmodern package

"""
import os
import shesha
import pandas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pylatex import Document, Section, Subsection, Tabular, Figure, NoEscape, Command
from subprocess import check_output

matplotlib.use('Agg')
SHESHA = os.environ.get('SHESHA_ROOT')
BENCH_SAVEPATH = SHESHA + "data/bench-results/"
plt.ion()


def depouillePerf(filename, version=None, mode="profile"):

    store = pandas.HDFStore(filename)
    if(version is None):
        version = shesha.__version__
    try:
        df = store.get(version)
    except KeyError:
        print "No results for git version : " + version + ", taking " + store.keys()[-1] + " version"
        version = store.keys()[-1]
        df = store.get(version)

    simulnames = df["simulname"].values
    wfs_type = np.unique(df["sensor_type"].values)
    nxsub = np.unique(df["nxsub"].values)
    npix = np.unique(df["npix"].values)
    controllers = np.unique(df["controller"].values)
    centroiders = np.unique(df["centroider"].values)

    times = ["move_atmos", "target_trace_atmos", "target_trace_dm", "sensor_trace_atmos",
             "sensor_trace_dm", "comp_img", "docentroids", "docontrol", "applycontrol"]

    plt.figure()

    colors = ["blue", "green", "red", "yellow", "orange",
              "cyan", "purple", "magenta", "darkcyan"]
    width = 0.9 / len(simulnames)
    ind = 0
    cc = 0
    pos = []
    lab = []
    for indx in df.index:
        ce = df.loc[indx, "centroider"]
        co = df.loc[indx, "controller"]
        ccc = 0
        if(mode == "full"):
            plt.barh(cc, df.loc[indx, times[0]], width, color=colors[ccc])
            timeb = df.loc[indx, times[0]]
            ccc += 1
            for i in times[1:]:
                plt.barh(cc, df.loc[indx, i], width,
                         color=colors[ccc], left=timeb)
                timeb += df.loc[indx, i]
                ccc += 1
        elif(mode == "profile"):
            tottime = 0
            for i in times:
                tottime += df.loc[indx, i]
            plt.barh(cc, df.loc[indx, times[0]] /
                     tottime * 100, width, color=colors[ccc])
            timeb = df.loc[indx, times[0]] / tottime * 100
            ccc += 1
            for i in times[1:]:
                plt.barh(cc, df.loc[indx, i] / tottime * 100,
                         width, color=colors[ccc], left=timeb)
                if(df.loc[indx, i] / tottime * 100 > 10.):
                    plt.text(timeb + df.loc[indx, i] / tottime * 100 / 2, cc +
                             width / 2., '%d' % int(df.loc[indx, i] / tottime * 100) + " %")
                timeb += df.loc[indx, i] / tottime * 100
                ccc += 1
        elif(mode == "framerate"):
            plt.barh(cc, 1. / df.loc[indx, "iter_time"]
                     * 1000., width, color=colors[ind])
            plt.text(1. / df.loc[indx, "iter_time"] * 1000. / 2., cc + width /
                     2., '%.1f' % float(1. / df.loc[indx, "iter_time"] * 1000.))
            ccc += 1

        pos.append(cc + width / 2.)
        infos = df.loc[indx, "simulname"].split('_')
        pyr = df.loc[indx, "simulname"].split('pyr')
        if(len(pyr) > 1):
            ssp = pyr[1]
            wfs = 'pyramid'
        else:
            ssp = infos[2]
            wfs = 'sh ' + infos[3]

        lab.append("%s %s \n %s \n %s ssp " % (infos[0], infos[1], wfs, ssp))
        ind += 1
        cc += 1. / len(simulnames)

    plt.yticks(pos, lab)
    if(mode == "full"):
        plt.title("Execution times")
        plt.xlabel("Execution time (ms)")
        plt.legend(times)
    elif(mode == "profile"):
        plt.title("Execution profile")
        plt.xlabel("Occupation time (%)")
        plt.legend(times)
    elif(mode == "framerate"):
        plt.title("Framerate")
        plt.xlabel("Framerate (frames/s)")

    store.close()


store = pandas.HDFStore("benchmarks.h5")
version = store.keys()[-1]
df = store.get(version)
simulnames = df["simulname"].values
date = df["date"].values[0]

geometry_options = {"tmargin": "1cm", "lmargin": "10cm"}
doc = Document(geometry_options=geometry_options)

doc.preamble.append(Command(
    'title', 'COMPASS Performance Report on %s' % str(df["device"].values[0])))
doc.preamble.append(Command('author', 'F. Ferreira'))
doc.preamble.append(Command('date', '%d/%d/%d' % (date[2], date[1], date[0])))
doc.append(NoEscape(r'\maketitle'))


with doc.create(Section('Environment')):
    doc.append(
        'Simulations have been performed within the following environment:\n ')
    with doc.create(Tabular('|l|c|')) as table:
        table.add_hline()
        table.add_row(("Platform", str(df["platform"].values[0])))
        table.add_hline()
        table.add_row(("Processors", str(
            str(df["ncpu"].values[0]) + "x " + df["processor"].values[0])))
        table.add_hline()
        table.add_row(("GPU devices", str(
            str(df["ndevices"].values[0]) + "x " + df["device"].values[0])))
        table.add_hline()
        table.add_row(("CUDA Version", str(df["cuda_version"].values[0])))
        table.add_hline()
        table.add_row(("MAGMA Version", str(df["magma_version"].values[0])))
        table.add_hline()

with doc.create(Section('Simulation parameters')):
    doc.append('We ran several simulation cases: SCAO cases with Shack-Hartman sensor and pyramid sensor,\
         MCAO cases with tomographic controller. The following tables show the relevant parameters for performance study.\n ')
    for indx in df.index:
        with doc.create(Subsection(df.loc[indx, "simulname"])):
            with doc.create(Tabular('|l|c|')) as table:
                table.add_hline()
                table.add_row(
                    ("Simulation case", str(df.loc[indx, "simulname"])))
                table.add_hline()
                table.add_row(
                    ("Telescope diameter [m]", str(df.loc[indx, "tel.diam"])))
                table.add_hline()
                table.add_row(("WFS type", str(df.loc[indx, "sensor_type"])))
                table.add_hline()
                table.add_row(("Number of subaps", str(df.loc[indx, "nxsub"])))
                table.add_hline()
                table.add_row(("Number of pixels/subap",
                               str(df.loc[indx, "npix"])))
                table.add_hline()
                table.add_row(
                    ("Controller type", str(df.loc[indx, "controller"])))
                table.add_hline()
                table.add_row(
                    ("Centroider type", str(df.loc[indx, "centroider"])))
                table.add_hline()
                # table.add_row(("Total number of slopes",str(df.loc[indx,"nslopes"])))
                # table.add_hline()
                # table.add_row(("Total number of actuators",str(df.loc[indx,"nactus"])))
                # table.add_hline()
            doc.append('\n')

with doc.create(Section('Results')):
    depouillePerf("benchmarks.h5")
    with doc.create(Subsection("Execution profile")):
        with doc.create(Figure(position='!htbp')) as plot:
            plot.add_plot(width=NoEscape(r'1\textwidth'))
    depouillePerf("benchmarks.h5", mode="framerate")
    with doc.create(Subsection("Framerate")):
        with doc.create(Figure(position='!htbp')) as plot:
            plot.add_plot(width=NoEscape(r'1\textwidth'))

doc.generate_pdf("Test_report", clean_tex=False)
