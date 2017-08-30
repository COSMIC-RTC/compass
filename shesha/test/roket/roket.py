"""
Created on Wed Apr 27 09:28:23 2016

@author: fferreira
"""

import cProfile
import pstats as ps

import sys, os
import numpy as np
from naga import naga_context
from shesha_sim.simulator import Simulator
from shesha_ao.tomo import create_nact_geom
from shesha_ao.basis import compute_Btt, compute_cmat_with_Btt
import shesha_constants as scons
import time
import matplotlib.pyplot as pl
pl.ion()
import hdf5_utils as h5u
import pandas
from scipy.sparse import csr_matrix


class Roket(Simulator):

    def __init__(self, str=None, N_preloop=1000, gamma=1.):
        super().__init__(str)
        self.N_preloop = N_preloop
        self.gamma = gamma

    def init_sim(self):
        super().init_sim()
        self.iter_number = 0
        self.n = self.config.p_loop.niter + self.N_preloop
        self.nfiltered = int(self.config.p_controllers[0].maxcond)
        self.nactu = self.rtc.get_com(0).size
        self.com = np.zeros((self.n, self.nactu), dtype=np.float32)
        self.noise_com = np.zeros((self.n, self.nactu), dtype=np.float32)
        self.alias_wfs_com = np.copy(self.noise_com)
        self.wf_com = np.copy(self.noise_com)
        self.tomo_com = np.copy(self.noise_com)
        self.trunc_com = np.copy(self.noise_com)
        self.H_com = np.copy(self.noise_com)
        self.mod_com = np.copy(self.noise_com)
        self.bp_com = np.copy(self.noise_com)
        self.fit = np.zeros(self.n)
        self.psf_ortho = self.tar.get_image(0, b'se') * 0.
        self.Ee = np.copy(self.noise_com)
        self.Ff = np.copy(self.Ee)
        #gamma = 1.0
        self.config.p_loop.set_niter(self.n)
        self.IFpzt = self.rtc.get_IFsparse(1)
        self.TT = self.rtc.get_IFtt(1)

        self.Btt, self.P = compute_Btt(self.IFpzt.T, self.TT)
        self.rtc.load_Btt(1, self.Btt.dot(self.Btt.T))
        self.D = self.rtc.get_imat(0)
        compute_cmat_with_Btt(self.rtc, self.Btt, self.nfiltered)
        self.cmat = self.rtc.get_cmat(0)
        self.D = self.rtc.get_imat(0)
        self.RD = np.dot(self.cmat, self.D)
        self.gRD = np.identity(
                self.RD.shape[0]
        ) - self.config.p_controllers[0].gain * self.gamma * self.RD

        self.Nact = create_nact_geom(self.config.p_dms[0])

    def next(self, **kwargs):
        """
            function next
            Iterates the AO loop, with optional parameters

        :param move_atmos (bool): move the atmosphere for this iteration, default: True
        :param nControl (int): Controller number to use, default 0 (single control configurations)
        :param tar_trace (None or list[int]): list of targets to trace. None equivalent to all.
        :param wfs_trace (None or list[int]): list of WFS to trace. None equivalent to all.
        :param apply_control (bool): (optional) if True (default), apply control on DMs
        """
        super().next(apply_control=False)
        self.error_breakdown()
        self.rtc.apply_control(0, self.dms)
        self.iter_number += 1

    def loop(self, monitoring_freq=100, **kwargs):
        """
            TODO: docstring
        """
        print("-----------------------------------------------------------------")
        print("iter# | SE SR | LE SR | FIT SR | PH SR | ETR (s) | Framerate (Hz)")
        print("-----------------------------------------------------------------")
        t0 = time.time()
        for i in range(self.n):
            self.next(**kwargs)
            if ((i + 1) % monitoring_freq == 0):
                framerate = (i + 1) / (time.time() - t0)
                strehltmp = self.tar.get_strehl(0)
                etr = (self.n - i) / framerate
                print("%d \t %.2f \t  %.2f\t %.2f \t %.2f \t    %.1f \t %.1f" %
                      (i + 1, strehltmp[0], strehltmp[1], np.exp(-strehltmp[2]),
                       np.exp(-strehltmp[3]), etr, framerate))
        t1 = time.time()

        print(" loop execution time:", t1 - t0, "  (", self.n, "iterations), ",
              (t1 - t0) / self.n, "(mean)  ", self.n / (t1 - t0), "Hz")
        self.SR2 = np.exp(-self.tar.get_strehl(0, comp_strehl=False)[3])
        self.SR = self.tar.get_strehl(0, comp_strehl=False)[1]

    def error_breakdown(self):
        """
        Compute the error breakdown of the AO simulation. Returns the error commands of
        each contributors. Suppose no delay (for now) and only 2 controllers : the main one, controller #0, (specified on the parameter file)
        and the geometric one, controller #1 (automatically added if error_budget is asked in the parameter file)
        Commands are computed by applying the loop filter on various kind of commands : (see schema_simulation_budget_erreur_v2)

            - Ageom : Aliasing contribution on WFS direction
                Obtained by computing commands from DM orthogonal phase (projection + slopes_geom)

            - B : Projection on the target direction
                Obtained as the commmands output of the geometric controller

            - C : Wavefront
                Obtained by computing commands from DM parallel phase (RD*B)

            - E : Wavefront + aliasing + ech/trunc + tomo
                Obtained by performing the AO loop iteration without noise on the WFS

            - F : Wavefront + aliasing + tomo
                Obtained by performing the AO loop iteration without noise on the WFS and using phase deriving slopes

            - G : tomo

        Note : rtc.get_err returns to -CMAT.slopes
        """
        g = self.config.p_controllers[0].gain
        Dcom = self.rtc.get_com(0)
        Derr = self.rtc.get_err(0)
        self.com[self.iter_number, :] = Dcom
        tarphase = self.tar.get_phase(0)
        ###########################################################################
        ## Noise contribution
        ###########################################################################
        if (self.config.p_wfss[0].type_wfs == scons.WFSType.SH):
            ideal_bincube = self.wfs.get_bincube_not_noisy(0)
            bincube = self.wfs.get_bincube(0)
            if (self.config.p_centroiders[0].type_centro == scons.CentroiderType.
                        TCOG):  # Select the same pixels with or without noise
                invalidpix = np.where(bincube <= self.config.p_centroiders[0].thresh)
                ideal_bincube[invalidpix] = 0
                self.rtc.set_thresh(0, -1e16)
            self.wfs.set_bincube(0, ideal_bincube)
        elif (self.config.p_wfss[0].type_wfs == scons.centroiderType.PYRHR):
            ideal_pyrimg = self.wfs.get_binimg_not_noisy(0)
            self.wfs.set_pyrimg(0, ideal_pyrimg)

        self.rtc.do_centroids(0)
        if (self.config.p_centroiders[0].type_centro == scons.CentroiderType.TCOG):
            self.rtc.set_thresh(0, config.p_centroiders[0].thresh)

        self.rtc.do_control(0)
        E = self.rtc.get_err(0)
        self.Ee[self.iter_number, :] = E
        # Apply loop filter to get contribution of noise on commands
        if (self.iter_number + 1 < self.config.p_loop.niter):
            self.noise_com[self.iter_number + 1, :] = self.gRD.dot(
                    self.noise_com[self.iter_number, :]) + g * (Derr - E)
        ###########################################################################
        ## Sampling/truncature contribution
        ###########################################################################
        self.rtc.do_centroids_geom(0)
        self.rtc.do_control(0)
        F = self.rtc.get_err(0)
        self.Ff[self.iter_number, :] = F
        # Apply loop filter to get contribution of sampling/truncature on commands
        if (self.iter_number + 1 < self.config.p_loop.niter):
            self.trunc_com[self.iter_number + 1, :] = self.gRD.dot(
                    self.trunc_com[self.iter_number, :]) + g * (E - self.gamma * F)

        ###########################################################################
        ## Aliasing contribution on WFS direction
        ###########################################################################
        self.rtc.do_control_geo_on_wfs(1, self.dms, self.wfs, 0)
        self.rtc.apply_control(1, self.dms)
        for w in range(len(self.config.p_wfss)):
            self.wfs.raytrace(w, b"dm", self.tel, self.atm, self.dms)
        """
            wfs.sensors_compimg(0)
        if(config.p_wfss[0].type_wfs == scons.WFSType.SH):
            ideal_bincube = wfs.get_bincubeNotNoisy(0)
            bincube = wfs.get_bincube(0)
            if(config.p_centroiders[0].type_centro == scons.CentroiderType.TCOG): # Select the same pixels with or without noise
                invalidpix = np.where(bincube <= config.p_centroiders[0].thresh)
                ideal_bincube[self.iter_numbernvalidpix] = 0
                rtc.setthresh(0,-1e16)
            wfs.set_bincube(0,ideal_bincube)
        elif(config.p_wfss[0].type_wfs == scons.centroiderType.PYRHR):
            ideal_pyrimg = wfs.get_binimg_notnoisy(0)
            wfs.set_pyrimg(0,ideal_pyrimg)
        """
        self.rtc.do_centroids_geom(0)
        self.rtc.do_control(0)
        Ageom = self.rtc.get_err(0)
        if (self.iter_number + 1 < self.config.p_loop.niter):
            self.alias_wfs_com[self.iter_number + 1, :] = self.gRD.dot(
                    self.alias_wfs_com[self.iter_number, :]) + self.gamma * g * (
                            Ageom)  # - (E-F))

        ###########################################################################
        ## Wavefront + filtered modes reconstruction
        ###########################################################################
        self.tar.raytrace(0, b"atmos", self.tel, self.atm)
        self.rtc.do_control_geo(1, self.dms, self.tar, 0)
        B = self.rtc.get_com(1)

        ###########################################################################
        ## Fitting
        ###########################################################################
        self.rtc.apply_control(1, self.dms)
        self.tar.raytrace(0, b"dm", self.tel, dms=self.dms, do_phase_var=0)
        self.fit[self.iter_number] = self.tar.get_strehl(0, comp_strehl=False)[2]
        if (self.iter_number >= self.N_preloop):
            self.psf_ortho += self.tar.get_image(0, b'se')

        ###########################################################################
        ## Filtered modes error & Commanded modes
        ###########################################################################
        modes = self.P.dot(B)
        modes_filt = modes.copy() * 0.
        modes_filt[-self.nfiltered - 2:-2] = modes[-self.nfiltered - 2:-2]
        self.H_com[self.iter_number, :] = self.Btt.dot(modes_filt)
        modes[-self.nfiltered - 2:-2] = 0
        self.mod_com[self.iter_number, :] = self.Btt.dot(modes)

        ###########################################################################
        ## Bandwidth error
        ###########################################################################
        C = self.mod_com[self.iter_number, :] - self.mod_com[self.iter_number - 1, :]

        self.bp_com[self.iter_number, :] = self.gRD.dot(
                self.bp_com[self.iter_number - 1, :]) - C

        ###########################################################################
        ## Tomographic error
        ###########################################################################
        #G = F - (mod_com[self.iter_number,:] + Ageom - np.dot(RDgeom,com[self.iter_number-1,:]))
        for w in range(len(self.config.p_wfss)):
            self.wfs.raytrace(w, b"atmos", self.tel, self.atm, self.dms)
        self.rtc.do_control_geo_on_wfs(1, self.dms, self.wfs, 0)
        G = self.rtc.get_com(1)
        modes = self.P.dot(G)
        modes[-self.nfiltered - 2:-2] = 0
        self.wf_com[self.iter_number, :] = self.Btt.dot(modes)

        G = self.mod_com[self.iter_number, :] - self.wf_com[self.iter_number, :]
        if (self.iter_number + 1 < self.config.p_loop.niter):
            self.tomo_com[self.iter_number + 1, :] = self.gRD.dot(
                    self.tomo_com[self.iter_number, :]) - g * self.gamma * self.RD.dot(G)

        # Without anyone noticing...
        self.tar.set_phase(0, tarphase)
        self.rtc.set_com(0, Dcom)

    def save_in_hdf5(self, savename):
        tmp = (self.config.p_geom._ipupil.shape[0] -
               (self.config.p_dms[0]._n2 - self.config.p_dms[0]._n1 + 1)) // 2
        tmp_e0 = self.config.p_geom._ipupil.shape[0] - tmp
        tmp_e1 = self.config.p_geom._ipupil.shape[1] - tmp
        pup = self.config.p_geom._ipupil[tmp:tmp_e0, tmp:tmp_e1]
        indx_pup = np.where(pup.flatten() > 0)[0].astype(np.int32)
        dm_dim = self.config.p_dms[0]._n2 - self.config.p_dms[0]._n1 + 1
        self.cov_cor()
        psf = self.tar.get_image(0, b"le")

        fname = "/home/fferreira/Data/" + savename
        pdict = {
                "noise":
                        self.noise_com[self.N_preloop:, :].T, "aliasing":
                                self.alias_wfs_com[self.N_preloop:, :].T, "tomography":
                                        self.tomo_com[self.N_preloop:, :].T,
                "filtered modes":
                        self.H_com[self.N_preloop:, :].T, "non linearity":
                                self.trunc_com[self.N_preloop:, :].T, "bandwidth":
                                        self.bp_com[self.N_preloop:, :].T, "wf_com":
                                                self.wf_com[self.N_preloop:, :].T, "P":
                                                        self.P, "Btt":
                                                                self.Btt, "IF.data":
                                                                        self.IFpzt.data,
                "IF.indices":
                        self.IFpzt.indices, "IF.indptr":
                                self.IFpzt.indptr, "TT":
                                        self.TT, "dm_dim":
                                                dm_dim, "indx_pup":
                                                        indx_pup,
                "fitting":
                        self.fit[self.N_preloop:], "SR":
                                self.SR, "SR2":
                                        self.SR2, "cov":
                                                self.cov, "cor":
                                                        self.cor,
                "psfortho":
                        np.fft.fftshift(self.psf_ortho) /
                        (self.config.p_loop.niter - self.N_preloop), "E":
                                self.Ee, "F":
                                        self.Ff, "dm.xpos":
                                                self.config.p_dms[0]._xpos, "dm.ypos":
                                                        self.config.p_dms[0]._ypos, "R":
                                                                self.rtc.get_cmat(0),
                "D":
                        self.rtc.get_imat(0), "Nact":
                                self.Nact
        }
        h5u.save_h5(fname, "psf", self.config, psf)
        for k in list(pdict.keys()):
            h5u.save_hdf5(fname, k, pdict[k])

    def cov_cor(self):
        self.cov = np.zeros((6, 6))
        self.cor = np.zeros((6, 6))
        bufdict = {
                "0": self.noise_com.T, "1": self.trunc_com.T, "2": self.alias_wfs_com.T,
                "3": self.H_com.T, "4": self.bp_com.T, "5": self.tomo_com.T
        }
        for i in range(self.cov.shape[0]):
            for j in range(self.cov.shape[1]):
                if (j >= i):
                    tmpi = self.P.dot(bufdict[str(i)])
                    tmpj = self.P.dot(bufdict[str(j)])
                    self.cov[i, j] = np.sum(
                            np.mean(tmpi * tmpj, axis=1) -
                            np.mean(tmpi, axis=1) * np.mean(tmpj, axis=1))
                else:
                    self.cov[i, j] = self.cov[j, i]

        s = np.reshape(np.diag(self.cov), (self.cov.shape[0], 1))
        sst = np.dot(s, s.T)
        ok = np.where(sst)
        self.cor[ok] = self.cov[ok] / np.sqrt(sst[ok])


###############################################################################
#
#                                 MAIN
#
###############################################################################
if __name__ == "__main__":
    if (len(sys.argv) < 2):
        error = 'command line should be at least:"python -i test.py parameters_filename"\n with "parameters_filename" the path to the parameters file'
        raise Exception(error)

    #get parameters from file
    param_file = sys.argv[1]

    if (len(sys.argv) > 2):
        savename = sys.argv[2]
    else:
        savename = "roket_default.h5"

    roket = Roket(param_file)
    roket.init_sim()
    roket.loop()
    roket.save_in_hdf5(savename)
