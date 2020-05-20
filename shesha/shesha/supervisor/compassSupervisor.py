## @package   shesha.supervisor.compassSupervisor
## @brief     Initialization and execution of a COMPASS supervisor
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   5.0.0
## @date      2020/05/18
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
#  All rights reserved.
#  Distributed under GNU - LGPL
#
#  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
#  General Public License as published by the Free Software Foundation, either version 3 of the License,
#  or any later version.
#
#  COMPASS: End-to-end AO simulation tool using GPU acceleration
#  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
#
#  The final product includes a software package for simulating all the critical subcomponents of AO,
#  particularly in the context of the ELT and a real-time core based on several control approaches,
#  with performances consistent with its integration into an instrument. Taking advantage of the specific
#  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
#  conduct large simulation campaigns called to the ELT.
#
#  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
#  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
#  various systems configurations such as multi-conjugate AO.
#
#  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
#  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
#  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.

from shesha.supervisor.aoSupervisor import AoSupervisor
import numpy as np

import shesha.constants as scons
from shesha.constants import CONST
import shesha.ao.basis as basis
import astropy.io.fits as pfits
from tqdm import trange, tqdm
import time

from typing import List


class CompassSupervisor(AoSupervisor):

    def __init__(self, config_file: str = None, cacao: bool = False,
                 use_DB: bool = False):
        """ Init the COMPASS supervisor

        Parameters:
            config_file: (str): (optionnal) Path to the parameter file

            cacao: (bool): (optionnal) Flag to enable cacao

            use_DB: (bool): (optionnal) Flag to enable database
        """
        self._sim = None
        self.config = None
        self.cacao = cacao
        self.use_DB = use_DB
        self.P = None
        self.modal_basis = None
        if config_file is not None:
            self.load_config(config_file=config_file)

    def __repr__(self):
        return object.__repr__(self) + str(self._sim)

    #     _    _         _                  _
    #    / \  | |__  ___| |_ _ __ __ _  ___| |_
    #   / _ \ | '_ \/ __| __| '__/ _` |/ __| __|
    #  / ___ \| |_) \__ \ |_| | | (_| | (__| |_
    # /_/   \_\_.__/|___/\__|_|  \__,_|\___|\__|
    #
    #  __  __      _   _               _
    # |  \/  | ___| |_| |__   ___   __| |___
    # | |\/| |/ _ \ __| '_ \ / _ \ / _` / __|
    # | |  | |  __/ |_| | | | (_) | (_| \__ \
    # |_|  |_|\___|\__|_| |_|\___/ \__,_|___/

    def single_next(self, move_atmos: bool = True, show_atmos: bool = True,
                   get_tar_image: bool = False, get_residual: bool = False) -> None:
        """ Performs a single loop iteration 

        Parameters:
            move_atmos : (bool, optional) : Move the atmosphere layers. Default is True

            show_atmos : (bool, optional) : WFS and targets see the atmosphere layers. Default is True

            get_tar_image : (bool, optional) : 

            get_residual : (bool, optional) :
        """
        self._sim.next(see_atmos=show_atmos)  # why not self.atmos.is_enable?


    def set_command(self, nctrl: int, command: np.ndarray) -> None:
        """ Set the RTC command vector. It is set as it has just been computed
        from the do_control() method. Then, this is not the vector that will
        be send to the DM directly (latency will be computed, clipping, etc...)

        See set_voltages() to set directly the DM vector

        Parameters:
            nctrl : (int) : controller index

            command : (np.ndarray) : command vector
        """
        self._sim.rtc.d_control[nctrl].set_com(command, command.size)

    def get_command(self, nctrl: int):
        """ Get the command vector computed after do_control() method from ncontrol controller.
        This is not the vector that is applied on the DM

        See get_voltages() to get the vector applied on the DM

        Parameters:
            nctrl : (int) : controller index

        Return:
            com : (np.ndarray) : current command vector
        """
        return np.array(self._sim.rtc.d_control[ncontrol].d_com)

    #  ____                  _ _   _        __  __      _   _               _
    # / ___| _ __   ___  ___(_) |_(_) ___  |  \/  | ___| |_| |__   ___   __| |___
    # \___ \| '_ \ / _ \/ __| | __| |/ __| | |\/| |/ _ \ __| '_ \ / _ \ / _` / __|
    #  ___) | |_) |  __/ (__| | |_| | (__  | |  | |  __/ |_| | | | (_) | (_| \__ \
    # |____/| .__/ \___|\___|_|\__|_|\___| |_|  |_|\___|\__|_| |_|\___/ \__,_|___/
    #       |_|


    def set_pyr_modulation_ampli(self, wfs_index: int, pyr_mod: float) -> None:
        """ Set pyramid circular modulation amplitude value - in lambda/D units.

        Compute new modulation points corresponding to the new amplitude value
        and upload them

        Parameters:
            wfs_index : (int) : WFS index

            pyr_mod : (float) : new pyramid modulation amplitude value
        """
        from shesha.ao.wfs import comp_new_pyr_ampl
        p_wfs = self._sim.config.p_wfss[wfs_index]

        cx, cy, scale, pyr_npts = comp_new_pyr_ampl(wfs_index, pyr_mod,
                                                    self._sim.config.p_wfss,
                                                    self._sim.config.p_tel)
        p_wfs.set_pyr_ampli(pyr_mod)
        self.set_pyr_modulation_points(wfs_index, cx, cy)
        self._sim.rtc.d_centro[wfs_index].set_scale(scale)

        if (len(p_wfs._halfxy.shape) == 2):
            print("PYR modulation set to: %f L/D using %d points" % (pyr_mod, pyr_npts))
        elif (len(p_wfs._halfxy.shape) == 3):
            newhalfxy = np.tile(p_wfs._halfxy[0, :, :], (pyr_npts, 1, 1))
            print("Loading new modulation arrays")
            self._sim.wfs.d_wfs[wfs_index].set_phalfxy(
                    np.exp(1j * newhalfxy).astype(np.complex64).T)
            print("Done. PYR modulation set to: %f L/D using %d points" % (pyr_mod,
                                                                           pyr_npts))
        else:
            raise ValueError("Error unknown p_wfs._halfxy shape")
        self._sim.rtc.do_centroids(0)  # To be ready for the next get_slopess




    def loop(self, n: int = 1, monitoring_freq: int = 100, **kwargs):
        """ Perform the AO loop for n iterations

        Parameters:
            n: (int, optional) : Number of iteration that will be done

            monitoring_freq: (int, optional) : Monitoring frequency [frames]
        """
        self._sim.loop(n, monitoring_freq=monitoring_freq, **kwargs)

    def force_context(self) -> None:
        """ Active all the GPU devices specified in the parameters file
        """
        self._sim.force_context()





    def reset_simu(self):
        """ Reset the simulation to return to its original state
        """
        self.atmos.reset_turbu()
        self.wfs.reset_noise()
        for tar_index in range(self._sim.tar.ntargets):
            self.target.reset_strehl(tar_index)
        self.dms.reset_dm()
        self.rtc.open_loop()
        self.rtc.close_loop()




    def load_config(self, config_file: str = None, sim=None) -> None:
        """ Init the COMPASS simulator wih the config_file

        Parameters:
            config_file : (str, optional) : path to the configuration file

            sim : (Simulator, optional) : Simulator instance
        """
        if self._sim is None:
            if sim is None:
                if self.cacao:
                    from shesha.sim.simulatorCacao import SimulatorCacao as Simulator
                else:
                    from shesha.sim.simulator import Simulator
                self._sim = Simulator(filepath=config_file, use_DB=self.use_DB)
            else:
                self._sim = sim
        else:
            self._sim.clear_init()
            self._sim.load_from_file(config_file)
        self.config = self._sim.config

    def is_init(self) -> bool:
        """ Return the status of COMPASS init

        Return:
            is_init : (bool) : Status of the initialisation
        """
        return self._sim.is_init

    def clear_init_sim(self) -> None:
        """ Clear the initialization of the simulation
        """
        self._sim.clear_init()

    def init_config(self) -> None:
        """ Initialize the simulation
        """
        self._sim.init_sim()
        self.rtc = self._sim.rtc
        self.iter = self._sim.iter
        self.enable_atmos(True)
        self.is_init = True









    def get_influ_basis_sparse(self, ncontrol: int):
        """ Return the influence function basis of all the DM as a sparse matrix
        Controller GEO only

        Parameters:
            ncontrol : (int) : Index of the controller GEO

        Return:
            influ_sparse : (scipy csr_matrix) : influence function phases
        """
        return self._sim.rtc.d_control[ncontrol].d_IFsparse.get_csr()

    def get_tt_influ_basis(self, ncontrol: int) -> np.ndarray:
        """ Return the influence function basisof a tip-tilt mirror
        Controller GEO only

        Parameters:
            ncontrol : (int) : Index of the controller GEO

        Return:
            influ : (np.ndarray) : influence function phases
        """
        return np.array(self._sim.rtc.d_control[ncontrol].d_TT)

    def compute_influ_basis(self, dm_index: int):
        """ Computes and return the influence function phase basis of the specified DM
        as a sparse matrix

        Parameters:
            dm_index : (int) : Index of the DM

        Return:
            influ_sparse : (scipy csr_matrix) : influence function phases
        """
        from shesha.ao import basis
        influ_sparse = basis.compute_dm_basis(self._sim.dms.d_dms[dm_index],
                                              self._sim.config.p_dms[dm_index],
                                              self._sim.config.p_geom)

        return influ_sparse





    def get_pupil(self, pupil_type: str = "ipupil") -> np.ndarray:
        """ Return the selected pupil.
        Available pupils are :
            spupil : smallest one (use for target phase screen)
            mpupil : medium one (use for wfs phase screen)
            ipupil : biggest one (use for FFT)

        Parameters:
            pupil_type : (str, optional) : Pupil to return ("spupil","mpupil" or "ipupil")

        Return:
            pupil : (np.ndarray) : pupil
        """
        if pupil_type == "spupil":
            pupil = self._sim.config.p_geom._spupil
        elif pupil_type == "mpupil":
            pupil = self._sim.config.p_geom._mpupil
        elif pupil_type == "ipupil":
            pupil = self._sim.config.p_geom._ipupil
        else:
            raise ArgumentError("Unknown pupil")

        return pupil



    def get_selected_pix(self) -> np.ndarray:
        """ Return the pyramid image with only the selected pixels used by the full pixels centroider

        Return:
            selected_pix : (np.ndarray) : PWFS image with only selected pixels
        """
        if (self.config.p_centroiders[0].type != scons.CentroiderType.MASKEDPIX):
            raise TypeError("Centroider must be maskedPix")

        carma_centroids = self._sim.rtc.d_control[0].d_centroids
        self._sim.rtc.d_centro[0].fill_selected_pix(carma_centroids)

        return np.array(self._sim.rtc.d_centro[0].d_selected_pix)

    def first_non_zero(self, array: np.ndarray, axis: int,
                       invalid_val: int = -1) -> np.ndarray:
        """ Find the first non zero element of an array

        Parameters:
            array : (np.ndarray) : input array

            axis : (int) : axis index

            invalid_val : (int, optional) : Default is -1

        Return:
            non_zeros_pos : (np.ndarray) : Index of the first non-zero element
                                            for each line or column following the axis
        """
        mask = array != 0
        return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)

    def compute_modes_to_volts_basis(self, modal_basis_type: str, merged: bool = False,
                                     nbpairs: int = None, return_delta: bool = False):
        """ Computes a given modal basis ("KL2V", "Btt", "Btt_petal") and return the 2 transfer matrices

        Parameters:
            modal_basis_type : (str) : modal basis to compute ("KL2V", "Btt", "Btt_petal")

            merged : (bool, optional) :

            nbpairs : (int, optional) :

        Return:
            modal_basis : (np.ndarray) : modes to volts matrix

            P : (np.ndarray) : volts to modes matrix (None if "KL")
        """
        if (modal_basis_type == "KL2V"):
            print("Computing KL2V basis...")
            self.modal_basis = basis.compute_KL2V(
                    self._sim.config.p_controllers[0], self._sim.dms,
                    self._sim.config.p_dms, self._sim.config.p_geom,
                    self._sim.config.p_atmos, self._sim.config.p_tel)
            fnz = self.first_non_zero(self.modal_basis, axis=0)
            # Computing the sign of the first non zero element
            #sig = np.sign(self.modal_basis[[fnz, np.arange(self.modal_basis.shape[1])]])
            sig = np.sign(self.modal_basis[tuple([
                    fnz, np.arange(self.modal_basis.shape[1])
            ])])  # pour remove le future warning!
            self.modal_basis *= sig[None, :]
            self.P = None
        elif (modal_basis_type == "Btt"):
            print("Computing Btt basis...")
            self.modal_basis, self.P = self.compute_btt(inv_method="cpu_svd",
                                                        merged=merged, nbpairs=nbpairs,
                                                        return_delta=return_delta)
            fnz = self.first_non_zero(self.modal_basis, axis=0)
            # Computing the sign of the first non zero element
            #sig = np.sign(self.modal_basis[[fnz, np.arange(self.modal_basis.shape[1])]])
            sig = np.sign(self.modal_basis[tuple([
                    fnz, np.arange(self.modal_basis.shape[1])
            ])])  # pour remove le future warning!
            self.modal_basis *= sig[None, :]
        elif (modal_basis_type == "Btt_petal"):
            print("Computing Btt with a Petal basis...")
            self.modal_basis, self.P = self.compute_btt_petal()
        else:
            raise ArgumentError("Unsupported modal basis")

        return self.modal_basis, self.P

    def compute_btt_basis(self, merged: bool = False, nbpairs: int = None,
                          return_delta: bool = False):
        """ Computes the so-called Btt modal basis. The <merged> flag allows merto merge
        2x2 the actuators influence functions for actuators on each side of the spider (ELT case)

        Parameters:
            merged : (bool, optional) : If True, merge 2x2 the actuators influence functions for
                                        actuators on each side of the spider (ELT case). Default
                                        is False

            nbpairs : (int, optional) : Default is None. TODO : description

            return_delta : (bool, optional) : If False (default), the function returns
                                              Btt (modes to volts matrix),
                                              and P (volts to mode matrix).
                                              If True, returns delta = IF.T.dot(IF) / N
                                              instead of P

        Return:
            Btt : (np.ndarray) : Btt modes to volts matrix

            P : (np.ndarray) : volts to Btt modes matrix
        """
        influ_basis = self.get_influ_basis_sparse(1)
        tt_basis = self.get_tt_influ_basis(1)
        if (merged):
            couples_actus, index_under_spiders = self.compute_merged_influ(
                    nbpairs=nbpairs)
            influ_basis2 = influ_basis.copy()
            index_remove = index_under_spiders.copy()
            index_remove += list(couples_actus[:, 1])
            print("Pairing Actuators...")
            for i in tqdm(range(couples_actus.shape[0])):
                influ_basis2[couples_actus[i, 0], :] += influ_basis2[
                        couples_actus[i, 1], :]
            print("Pairing Done")
            boolarray = np.zeros(influ_basis2.shape[0], dtype=np.bool)
            boolarray[index_remove] = True
            self.slaved_actus = boolarray
            self.selected_actus = ~boolarray
            self.couples_actus = couples_actus
            self.index_under_spiders = index_under_spiders
            influ_basis2 = influ_basis2[~boolarray, :]
            influ_basis = influ_basis2
        else:
            self.slaved_actus = None
            self.selected_actus = None
            self.couples_actus = None
            self.index_under_spiders = None

        Btt, P = basis.compute_btt(influ_basis.T, tt_basis, return_delta=return_delta)

        if (merged):
            Btt2 = np.zeros((len(boolarray) + 2, Btt.shape[1]))
            Btt2[np.r_[~boolarray, True, True], :] = Btt
            Btt2[couples_actus[:, 1], :] = Btt2[couples_actus[:, 0], :]

            P2 = np.zeros((Btt.shape[1], len(boolarray) + 2))
            P2[:, np.r_[~boolarray, True, True]] = P
            P2[:, couples_actus[:, 1]] = P2[:, couples_actus[:, 0]]
            Btt = Btt2
            P = P2
        if (return_delta):
            P = delta
        return Btt, P

    def compute_merged_influ(self, nbpairs: int = None):
        """ Used to compute merged IF from each side of the spider
        for an ELT case (Petalling Effect)

        Parameters:
            nbpairs : (int, optional) : Default is None. TODO : description

        Return:
            pairs : (np.ndarray) : TODO description

            discard : (list) : TODO description
        """
        p_geom = self._sim.config.p_geom
        import shesha.util.make_pupil as mkP
        import shesha.util.utilities as util
        import scipy.ndimage

        cent = p_geom.pupdiam / 2. + 0.5
        p_tel = self._sim.config.p_tel
        p_tel.t_spiders = 0.51
        spup = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent,
                              cent).astype(np.float32).T

        p_tel.t_spiders = 0.
        spup2 = mkP.make_pupil(p_geom.pupdiam, p_geom.pupdiam, p_tel, cent,
                               cent).astype(np.float32).T

        spiders = spup2 - spup

        (spidersID, k) = scipy.ndimage.label(spiders)
        spidersi = util.pad_array(spidersID, p_geom.ssize).astype(np.float32)
        px_list_spider = [np.where(spidersi == i) for i in range(1, k + 1)]

        # DM positions in iPupil:
        dm_posx = self._sim.config.p_dm0._xpos - 0.5
        dm_posy = self._sim.config.p_dm0._ypos - 0.5
        dm_pos_mat = np.c_[dm_posx, dm_posy].T  # one actu per column

        pitch = self._sim.config.p_dm0._pitch
        discard = np.zeros(len(dm_posx), dtype=np.bool)
        pairs = []

        # For each of the k pieces of the spider
        for k, px_list in enumerate(px_list_spider):
            pts = np.c_[px_list[1],
                        px_list[0]]  # x,y coord of pixels of the spider piece
            # line_eq = [a, b]
            # Which minimizes leqst squares of aa*x + bb*y = 1
            line_eq = np.linalg.pinv(pts).dot(np.ones(pts.shape[0]))
            aa, bb = line_eq[0], line_eq[1]

            # Find any point of the fitted line.
            # For simplicity, the intercept with one of the axes x = 0 / y = 0
            if np.abs(bb) < np.abs(aa):  # near vertical
                one_point = np.array([1 / aa, 0.])
            else:  # otherwise
                one_point = np.array([0., 1 / bb])

            # Rotation that aligns the spider piece to the horizontal
            rotation = np.array([[-bb, aa], [-aa, -bb]]) / (aa**2 + bb**2)**.5

            # Rotated the spider mask
            rotated_px = rotation.dot(pts.T - one_point[:, None])
            # Min and max coordinates along the spider length - to filter actuators that are on
            # 'This' side of the pupil and not the other side
            min_u, max_u = rotated_px[0].min() - 5. * pitch, rotated_px[0].max(
            ) + 5. * pitch

            # Rotate the actuators
            rotated_actus = rotation.dot(dm_pos_mat - one_point[:, None])
            sel_good_side = (rotated_actus[0] > min_u) & (rotated_actus[0] < max_u)
            threshold = 0.05
            # Actuators below this piece of spider
            sel_discard = (np.abs(rotated_actus[1]) < threshold * pitch) & sel_good_side
            discard |= sel_discard

            # Actuator 'near' this piece of spider
            sel_pairable = (np.abs(rotated_actus[1]) > threshold  * pitch) & \
                            (np.abs(rotated_actus[1]) < 1. * pitch) & \
                            sel_good_side

            pairable_index = np.where(sel_pairable)[0]  # Indices of these actuators
            u_coord = rotated_actus[
                    0, sel_pairable]  # Their linear coord along the spider major axis

            order = np.sort(u_coord)  # Sort by linear coordinate
            order_index = pairable_index[np.argsort(
                    u_coord)]  # And keep track of original indexes

            # i = 0
            # while i < len(order) - 1:
            if (nbpairs is None):
                i = 0
                ii = len(order) - 1
            else:
                i = len(order) // 2 - nbpairs
                ii = len(order) // 2 + nbpairs
            while (i < ii):
                # Check if next actu in sorted order is very close
                # Some lonely actuators may be hanging in this list
                if np.abs(order[i] - order[i + 1]) < .2 * pitch:
                    pairs += [(order_index[i], order_indx[i + 1])]
                    i += 2
                else:
                    i += 1
        print('To discard: %u actu' % np.sum(discard))
        print('%u pairs to slave' % len(pairs))
        if np.sum(discard) == 0:
            discard = []
        else:
            list(np.where(discard)[0])
        return np.asarray(pairs), list(np.where(discard)[0])

    def compute_btt_petal(self):
        """ Computes a Btt modal basis with Pistons filtered

        Return:
            Btt : (np.ndarray) : Btt modes to volts matrix

            P : (np.ndarray) : volts to Btt modes matrix
        """
        influ_pzt = self.get_influ_basis_sparse(1)
        petal_dm_index = np.where([
                d.influ_type is scons.InfluType.PETAL for d in self.config.p_dms
        ])[0][0]
        influ_petal = self.compute_influ_basis(petal_dm_index)
        tt_index = np.where([d.type is scons.DmType.TT for d in self.config.p_dms])[0][0]

        influ_tt = self.compute_influ_basis(tt_index)

        return basis.compute_btt(influ_pzt.T, influ_tt, influ_petal=influ_petal)

    def compute_phase_to_modes(self, modal_basis: np.ndarray):
        """ Return the phase to modes matrix by using the given modal basis

        Parameters:
            modal_basis : (np.ndarray) : Modal basis matrix

        Return:
            phase_to_modes : (np.ndarray) : phase to modes matrix
        """
        old_noise = self.config.p_wfss[0].noise
        self.set_noise(0, -1)

        nbmode = modal_basis.shape[1]
        phase = self.get_tar_phase(0)
        phase_to_modes = np.zeros((nbmode, phase.shape[0], phase.shape[1]))
        S = np.sum(pup)
        for i in range(nbmode):
            self.target.reset_tar_phase(0)
            self.set_command((modal_basis[:, i]).copy())
            self._sim.next(see_atmos=False)
            phase = self.get_tar_phase(0, pupil=true)
            # Normalisation pour les unites rms en microns !!!
            norm = np.sqrt(np.sum((phase)**2) / S)
            phase_to_modes[i] = phase / norm
        self.phase_to_modes = phase_to_modes
        self.set_noise(0, old_noise)
        return phase_to_modes


    def do_imat_modal(self, controller_index, ampli, modal_basis, noise=False,
                      nmodes_max=0, with_turbu=False, push_pull=False):
        """ Computes an interaction matrix from provided modal basis

        Parameters:
            controller_index : (int) : Controller index

            ampli : (np.ndarray) : amplitude to apply on each mode

            modal_basis : (np.ndarray) : modal basis matrix

            noise : (bool, optional) : Flag to enable noise for WFS image compuation. Default is False

            nmodes_max : (int, optional) : Default is 0. TODO : description

            with_turbu : (bool, optional) : Flag to enable atmosphere for WFS phase screen raytracing.
                                            Default is False

            push_pull : (bool, optional) : If True, imat is computed as an average of push and pull ampli
                                            on each mode

        Return:
            modal_imat : (np.ndarray) : Modal interaction matrix
        """
        nslopes = self.config.p_controllers[controller_index].nslope
        modal_imat = np.zeros((nslopes, modal_basis.shape[1]))

        if (nmodes_max == 0):
            nmodes_max = modal_basis.shape[1]
        v_old = self.get_command(controller_index)
        self.rtc.open_loop(rst=False)
        for m in range(nmodes_max):
            # v = ampli[m] * modal_basis[:, m:m + 1].T.copy()
            v = ampli[m] * modal_basis[:, m]
            if ((push_pull is True) or
                (with_turbu is True)):  # with turbulence/aberrations => push/pull
                self.rtc.set_perturbation_voltage(
                        0, "imat_modal",
                        v_old + v)  # Adding Perturbation voltage on current iteration
                devpos = self.apply_volts_and_get_slopes(controller_index,
                                                         turbu=with_turbu, noise=noise)
                self.rtc.set_perturbation_voltage(controller_index, "imat_modal", v_old - v)
                devmin = self.apply_volts_and_get_slopes(controller_index,
                                                         turbu=with_turbu, noise=noise)
                modal_imat[:, m] = (devpos - devmin) / (2. * ampli[m])
                #imat[:-2, :] /= pushDMMic
                #if(nmodes_max == 0):# i.e we measured all modes including TT
                #imat[-2:, :] /= pushTTArcsec
            else:  # No turbulence => push only
                self.rtc.open_loop()  # open_loop
                self.rtc.set_perturbation_voltage(controller_index, "imat_modal", v)
                modal_imat[:, m] = self.apply_volts_and_get_slopes(controller_index, noise=noise) / ampli[m]
        self.remove_perturbation_voltage(controller_index, "imat_modal")
        if ((push_pull is True) or (with_turbu is True)):
            self.rtc.close_loop()  # We are supposed to be in close loop now
        return modal_imat

    def do_imat_phase(self, controller_index: int, cube_phase: np.ndarray, noise=False,
                      nmodes_max=0, with_turbu=False, push_pull=False, wfs_index=0):
        """ Computes an interaction matrix with the provided cube phase

        Parameters:
            controller_index : (int) : Controller index

            cube_phase : (np.ndarray) : Cube of phase to insert as NCPA

            noise : (bool, optional) : Flag to enable noise for WFS image compuation. Default is False

            nmodes_max : (int, optional) : Default is 0. TODO : description

            with_turbu : (bool, optional) : Flag to enable atmosphere for WFS phase screen raytracing.
                                            Default is False

            push_pull : (bool, optional) : If True, imat is computed as an average of push and pull ampli
                                            on each mode

            wfs_index : (int, optional) : WFS index. Default is 0

        Return:
            phase_imat : (np.ndarray) : Phase interaction matrix
        """
        nslopes = self.config.p_controllers[controller_index].nslope
        imat_phase = np.zeros((cube_phase.shape[0], nslopes))
        for nphase in range(cube_phase.shape[0]):
            if ((push_pull is True) or (with_turbu is True)
                ):  # with turbulence/aberrations => push/pullADOPT/projects/cosmic/
                self.set_ncpa_wfs(wfs_index, cube_phase[nphase, :, :])
                devpos = self.apply_volts_and_get_slopes(controller_index,
                                                         turbu=with_turbu, noise=noise)
                self.set_ncpa_wfs(wfs_index, -cube_phase[nphase, :, :])
                devmin = self.apply_volts_and_get_slopes(controller_index,
                                                         turbu=with_turbu, noise=noise)
                imat_phase[nphase, :] = (devpos - devmin) / 2
            else:  # No turbulence => push only
                self.rtc.open_loop()  # open_loop
                self.set_ncpa_wfs(wfs_index, cube_phase[nphase, :, :])
                imat_phase[nphase, :] = self.apply_volts_and_get_slopes(
                        controller_index, noise=noise)
        self.set_ncpa_wfs(wfs_index,
                          cube_phase[nphase, :, :] * 0.)  # Remove the Phase on WFS
        _ = self.apply_volts_and_get_slopes(controller_index, turbu=with_turbu,
                                            noise=noise)

        return imat_phase

    def compute_modal_residuals(self):
        """ Computes the modal residual coefficients of the residual phase.

        Uses the P matrix computed from compute_modes_to_volts_basis

        Return:
            ai : (np.ndarray) : Modal coefficients
        """
        try:
            self._sim.do_control(1, 0)
        except:
            return [0]
        v = self.get_command(
                1
        )  # We compute here the residual phase on the DM modes. Gives the Equivalent volts to apply/
        if (self.P is None):
            return [0]
            # self.modal_basis, self.P = self.compute_modes_to_volts_basis("Btt")
        if (self.selected_actus is None):
            ai = self.P.dot(v) * 1000.  # np rms units
        else:  # Slaving actus case
            v2 = v[:-2][list(
                    self.selected_actus)]  # If actus are slaved then we select them.
            v3 = v[-2:]
            ai = self.P.dot(np.concatenate((v2, v3))) * 1000.
        return ai

    def set_pyr_multiple_stars_source(self, wfs_index: int, coords: List,
                                      weights: List = None, pyr_mod: float = 3.,
                                      niters: int = None):
        """ Sets the Pyramid modulation points with a multiple star system

        Parameters:
            wfs_index : (int) : WFS index

            coords : (list) : list of couples of length n, coordinates of the n stars in lambda/D

            weights : (list, optional) : list of weights to apply on each modulation points. Default is None

            pyr_mod : (float, optional): modulation amplitude of the pyramid in lambda/D. Default is 3

            niters : (int, optional) : number of iteration. Default is None
        """
        if niters is None:
            perim = pyr_mod * 2 * np.pi
            niters = int((perim // 4 + 1) * 4)
            print(niters)
        scale_circ = self.config.p_wfss[wfs_index]._pyr_scale_pos * pyr_mod
        temp_cx = []
        temp_cy = []
        for k in coords:
            temp_cx.append(scale_circ * \
                np.sin((np.arange(niters)) * 2. * np.pi / niters) + \
                k[0] * self.config.p_wfss[wfs_index]._pyr_scale_pos)
            temp_cy.append(scale_circ * \
                np.cos((np.arange(niters)) * 2. * np.pi / niters) + \
                k[1] * self.config.p_wfss[wfs_index]._pyr_scale_pos)
        cx = np.concatenate(np.array(temp_cx))
        cy = np.concatenate(np.array(temp_cy))
        #Gives the arguments to the simulation
        if weights is not None:
            w = []
            for k in weights:
                w += niters * [k]
            weights = np.array(w)
        self.set_pyr_modulation_points(wfs_index, cx, cy, weights)

    def set_pyr_disk_source_hexa(self, wfs_index: int, radius: float):
        """ Create disk object by packing PSF in a given radius, using hexagonal packing

        /!\ There is no modulation

        Parameters:
            wfs_index  : (int) : WFS index

            radius : (float) : radius of the disk object in lambda/D
        """
        #Vectors used to generate the hexagonal paving
        gen_xp, gen_yp = np.array([1,
                                   0.]), np.array([np.cos(np.pi / 3),
                                                   np.sin(np.pi / 3)])
        n = 1 + int(1.2 * radius)
        mat_circ = []
        for k in range(-n, n):
            for l in range(-n, n):
                coord = k * gen_xp + l * gen_yp
                if np.sqrt(coord[0]**2 + coord[1]**2) <= radius:
                    mat_circ.append(coord)
        mat_circ = np.array(mat_circ)
        cx, cy = mat_circ[:, 0], mat_circ[:, 1]
        self.set_pyr_modulation_points(wfs_index, cx, cy)

    def generate_square(self, radius: float, density: float = 1.):
        """ Generate modulation points positions following a square pattern

        Parameters:
            radius : (float) : half the length of a side in lambda/D

            density : (float), optional) : number of psf per lambda/D. Default is 1

        Return:
            cx : (np.ndarray) : X-positions of the modulation points

            cy : (np.ndarray) : Y-positions of the modulation points
        """
        x = np.linspace(-radius, radius, 1 + 2 * int(radius * density))
        cx, cy = np.meshgrid(x, x, indexing='ij')
        cx = cx.flatten()
        cy = cy.flatten()
        return (cx, cy)

    def generate_circle(self, radius: float, density: float = 1.):
        """ Generate modulation points positions following a circular pattern

        Parameters:
            radius : (float) : half the length of a side in lambda/D

            density : (float), optional) : number of psf per lambda/D. Default is 1

        Return:
            cx : (np.ndarray) : X-positions of the modulation points

            cy : (np.ndarray) : Y-positions of the modulation points
        """
        cx, cy = self.generate_square(radius, density)
        r = cx * cx + cy * cy <= radius**2
        return (cx[r], cy[r])

    def set_pyr_disk_source(self, wfs_index: int, radius: float, density: float = 1.):
        """ Create disk object by packing PSF in a given radius, using square packing

        /!\ There is no modulation

        Parameters:
            wfs_index  : (int) : WFS index

            radius : (float) : radius of the disk object in lambda/D

            density : (float, optional) : Spacing between the packed PSF in the disk object, in lambda/D.
                                          Default is 1
        """
        cx, cy = self.generate_circle(radius, density)
        cx = cx.flatten() * self.config.p_wfss[wfs_index]._pyr_scale_pos
        cy = cy.flatten() * self.config.p_wfss[wfs_index]._pyr_scale_pos
        self.set_pyr_modulation_points(wfs_index, cx, cy)

    def set_pyr_square_source(self, wfs_index: int, radius: float, density: float = 1.):
        """ Create a square object by packing PSF in a given radius, using square packing

        /!\ There is no modulation

        Parameters:
            wfs_index  : (int) : WFS index

            radius : (float) : radius of the disk object in lambda/D

            density : (float, optional) : Spacing between the packed PSF in the disk object, in lambda/D.
                                          Default is 1
        """
        cx, cy = self.generate_square(radius, density)
        cx = cx.flatten() * self.config.p_wfss[wfs_index]._pyr_scale_pos
        cy = cy.flatten() * self.config.p_wfss[wfs_index]._pyr_scale_pos
        self.set_pyr_modulation_points(wfs_index, cx, cy)

    def generate_pseudo_source(self, radius: float, additional_psf=0, density=1.):
        """ Used to generate a pseudo source for PYRWFS

        Parameters:
            radius : (float) : TODO description

            additional_psf : (int) :TODO description

            density : (float, optional) :TODO description

        Return:
            ox : TODO description & explicit naming

            oy : TODO description & explicit naming

            w : TODO description & explicit naming

            xc : TODO description & explicit naming

            yc : TODO description & explicit naming
        """
        struct_size = (1 + 2 * additional_psf)**2
        center_x, center_y = self.generate_square(additional_psf, density)
        center_weight = (1 + 2 * int(additional_psf * density))**2 * [1]
        center_size = 1 + 2 * int(additional_psf * density)

        weight_edge = [(1 + 2 * int(radius * density) - center_size) // 2]
        xc, yc = self.generate_circle(radius, density)
        for k in range(additional_psf):
            line_length = np.sum(yc == (k + 1))
            print(line_length)
            weight_edge.append((line_length - center_size) // 2)

        edge_dist = (radius + additional_psf) // 2
        V_edge_x = []
        V_edge_y = []
        V_edge_weight = []
        for m in [-1, 1]:
            V_edge_x.append(0)
            V_edge_y.append(m * edge_dist)
            V_edge_weight.append(weight_edge[0])
        for k, val in enumerate(weight_edge[1:]):
            for l in [-1, 1]:
                for m in [-1, 1]:
                    V_edge_x.append(l * (k + 1) * density)
                    V_edge_y.append(m * edge_dist)
                    V_edge_weight.append(val)
        H_edge_x = []
        H_edge_y = []
        H_edge_weight = []
        for m in [-1, 1]:
            H_edge_x.append(m * edge_dist)
            H_edge_y.append(0)
            H_edge_weight.append(weight_edge[0])
        for k, val in enumerate(weight_edge[1:]):
            for l in [-1, 1]:
                for m in [-1, 1]:
                    H_edge_x.append(m * edge_dist)
                    H_edge_y.append(l * (k + 1) * density)
                    H_edge_weight.append(val)
        pup_cent_x = []
        pup_cent_y = []
        pup_cent_weight = 4 * [(len(xc) - 2 * np.sum(H_edge_weight) - struct_size) / 4]
        pup_cent_dist = int(edge_dist // np.sqrt(2))
        for l in [-1, 1]:
            for m in [-1, 1]:
                pup_cent_x.append(l * pup_cent_dist)
                pup_cent_y.append(m * pup_cent_dist)
        ox = np.concatenate((center_x, V_edge_x, H_edge_x, pup_cent_x))
        oy = np.concatenate((center_y, V_edge_y, H_edge_y, pup_cent_y))
        w = np.concatenate((center_weight, V_edge_weight, H_edge_weight,
                            pup_cent_weight))
        return (ox, oy, w, xc, yc)

    def set_pyr_pseudo_source(self, wfs_index: int, radius: float,
                              additional_psf: int = 0, density: float = 1.):
        """ TODO : DESCRIPTION

        Parameters:
            wfs_index : (int) : WFS index

            radius : (float) : TODO : DESCRIPTION

            additional_psf : (int, optional) : TODO : DESCRIPTION

            density : (float, optional) :TODO : DESCRIPTION
        """
        cx, cy, weights, _, _ = self.generate_pseudo_source(radius, additional_psf,
                                                            density)
        cx = cx.flatten() * self.config.p_wfss[wfs_index]._pyr_scale_pos
        cy = cy.flatten() * self.config.p_wfss[wfs_index]._pyr_scale_pos
        self.set_pyr_modulation_points(wfs_index, cx, cy, weights)

    def record_ao_circular_buffer(
            self, cb_count: int, sub_sample: int = 1, controller_index: int = 0,
            tar_index: int = 0, see_atmos: bool = True, cube_data_type: str = None,
            cube_data_file_path: str = "", ncpa: int = 0, ncpa_wfs: np.ndarray = None,
            ref_slopes: np.ndarray = None, ditch_strehl: bool = True):
        """ Used to record a synchronized circular buffer AO loop data.

        Parameters:
            cb_count: (int) : the number of iterations to record.

            sub_sample:  (int) : sub sampling of the data (default=1, I.e no subsampling)

            controller_index:  (int) :

            tar_index:  (int) : target number

            see_atmos:  (int) : used for the next function to enable or not the Atmos

            cube_data_type:   (int) : if  specified ("tarPhase" or "psfse") returns the target phase or short exposure PSF data cube in the output variable

            cube_data_file_path:  (int) : if specified it will also save the target phase cube data (full path on the server)

            ncpa:  (int) : !!experimental!!!: Used only in the context of PYRWFS + NCPA compensation on the fly (with optical gain)
            defines how many iters the NCPA refslopes are updates with the proper optical gain. Ex: if NCPA=10 refslopes will be updates every 10 iters.

            ncpa_wfs:  (int) : the ncpa phase as seen from the wfs array with dims = size of Mpupil

            ref_slopes:  (int) : the reference slopes to use.

            ditch_strehl:  (int) : resets the long exposure SR computation at the beginning of the Circular buffer (default= True)

        Return:
            slopes:  (int) : the slopes CB

            volts:  (int) : the volts applied to the DM(s) CB

            ai:  (int) : the modal coefficient of the residual phase projected on the currently used modal Basis

            psf_le:  (int) : Long exposure PSF over the <cb_count> iterations (I.e SR is reset at the begining of the CB if ditch_strehl=True)

            sthrel_se_list:  (int) : The SR short exposure evolution during CB recording

            sthrel_le_list:  (int) : The SR long exposure evolution during CB recording

            g_ncpa_list:  (int) : the gain applied to the NCPA (PYRWFS CASE) if NCPA is set to True

            cube_data:  (int) : the tarPhase or psfse cube data (see cube_data_type)
        """
        slopes_data = None
        volts_data = None
        cube_data = None
        ai_data = None
        k = 0
        sthrel_se_list = []
        sthrel_le_list = []
        g_ncpa_list = []

        # Resets the target so that the PSF LE is synchro with the data
        # Doesn't reset it if Ditch_strehl == False (used for real time gain computation)
        if ditch_strehl:
            for i in range(len(self._sim.config.p_targets)):
                self.reset_strehl(i)

        # Starting CB loop...
        for j in range(cb_count):
            print(j, end="\r")
            if (ncpa):
                if (j % ncpa == 0):
                    ncpa_diff = ref_slopes[None, :]
                    ncpa_turbu = self.do_imat_phase(controller_index,
                                                    -ncpa_wfs[None, :, :], noise=False,
                                                    with_turbu=True)
                    g_ncpa = float(
                            np.sqrt(
                                    np.dot(ncpa_diff, ncpa_diff.T) / np.dot(
                                            ncpa_turbu, ncpa_turbu.T)))
                    if (g_ncpa > 1e18):
                        g_ncpa = 0
                        print('Warning NCPA ref slopes gain too high!')
                        g_ncpa_list.append(g_ncpa)
                        self.rtc.set_ref_slopes(-ref_slopes * g_ncpa)
                    else:
                        g_ncpa_list.append(g_ncpa)
                        print('NCPA ref slopes gain: %4.3f' % g_ncpa)
                        self.rtc.set_ref_slopes(-ref_slopes / g_ncpa)

            self._sim.next(see_atmos=see_atmos)
            for t in range(len(self._sim.config.p_targets)):
                self._sim.comp_tar_image(t)

            srse, srle, _, _ = self.get_strehl(tar_index)
            sthrel_se_list.append(srse)
            sthrel_le_list.append(srle)
            if (j % sub_sample == 0):
                ai_vector = self.compute_modal_residuals()
                if (ai_data is None):
                    ai_data = np.zeros((len(ai_vector), int(cb_count / sub_sample)))
                ai_data[:, k] = ai_vector

                slopes_vector = self.rtc.get_slopes(0)
                if (slopes_data is None):
                    slopes_data = np.zeros((len(slopes_vector),
                                            int(cb_count / sub_sample)))
                slopes_data[:, k] = slopes_vector

                volts_vector = self.get_command(0)
                if (volts_data is None):
                    volts_data = np.zeros((len(volts_vector),
                                           int(cb_count / sub_sample)))
                volts_data[:, k] = volts_vector

                if (cube_data_type):
                    if (cube_data_type == "tarPhase"):
                        dataArray = self.get_tar_phase(tar_index, pupil=True)
                    elif (cube_data_type == "psfse"):
                        dataArray = self.get_tar_image(tar_index, "se")
                    else:
                        raise ValueError("unknown dataData" % cube_data_type)
                    if (cube_data is None):
                        cube_data = np.zeros((*dataArray.shape,
                                              int(cb_count / sub_sample)))
                    cube_data[:, :, k] = dataArray
                k += 1
        if (cube_data_file_path != ""):
            print("Saving tarPhase cube at: ", cube_data_file_path)
            pfits.writeto(cube_data_file_path, cube_data, overwrite=True)

        psf_le = self.get_tar_image(tar_index, "le")
        return slopes_data, volts_data, ai_data, psf_le, sthrel_se_list, sthrel_le_list, g_ncpa_list, cube_data
