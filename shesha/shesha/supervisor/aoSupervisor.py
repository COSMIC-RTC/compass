## @package   shesha.supervisor.aoSupervisor
## @brief     Abstract layer for initialization and execution of a AO supervisor
## @author    COMPASS Team <https://github.com/ANR-COMPASS>
## @version   4.4.1
## @date      2011/01/28
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

from shesha.supervisor.abstractSupervisor import AbstractSupervisor
from shesha.constants import CentroiderType, CONST, WFSType
import numpy as np
from tqdm import trange
import os
from collections import OrderedDict
import astropy.io.fits as pfits
import shesha.constants as scons

class AoSupervisor(AbstractSupervisor):

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

    def __init__(self):
        self.use_close = False
        self.phase_to_modes = None
        self.P = None
        self.current_buffer = 1


    def get_config(self):
        """ Returns the configuration in use, in a supervisor specific format ? 

        Return:
            config : (config module) : Current supervisor configuration
        """
        return self.config

    def set_flat(self, flat: np.ndarray, centro_index: int = 0):
        """ Load flat field for the given wfs

        Parameters:
            flat : (np.ndarray) : New WFS flat to use

            centro_index : (int, optional) : index of the centroider handling the WFS
        """
        self.rtc.d_centro[centro_index].set_flat(flat, flat.shape[0])

    def set_dark(self, background: np.ndarray, centro_index: int = 0):
        """ Load background for the given wfs

        Parameters:
            dark : (np.ndarray) : New WFS dark to use

            centro_index : (int, optional) : index of the centroider handling the WFS
        """
        self.rtc.d_centro[centro_index].set_dark(background, background.shape[0])

    def compute_slopes(self, controller_index: int = 0):
        """ Compute the slopes handled by a controller, and returns it

        Parameters:
            controller_index : (int, optional) : Controller index that will compute its slopes. Default is 0
        
        Return:
            slopes : (np.ndarray) : Slopes vector
        """
        self.rtc.do_centroids(controller_index)
        return self.get_slopes(controller_index)

    def get_wfs_image(self, wfs_index: int = 0, cal_pix=False) -> np.ndarray:
        """ Get an image from the WFS (wfs[0] by default), or from the centroider handling the WFS
        to get the calibrated image

        Parameters:
            wfs_index : (int) : index of the WFS (or the centroider) to request an image

            cal_pix : (bool, optional) : Flag to get calibrated image instead of raw one. Default is False
        
        Return:
            image : (np.ndarray) : WFS image
        """
        if (cal_pix):
            if self.rtc.d_centro[numWFS].d_dark is None: # Simulation case
                return np.array(self._sim.wfs.d_wfs[numWFS].d_binimg)
            if self.rtc.d_centro[numWFS].type == CentroiderType.MASKEDPIX: # Standalone case
                #     self.rtc.d_centro[numWFS].fill_selected_pix(
                #             self.rtc.d_control[0].d_centroids)
                #     return np.array(self.rtc.d_centro[numWFS].d_selected_pix)
                # else:
                self.rtc.d_centro[numWFS].fill_selected_pix(
                        self.rtc.d_control[0].d_centroids)
                mask = np.array(self.rtc.d_centro[numWFS].d_selected_pix)
                mask[np.where(mask)] = np.array(
                        self.rtc.d_centro[numWFS].d_img)[np.where(mask)]
                return mask
            return np.array(self.rtc.d_centro[numWFS].d_img)
        else: 
            if self.rtc.d_centro[numWFS].d_dark is None: # Simulation case
                return np.array(self._sim.wfs.d_wfs[numWFS].d_binimg)
            return np.array(self.rtc.d_centro[numWFS].d_img_raw) # Standalone case

    def get_intensities(self) -> np.ndarray:
        """ Return sum of intensities in subaps. Size nSubaps, same order as slopes
        """
        raise NotImplementedError("Not implemented")
        # return np.empty(1)

    def set_perturbation_voltage(self, controller_index: int, name: str,
                               command: np.ndarray) -> None:
        """ Add circular buffer of offset values to integrator (will be applied at the end of next iteration)

        Parameters:
            nControl : (int) : Controller index

            name : (str) : Buffer name

            command : (np.ndarray) : perturbation voltage circular buffer
        """
        if len(command.shape) == 1:
            self.rtc.d_control[controller_index].set_perturb_voltage(name, command, 1)
        elif len(command.shape) == 2:
            self.rtc.d_control[controller_index].set_perturb_voltage(name, command,
                                                             command.shape[0])
        else:
            raise AttributeError("command should be a 1D or 2D array")

    def reset_perturbation_voltage(self, controller_index: int) -> None:
        """ Reset the perturbation voltage of the controller_index controller
        (i.e. will remove ALL perturbation voltages.)
        If you want to reset just one, see the function removePerturbationVoltage()

        Parameters:
            controller_index : (int) : controller index from where to remove the buffer
        """
        self.rtc.d_control[controller_index].reset_perturb_voltage()

    def removePerturbationVoltage(self, controller_index: int, name: str) -> None:
        """ Remove the perturbation voltage called <name>, from the controller number <controller_index>.
        If you want to remove all of them, see function reset_perturbation_voltage()

        Parameters:
            controller_index : (int) : controller index from where to remove the buffer
        """
        self.rtc.d_control[controller_index].remove_perturb_voltage(name)

    def get_slopes(self, controller_index: int = 0):
        """ Return the current slopes vector of the controller_index controller

        Parameters:
            controller_index : (int) : controller index handling the slopes

        Return:
            slopes : (np.ndarray) : Current slopes vector containing slopes of all
                                    the WFS handled by the specified controller
        """
        return np.array(self.rtc.d_control[controller_index].d_centroids)

    def get_err(self, controller_index: int = 0) -> np.ndarray:
        """ Get integrator increment from controller_index controller

        Parameters:
            controller_index : (int) : controller index           
        """
        return np.array(self.rtc.d_control[controller_index].d_err)

    def get_voltages(self, controller_index: int = 0) -> np.ndarray:
        """ Get voltages vector (i.e. vector sent to the DM) from controller_index controller

        Parameters:
            controller_index : (int) : controller index  

        Return:
            voltages : (np.ndarray) : current voltages vector         

        """
        return np.array(self.rtc.d_control[controller_index].d_voltage)

    def set_integrator_law(self, controller_index: int = 0):
        """ Set the control law to integrator (controller generic only)
            v[k] = v[k-1] + g.R.s[k]

        Parameters:
            controller_index: (int): controller index
        """
        self.rtc.d_control[controller_index].set_commandlaw("integrator")

    def set_2matrices_law(self, controller_index: int = 0):
        """ Set the control law to 2matrices (controller generic only)
        v[k] = decayFactor.E.v[k-1] + g.R.s[k]

        Parameters:
            controller_index: (int): controller index
        """
        self.rtc.d_control[controller_index].set_commandlaw("2matrices")

    def set_modal_integrator_law(self, controller_index: int = 0):
        """ Set the control law to 2matrices (controller generic only)
        v[k] = v[k-1] + E.g.R.s[k]

        Parameters:
            controller_index: (int): controller index
        """
        self.rtc.d_control[controller_index].set_commandlaw("modal_integrator")

    def set_decay_factor(self, decay, controller_index: int = 0):
        """ Set the decay factor used in 2matrices command law (controller generic only)

        Parameters :
            decay : (np.ndarray) : decay factor vector

            controller_index: (int): controller index
        """
        self.rtc.d_control[controller_index].set_decayFactor(decay)

    def set_E_matrix(self, e_matrix, controller_index: int = 0):
        """ Set the E matrix used in 2matrices or modal command law (controller generic only)

        Parameters :
            e_matrix : (np.ndarray) : E matrix to set

            controller_index: (int): controller index
        """
        self.rtc.d_control[controller_index].set_matE(e_matrix)

    def do_ref_slopes(self, controller_index: int = 0):
        """ Computes and set a new reference slopes for each WFS handled by
        the specified controller

        Parameters:
            controller_index: (int, optional): controller index. Default is 0
        """
        print("Doing reference slopes...")
        self.rtc.do_centroids_ref(controller_index)
        print("Reference slopes done")

    def reset_ref_slopes(self, controller_index: int = 0):
        """ Reset the reference slopes of each WFS handled by the specified controller

        Parameters:
            controller_index: (int, optional): controller index. Default is 0            
        """
        for centro in self.rtc.d_centro:
            centro.d_centroids_ref.reset()

    def close_loop(self, controller_index: int = 0) -> None:
        """ DM receives controller output + pertuVoltage

        Parameters:
            controller_index: (int, optional): controller index. Default is 0     
        """
        self.rtc.d_control[controller_index].set_open_loop(0)  # close_loop

    def open_loop(self, rst=True, controller_index: int = 0) -> None:
        """ Integrator computation goes to /dev/null but pertuVoltage still applied

        Parameters:
            rst : (bool, optional) : If True (default), integrator is reset

            controller_index: (int, optional): controller index. Default is 0    
        """
        self.rtc.d_control[controller_index].set_open_loop(1, rst)  # open_loop

    def set_ref_slopes(self, ref_slopes: np.ndarray, centro_index=None) -> None:
        """ Set given ref slopes in controller

        Parameters:
            ref_slopes : (ndarray) : Reference slopes vectoronly set the reference slop

            centro_index : (int, optionnal) : If given, only set the reference slopes vector
                                             used by the specified centroider. If None, the reference
                                             slopes vector must be a concatenation of all the reference
                                             slopes to use for each centroiders handled by the controller
        """
        if (centro_index is None):
            self.rtc.set_centroids_ref(ref_slopes)
        else:
            self.rtc.d_centro[centro_index].set_centroids_ref(ref_slopes)

    def get_ref_slopes(self, centro_index=None) -> np.ndarray:
        """ Get the currently used reference slopes

        Parameters:
            centro_index : (int, optionnal) : If given, only get the reference slopes vector
                                             used by the specified centroider. If None, the reference
                                             slopes vector returned is a concatenation of all the reference
                                             slopes used for by centroiders in the RTC
        
        Return:
            ref_slopes : (np.ndarray) : Reference slopes vector
        """
        ref_slopes = np.empty(0)
        if (centro_index is None):
            for centro in self.rtc.d_centro:
                ref_slopes = np.append(ref_slopes, np.array(centro.d_centroids_ref))
            return ref_slopes
        else:
            return np.array(self.rtc.d_centro[centro_index].d_centroids_ref)

    def get_interaction_matrix(self, controller_index: int = 0):
        """ Return the interaction matrix of the controller

        Parameters :
            controller_index: (int): controller index

        Return:
            imat : (np.ndarray) : Interaction matrix currently set in the controller
        """
        return np.array(self.rtc.d_control[controller_index].d_imat)

    def get_command_matrix(self, controller_index: int = 0):
        """ Return the command matrix of the controller

        Parameters:
            controller_index: (int): controller index

        Return:
            cmat : (np.ndarray) : Command matrix currently used by the controller
        """
        return np.array(self.rtc.d_control[controller_index].d_cmat)

    def set_centroider_threshold(self, centro_index: int = 0, thresh: float = 0.):
        """ Set the threshold value of a thresholded COG

        Parameters:
            centro_index: (int): centroider index

            thresh: (float): new threshold value
        """
        self.rtc.d_centro[centro_index].set_threshold(thresh)

    def get_pyr_method(self, centro_index):
        """ Get pyramid compute method currently used

        Parameters:
            centro_index: (int): centroider index

        Return:
            method : (str) : Pyramid compute method currently used
        """
        return self.rtc.d_centro[centro_index].pyr_method

    def set_pyr_method(self, pyr_method, centro_index: int = 0):
        """ Set the pyramid method for slopes computation

        Parameters :
            pyr_method : (int) : new centroiding method (0: nosinus global
                                                1: sinus global
                                                2: nosinus local
                                                3: sinus local)
        """
        self.rtc.d_centro[centro_index].set_pyr_method(pyr_method)  # Sets the pyr method
        self.rtc.do_centroids(0)  # To be ready for the next get_slopess
        print("PYR method set to " + self.rtc.d_centro[centro_index].pyr_method)

    def set_gain(self, controller_index : int, gain : float) -> None:
        """ Set the scalar gain

        Parameters:
            controller_index : (int) : Index of the controller to modify

            gain : (float) : scalar gain of modal gain to set
        """
        self.rtc.d_control[controller_index].set_gain(gain)

    def set_modal_gains(self, mgain : np.ndarray, controller_index : int=0):
        """ Sets the modal gain (when using modal integrator control law)

        Parameters:
            mgain : (np.ndarray) : Modal gains to set

            controller_index : (int, optional) : Controller index to modify. Default is 0
        """
        self.rtc.d_control[control].set_modal_gains(mgain)

    def get_modal_gains(self, controller_index : int=0):
        """ Returns the modal gains (when using modal integrator control law)  

        Parameters:
            controller_index : (int, optional) : Controller index to modify. Default is 0

        Return:
            mgain : (np.ndarray) : Modal gains vector currently used
        """
        return np.array(self.rtc.d_control[control].d_gain)

    def set_command_matrix(self, cmat: np.ndarray, controller_index : int=0) -> None:
        """ Set the command matrix for the controller to use

        Parameters:
            cmat : (np.ndarray) : command matrix to set
            controller_index : (int, optional) : Controller index to modify. Default is 0
        """
        self.rtc.d_control[controller_index].set_cmat(cmat)

    def get_frame_counter(self) -> int:
        """Return the current iteration number of the loop

        Return:
            framecounter : (int) : Number of iteration already performed
        """
        if not self.is_init:
            print('Warning - requesting frame counter of uninitialized BenchSupervisor.')
        return self.iter

    def get_masked_pix(self, centro_index: int = 0):
        """ Return the mask of valid pixels used by a maskedpix centroider

        Parameters:
            centro_index : (int): Centroider index. Must be a maskedpix centroider
        
        Return:
            mask : (np.ndarray) : Mask used
        """
        if (self.rtc.d_centro[centro_index].type != CentroiderType.MASKEDPIX):
            raise TypeError("Centroider must be a maskedpix one")
        self.rtc.d_centro[centro_index].fill_mask()
        return np.array(self.rtc.d_centro[centro_index].d_mask)