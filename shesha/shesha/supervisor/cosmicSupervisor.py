#
# This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
#
# COMPASS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COMPASS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COMPASS. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright (C) 2011-2024 COSMIC Team


from shesha.supervisor.compassSupervisor import CompassSupervisor
from shesha.supervisor.components import RtcCosmic
import time


from typing import Iterable


class CosmicSupervisor(CompassSupervisor):
    """This class implements generic supervisor to handle compass simulation

    Attributes inherited from GenericSupervisor:
        context : (CarmaContext) : a CarmaContext instance

        config : (config) : Parameters structure

        is_init : (bool) : Flag equals to True if the supervisor has already been initialized

        iter : (int) : Frame counter

    Attributes:
        telescope : (TelescopeComponent) : a TelescopeComponent instance

        atmos : (AtmosComponent) : An AtmosComponent instance

        target : (targetComponent) : A TargetComponent instance

        wfs : (WfsComponent) : A WfsComponent instance

        dms : (DmComponent) : A DmComponent instance

        rtc : (RtcComponent) : A Rtc component instance

        basis : (ModalBasis) : a ModalBasis instance (optimizer)

        calibration : (Calibration) : a Calibration instance (optimizer)

        modalgains : (ModalGains) : a ModalGain instance (optimizer) using CLOSE algorithm

        close_modal_gains : (list of floats) : list of the previous values of the modal gains
    """

    #     ___                  _      __  __     _   _            _
    #    / __|___ _ _  ___ _ _(_)__  |  \/  |___| |_| |_  ___  __| |___
    #   | (_ / -_) ' \/ -_) '_| / _| | |\/| / -_)  _| ' \/ _ \/ _` (_-<
    #    \___\___|_||_\___|_| |_\__| |_|  |_\___|\__|_||_\___/\__,_/__/

    def _init_rtc(self):
        """Initialize the rtc component of the supervisor as a RtcCompass"""
        super()._init_rtc()
        self.hrtc = RtcCosmic(self.config, self.wfs, self.dms)

    # def _init_components(self) -> None:
    #     """ Initialize all the components
    #     """
    #     super()._init_components()

    def next(
        self,
        *,
        move_atmos: bool = True,
        nControl: int = 0,
        tar_trace: Iterable[int] = None,
        wfs_trace: Iterable[int] = None,
        do_control: bool = True,
        apply_control: bool = True,
        compute_tar_psf: bool = True,
        compute_corono: bool = True,
    ) -> None:
        """Iterates the AO loop, with optional parameters.

        Overload the GenericSupervisor next() method to handle the GEO controller
        specific raytrace order operations

        Kwargs:
            move_atmos: (bool): move the atmosphere for this iteration. Default is True

            nControl: (int): Controller number to use. Default is 0 (single control configuration)

            tar_trace: (List): list of targets to trace. None is equivalent to all (default)

            wfs_trace: (List): list of WFS to trace. None is equivalent to all (default)

            do_control : (bool) : Performs RTC operations if True (Default)

            apply_control: (bool): if True (default), apply control on DMs

            compute_tar_psf : (bool) : If True (default), computes the PSF at the end of the iteration

            compute_corono: (bool): If True (default), computes the coronagraphic image
        """

        if tar_trace is None and self.target is not None:
            tar_trace = range(len(self.config.p_targets))
        if wfs_trace is None and self.wfs is not None:
            wfs_trace = range(len(self.config.p_wfss))

        if move_atmos and self.atmos is not None:
            self.atmos.move_atmos()
        # in case there is at least 1 controller GEO in the controller list : use this one only
        self.tel.update_input_phase()
        if tar_trace is not None:  # already checked at line 213?
            for t in tar_trace:
                if self.atmos.is_enable:
                    self.target.raytrace(t, tel=self.tel, atm=self.atmos, dms=self.dms)
                else:
                    self.target.raytrace(t, tel=self.tel, dms=self.dms)

        if wfs_trace is not None:  # already checked at line 215?
            for w in wfs_trace:
                if self.atmos.is_enable:
                    self.wfs.raytrace(w, tel=self.tel, atm=self.atmos)
                else:
                    self.wfs.raytrace(w, tel=self.tel)

                if not self.config.p_wfss[w].open_loop and self.dms is not None:
                    self.wfs.raytrace(w, dms=self.dms, ncpa=False, reset=False)
                self.wfs.compute_wfs_image(w)

        if do_control and self.hrtc is not None:
            self.hrtc.do_control()
        if apply_control:
            self.hrtc.apply_control()

        if compute_tar_psf:
            for tar_index in tar_trace:
                self.target.comp_tar_image(tar_index)
                self.target.comp_strehl(tar_index)

        if self.corono is not None and compute_corono:
            for coro_index in range(len(self.config.p_coronos)):
                self.corono.compute_image(coro_index)

        self.iter += 1

    def _print_strehl(
        self,
        monitoring_freq: int,
        iters_time: float,
        total_iters: int,
        *,
        tar_index: int = 0,
    ):
        """Print the Strehl ratio SE and LE from a target on the terminal, the estimated remaining time and framerate

        Args:
            monitoring_freq : (int) : Number of frames between two prints

            iters_time : (float) : time elapsed between two prints

            total_iters : (int) : Total number of iterations

        Kwargs:
            tar_index : (int) : Index of the target. Default is 0
        """
        framerate = monitoring_freq / iters_time
        strehl = self.target.get_strehl(tar_index)
        etr = (total_iters - self.iter) / framerate
        print(
            "%d \t %.3f \t  %.3f\t     %.1f \t %.1f"
            % (self.iter + 1, strehl[0], strehl[1], etr, framerate)
        )

    def loop(
        self,
        number_of_iter: int,
        *,
        monitoring_freq: int = 100,
        compute_tar_psf: bool = True,
        **kwargs,
    ):
        """Perform the AO loop for <number_of_iter> iterations

        Args:
            number_of_iter: (int) : Number of iteration that will be done

        Kwargs:
            monitoring_freq: (int) : Monitoring frequency [frames]. Default is 100

            compute_tar_psf : (bool) : If True (default), computes the PSF at each iteration
                                                 Else, only computes it each <monitoring_freq> frames
        """
        if not compute_tar_psf:
            print("WARNING: Target PSF will be computed (& accumulated) only during monitoring")

        print("----------------------------------------------------")
        print("iter# | S.E. SR | L.E. SR | ETR (s) | Framerate (Hz)")
        print("----------------------------------------------------")
        # self.next(**kwargs)
        t0 = time.time()
        t1 = time.time()
        if number_of_iter == -1:  # Infinite loop
            while True:
                self.next(compute_tar_psf=compute_tar_psf, **kwargs)
                if (self.iter + 1) % monitoring_freq == 0:
                    if not compute_tar_psf:
                        self.target.comp_tar_image(0)
                        self.target.comp_strehl(0)
                    self._print_strehl(monitoring_freq, time.time() - t1, number_of_iter)
                    t1 = time.time()

        for _ in range(number_of_iter):
            self.next(compute_tar_psf=compute_tar_psf, **kwargs)
            if (self.iter + 1) % monitoring_freq == 0:
                if not compute_tar_psf:
                    self.target.comp_tar_image(0)
                    self.target.comp_strehl(0)
                self._print_strehl(monitoring_freq, time.time() - t1, number_of_iter)
                t1 = time.time()
        t1 = time.time()
        print(
            " loop execution time:",
            t1 - t0,
            "  (",
            number_of_iter,
            "iterations), ",
            (t1 - t0) / number_of_iter,
            "(mean)  ",
            number_of_iter / (t1 - t0),
            "Hz",
        )

    def reset(self):
        """Reset the simulation to return to its original state"""
        self.atmos.reset_turbu()
        self.wfs.reset_noise()
        for tar_index in range(len(self.config.p_targets)):
            self.target.reset_strehl(tar_index)
        self.dms.reset_dm()
