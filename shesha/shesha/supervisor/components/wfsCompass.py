## @package   shesha.supervisor
## @brief     User layer for initialization and execution of a COMPASS simulation
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
from shesha.init.wfs_init import wfs_init
from shesha.supervisor.components.sourceCompass import SourceCompass
import numpy as np

class WfsCompass(SourceCompass):
    """ WFS handler for compass simulation

    Attributes:
        wfs : (sutraWrap.Wfs) : SutraSensors instance

        context : (carmaContext) : CarmaContext instance

        config : (config module) : Parameters configuration structure module

        sources : (List) : List of SutraSource instances used for raytracing
    """
    def __init__(self, context, config):
        """ Initialize a wfsCompass component for wfs related supervision

        Parameters:
            context : (carmaContext) : CarmaContext instance

            config : (config module) : Parameters configuration structure module
        """
        self.context = context
        self.config = config # Parameters configuration coming from supervisor init
        self.wfs = wfs_init(self.context, self.tel, self.config.p_wfss,
                                self.config.p_tel, self.config.p_geom, self.config.p_dms,
                                self.config.p_atmos)
        self.sources = [wfs.d_gs for wfs in self.wfs.d_wfss]
        
    def get_wfs_image(self, wfs_index : int) -> np.ndarray:
        """ Get an image from the WFS (wfs[0] by default), or from the centroider handling the WFS
        to get the calibrated image

        Parameters:
            wfs_index : (int) : index of the WFS (or the centroider) to request an image

        Return:
            image : (np.ndarray) : WFS image
        """
        return np.array(self.wfs.d_wfs[wfs_index].d_binimg)

    def set_pyr_modulation_points(self, wfs_index : int, cx: np.ndarray, cy: np.ndarray,
                                  weights: np.ndarray = None) -> None:
        """ Set pyramid modulation positions

        Parameters:
            wfs_index : (int) : WFS index

            cx : (np.ndarray) : X positions of the modulation points [arcsec]

            cy : (np.ndarray) : Y positions of the modulation points [arcsec]

            weights : (np.ndarray, optional) : Weights to apply on each modulation point contribution
        """
        pyr_npts = len(cx)
        pwfs = self.config.p_wfss[wfs_index]
        pwfs.set_pyr_npts(pyr_npts)
        pwfs.set_pyr_cx(cx)
        pwfs.set_pyr_cy(cy)
        if weights is not None:
            self.wfs.d_wfs[wfs_index].set_pyr_modulation_points(cx, cy, pyr_npts)
        else:
            self.wfs.d_wfs[wfs_index].set_pyr_modulation_points(
                    cx, cy, weights, pyr_npts)

    def set_fourier_mask(self, wfs_index : int, new_mask: np.ndarray) -> None:
        """ Set a mask in the Fourier Plane of the given WFS

        Parameters:
            wfs_index : (int, optional) : WFS index

            new_mask : (ndarray) : mask to set
        """
        if new_mask.shape != self.config.p_wfss[wfs_index].get_halfxy().shape:
            print('Error : mask shape should be {}'.format(
                    self.config.p_wfss[wfs_index].get_halfxy().shape))
        else:
            self.wfs.d_wfs[wfs_index].set_phalfxy(
                    np.exp(1j * np.fft.fftshift(new_mask)).astype(np.complex64).T)

    def set_noise(self, wfs_index : int, noise: float, seed: int = 1234) -> None:
        """ Set noise value of WFS wfs_index

        Parameters:
            wfs_index : (int, optional) : WFS index

            noise : (float) : readout noise value in e-

            seed : (int, optional) : RNG seed. The seed used will be computed as seed + wfs_index
                                     Default is 1234
        """
        self.wfs.d_wfs[wfs_index].set_noise(noise, int(seed + wfs_index))
        print("Noise set to: %f on WFS %d" % (noise, wfs_index))

    def set_gs_mag(self, wfs_index : int, mag : float) -> None:
        """ Change the guide star magnitude for the given WFS

        Parameters:
            wfs_index : (int, optional) : WFS index

            mag : (float) : New magnitude of the guide star
        """
        wfs = self.wfs.d_wfs[wfs_index]
        if (self.config.p_wfs0.type == "pyrhr"):
            r = wfs.comp_nphot(self.config.p_loop.ittime,
                               self.config.p_wfss[wfs_index].optthroughput,
                               self.config.p_tel.diam, self.config.p_tel.cobs,
                               self.config.p_wfss[wfs_index].zerop, mag)
        else:
            r = wfs.comp_nphot(self.config.p_loop.ittime,
                               self.config.p_wfss[wfs_index].optthroughput,
                               self.config.p_tel.diam, self.config.p_wfss[wfs_index].nxsub,
                               self.config.p_wfss[wfs_index].zerop, mag)
        if (r == 0):
            print("GS magnitude is now %f on WFS %d" % (mag, wfs_index))

    def compute_wfs_image(self, wfs_index : int, noise: bool = True) -> None:
        """ Computes the image produced by the WFS from its phase screen

        Parameters :
            wfs_index : (int): WFS index

            noise : (bool, optional) : Flag to enable noise for image computation. Default is True
        """
        self.wfs.d_wfs[wfs_index].comp_image(noise)

    def reset_noise(self) -> None:
        """ Reset all the WFS RNG to their original state
        """
        for wfs_index in range(len(self.config.p_wfss)):
            self.wfs.d_wfs[wfs_index].set_noise(
                    [p.noise for p in self.config.pwfss], 1234 + wfs_index)

    def get_ncpa_wfs(self, wfs_index : int) -> np.ndarray:
        """ Return the current NCPA phase screen of the WFS path

        Parameters:
            wfs_index : (int) : Index of the WFS

        Return:
            ncpa : (np.ndarray) : NCPA phase screen
        """
        return np.array(self.wfs.d_wfs[wfs_index].d_gs.d_ncpa_phase)

    def get_wfs_phase(self, wfs_index : int) -> np.ndarray:
        """ Return the WFS phase screen of WFS number wfs_index

        Parameters:
            wfs_index : (int) : Index of the WFS

        Return:
            phase : (np.ndarray) : WFS phase screen
        """
        return np.array(self.wfs.d_wfs[wfs_index].d_gs.d_phase)

    def get_pyrhr_image(self, wfs_index : int) -> np.ndarray:
        """ Get an high resolution image from the PWFS

        Parameters:
            wfs_index : (int) : Index of the WFS

        Return:
            image : (np.ndarray) : PWFS high resolution image

        """
        return np.array(self.wfs.d_wfs[wfs_index].d_hrimg)

    def set_ncpa_wfs(self, wfs_index : int, ncpa: np.ndarray) -> None:
        """ Set the additional fixed NCPA phase in the WFS path.
        ncpa must be of the same size of the mpupil support

        Parameters:
            wfs_index : (int) : WFS index

            ncpa : (ndarray) : NCPA phase screen to set [µm]
        """
        self.wfs.d_wfs[wfs_index].d_gs.set_ncpa(ncpa)

    def set_wfs_phase(self, wfs_index : int, phase : np.ndarray) -> None:
        """ Set the phase screen seen by the WFS

        Parameters:
            wfs_index : (int) : WFS index

            phase : (np.ndarray) : phase screen to set
        """
        self.wfs.d_wfs[wfs_index].d_gs.set_phase(phase)

    def set_wfs_pupil(self, wfs_index : int, pupil : np.ndarray) -> None:
        """ Set the pupil seen by the WFS
        Other pupils remain unchanged, i.e. DM and target can see an other
        pupil than the WFS after this call.
        <pupil> must have the same shape than p_geom._mpupil support

        Parameters:
            wfs_index : (int) : WFS index

            pupil : (np.ndarray) : new pupil to set
        """
        old_mpup = self.config.p_geom._mpupil
        dimx = old_mpup.shape[0]
        dimy = old_mpup.shape[1]
        if ((mpupil.shape[0] != dimx) or (mpupil.shape[1] != dimy)):
            print("Error mpupil shape on wfs %d must be: (%d,%d)" % (wfs_index, dimx,
                                                                     dimy))
        else:
            self.wfs.d_wfs[wfs_index].set_pupil(mpupil.copy())

    def get_pyr_focal_plane(self, wfs_index : int) -> np.ndarray:
        """ Returns the psf on the top of the pyramid.
        pyrhr WFS only

        Parameters:
            wfs_index : (int) : WFS index

        Return:
            focal_plane : (np.ndarray) : psf on the top of the pyramid
        """
        return np.fft.fftshift(np.array(self.wfs.d_wfs[wfs_index].d_pyrfocalplane))
